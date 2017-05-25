# -*- coding: utf-8 -*-
"""
CIF 1.0, 1.1 and 2.0 parser based on PyCifRW.

@author : Laurent P. Rene de Cotret 
"""
import numpy as np
from CifFile import CifFile, get_number_with_esd
import spglib

from . import Atom, Lattice, lattice_vectors_from_parameters, real_coords

# Useful relations
# Generated from spglib instead of cctbx, but should be equivalent 
# to the data in cif2cell
Hall2HM = dict( (dataset['hall_symbol'],dataset['international_short']) 
				for dataset in map(spglib.get_spacegroup_type, range(1, 531)))
HM2Hall = {v:k for (k,v) in Hall2HM.items()}

Hall2Number = dict( (dataset['hall_symbol'],dataset['number']) 
				for dataset in map(spglib.get_spacegroup_type, range(1, 531)))
Number2Hall = {v:k for (k,v) in Hall2Number.items()}

SymOpsHall = dict( map(lambda i: (spglib.get_spacegroup_type(i)['hall_symbol'], 
				   				  spglib.get_symmetry_from_database(i)), range(1, 531)) )


class CIFParser(object):
	"""
	Collection of methods that parses CIF files based on PyCifRW.
    
	Attributes
	----------
	file : str
		Absolute path to the CIF file associated with the parser's target structure.

	Methods
	-------    
	lattice_vectors
		Lattice vectors of the crystal structure.

	atoms
		Returns a list of atoms making up the asymmetric unit cell

	symmetry_operators
		Symmetry operators that relate the asymmetric unit cell and the unit cell.
	"""
	def __init__(self, filename):
		self.file = CifFile(datasource = filename)

	def __enter__(self):
		return self
	
	def __exit__(self, type, value, traceback):
		del self.file
	
	def block_containing(self, key):
		""" Returns the CifFile block containing a kay """
		block_name = filter(lambda block: key in self.file[block], self.file.keys())
		try:
			block = self.file[next(block_name)]
		except StopIteration:
			raise KeyError
		return block

	def lattice_vectors(self):
		""" 
		Returns the lattice vectors associated to a CIF structure.
        
		Returns
		-------
		lv : list of ndarrays, shape (3,)
		"""
		block = self.block_containing('_cell_length_a')
		
		a, _ = get_number_with_esd(block["_cell_length_a"])
		b, _ = get_number_with_esd(block["_cell_length_b"])
		c, _ = get_number_with_esd(block["_cell_length_c"])
		alpha, _ = get_number_with_esd(block["_cell_angle_alpha"])
		beta, _ = get_number_with_esd(block["_cell_angle_beta"])
		gamma, _ = get_number_with_esd(block["_cell_angle_gamma"])

		return lattice_vectors_from_parameters(a, b, c, alpha, beta, gamma)
	
	def symmetry_operators(self):
		"""
		Returns the symmetry operators that map the atomic positions in a
		CIF file to the crystal unit cell.

		Returns
		-------
		sym_ops : list of ndarrays, shape (4,4)
			Transformation matrices. Since translations and rotation are combined,
			the transformation matrices are 4x4.
		"""
		# TODO: parse '_space_group_symop_operation_xyz' key
		return [np.eye(4)]
	
	def atoms(self):
		"""
		Asymmetric unit cell.

		Returns
		-------
		atoms : iterable of Atom
		"""
		block = self.block_containing('_atom_site_label')
		elements = map(lambda label: label[0:2].lower().title(), block['_atom_site_label'])
		xs, ys, zs = block['_atom_site_fract_x'], block['_atom_site_fract_y'], block['_atom_site_fract_y']
		coords = [np.array([get_number_with_esd(x)[0], get_number_with_esd(y)[0], get_number_with_esd(z)[0]]) 
				  for x,y,z in zip(xs, ys, zs)]

		lv = self.lattice_vectors()
		return [Atom(element = element, coords = real_coords(coord, lv))
				for element, coord in zip(elements, coords)]
