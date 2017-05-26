# -*- coding: utf-8 -*-
"""
CIF 1.0, 1.1 and 2.0 parser based on cif2cell.

References
----------
..[1] Torbjorn Bjorkman, "CIF2Cell: Generating geometries for electronic structure programs", 
Computer Physics Communications 182, 1183-1186 (2011)
doi: 10.1016/j.cpc.2011.01.013
"""
import string
import warnings
from contextlib import suppress
from functools import lru_cache

import numpy as np
from CifFile import CifFile, get_number_with_esd

from . import Atom, Lattice, lattice_vectors_from_parameters, real_coords
from .. import affine_map, transform
from .spg_data import HM2Hall, Number2Hall, SymOpsHall

class ParseError(Exception):
	pass

def sym_ops(equiv_site):
	""" Parse a symmetry operator from an equivalent-site representation 
	
	Parameters
	----------
	equiv_site : str or iterable of strings
		Either comma-separated string e.g. "+y, +x, -z + 1/2" or an
		iterable of the comma-separated values, e.g. ["+y", "+x", "-z + 1/2"] 
	
	Returns
	-------
	sym_ops : ndarray, shape (4,4)
		Symmetry operator as a 4x4 affine transformation matrix.
	"""
	rotation = np.zeros( (3,3) )
	trans_vec = np.zeros( (3,) )

	if isinstance(equiv_site, str):
		equiv_site = equiv_site.split(',')
		
	equiv_site = tuple(map(lambda s: s.strip().lower(), equiv_site))
	for j in range(3):
		xyz = equiv_site[j].replace('+',' +').replace('-',' -').split()
		for i in xyz:
			if i.strip("+-") == 'x':
				rotation[0,j] = float(i.strip('x')+"1")
			elif i.strip("+-") == 'y':
				rotation[1, j] = float(i.strip('y')+"1")
			elif i.strip("+-") == 'z':
				rotation[2, j] = float(i.strip('z')+"1")
			
			if i.strip("+-xyz") != "":
				trans_vec[j] = eval(i)
	
	# Combination of rotation and translation into a single transformation
	# is done in a 4x4 affine transformation matrix
	symmetry_operation = affine_map(rotation)
	symmetry_operation[:3,3] = trans_vec
	return symmetry_operation

class CIFParser(object):
	"""
	Collection of methods that parses CIF files based on PyCifRW.
    
	Attributes
	----------
	file : str
		Absolute path to the CIF file associated with the parser's target structure.

	Methods
	-------
	hall_symbol
		Hall symbol taken directly from the file or inferred from other information.

	lattice_vectors
		Lattice vectors of the crystal structure.

	atoms
		Generator of atoms making up the asymmetric unit cell

	symmetry_operators
		Generator of symmetry operators that relate the asymmetric unit cell and the unit cell.
	"""
	def __init__(self, filename):
		self.file = CifFile(datasource = filename)

	def __enter__(self):
		return self
	
	def __exit__(self, type, value, traceback):
		del self.file
	
	@lru_cache(maxsize = 1)
	def hall_symbol(self):
		""" Returns the Hall symbol """
		block = self.file[self.file.keys()[0]]	# TODO: find right block
		hall_number = None
		for tag in ['_symmetry_space_group_name_Hall','_space_group_name_Hall']:
			try:
				hall_symbol = block[tag]
			except:
				for tag in ['_symmetry_Int_Tables_number','_space_group_IT_number']:
					try:
						hall_symbol = Number2Hall[block[tag]]
					except:
						for tag in ['_symmetry_space_group_name_H-M','_space_group_name_H-M_alt']:
							try:
								h_m_symbol = block[tag].translate(string.maketrans("", ""),string.whitespace)
								hall_number = HM2Hall[h_m_symbol]
							except: 
								pass
		
		if hall_symbol is None:
			raise ParseError('Hall number could not be inferred')
	
		if hall_symbol[0] == "-":
			hall_symbol = "-" + hall_symbol[1].upper() + hall_symbol[2:].lower()
		else:
			hall_symbol = hall_symbol[0].upper() + hall_symbol[1:].lower()
		
		return hall_symbol

	@lru_cache(maxsize = 1)
	def lattice_vectors(self):
		""" 
		Returns the lattice vectors associated to a CIF structure.
        
		Returns
		-------
		lv : list of ndarrays, shape (3,)
		"""
		block = self.file[self.file.keys()[0]]
		
		try:
			a, _ = get_number_with_esd(block["_cell_length_a"])
			b, _ = get_number_with_esd(block["_cell_length_b"])
			c, _ = get_number_with_esd(block["_cell_length_c"])
			alpha, _ = get_number_with_esd(block["_cell_angle_alpha"])
			beta, _ = get_number_with_esd(block["_cell_angle_beta"])
			gamma, _ = get_number_with_esd(block["_cell_angle_gamma"])
		except:
			raise ParseError('Lattice vectors could not be determined.')

		return lattice_vectors_from_parameters(a, b, c, alpha, beta, gamma)
	
	def symmetry_operators(self):
		"""
		Returns the symmetry operators that map the atomic positions in a
		CIF file to the crystal unit cell.

		Yields
		------
		sym_ops : ndarray, shape (4,4)
			Transformation matrices. Since translations and rotation are combined,
			the transformation matrices are 4x4.
		"""
		block = self.file[self.file.keys()[0]]

		equivalent_sites_str = None
		for tag in ['_symmetry_equiv_pos_as_xyz','_space_group_symop_operation_xyz']:
			with suppress(KeyError):
				equivalent_sites_str = block.GetLoop(tag).get(tag) 

		if not equivalent_sites_str:
			warnings.warn('Equivalent sites note stored in the file. \
						   Pulling from database', UserWarning)
			equivalent_sites_str = SymOpsHall[self.hall_symbol()]
		
		# P1 space group only has a single equivalent site
		if isinstance(equivalent_sites_str, str):
			equivalent_sites_str = [equivalent_sites_str]

		yield from map(sym_ops, equivalent_sites_str)
	
	def atoms(self):
		"""
		Asymmetric unit cell. Combine with CIFParser.symmetry_operators() for a full unit cell.

		Yields
		------
		atoms : skued.structure.Atom instance
		"""
		block = self.file[self.file.keys()[0]]
		try:
			tmpdata = block.GetLoop('_atom_site_fract_x')
			cartesian = False
		except:
			try:
				tmpdata = block.GetLoop('_atom_site_Cartn_x')
				cartesian = True
			except:
				raise ParseError('Atomic positions could not be found or inferred.')
			
			t11 = block.get('_atom_sites_Cartn_tran_matrix_11')
			t12 = block.get('_atom_sites_Cartn_tran_matrix_12')
			t13 = block.get('_atom_sites_Cartn_tran_matrix_13')
			t21 = block.get('_atom_sites_Cartn_tran_matrix_21')
			t22 = block.get('_atom_sites_Cartn_tran_matrix_22')
			t23 = block.get('_atom_sites_Cartn_tran_matrix_23')
			t31 = block.get('_atom_sites_Cartn_tran_matrix_13')
			t32 = block.get('_atom_sites_Cartn_tran_matrix_23')
			t33 = block.get('_atom_sites_Cartn_tran_matrix_33')
			cart_trans_matrix_inv = np.array([[float(t11),float(t12),float(t13)],
										      [float(t21),float(t22),float(t23)],
										      [float(t31),float(t32),float(t33)]])
			cart_trans_matrix = np.linalg.inv(cart_trans_matrix_inv)

			if not all([t11, t12, t13, t21, t22, t23, t31, t32, t33]):
				raise ParseError('Cartesian coordinates in CIF but no transformation matrix given')
		
		if cartesian:
			xs = tmpdata.get('_atom_site_Cartn_x')
			ys = tmpdata.get('_atom_site_Cartn_y')
			zs = tmpdata.get('_atom_site_Cartn_z')
		else:
			xs = tmpdata.get('_atom_site_fract_x')
			ys = tmpdata.get('_atom_site_fract_y')
			zs = tmpdata.get('_atom_site_fract_z')
		# TODO: handle wildcards like '?', '.' in xs, ys, zs

		elements = tmpdata.get('_atom_site_type_symbol')
		if not elements:
			elements = tmpdata.get('_atom_site_label')
			if not elements:
				raise ParseError('Atom symbols could not be found or inferred.')
		elements = map(lambda s: s.strip(string.punctuation + string.digits).title(), elements)
		
		lv = self.lattice_vectors()
		for e, x, y, z in zip(elements, xs, ys, zs):
			coords = np.array([get_number_with_esd(x)[0], 
							   get_number_with_esd(y)[0], 
							   get_number_with_esd(z)[0]])

			if cartesian:
				coords = transform(cart_trans_matrix, coords)
			else:
				coords = real_coords(coords, lv)
			
			yield Atom(element = e, coords = coords)
