"""
CIF 1.0, 1.1 and 2.0 parser based on PyCifRW.

@author : Laurent P. Rene de Cotret 
"""

from . import lattice_vectors_from_parameters, Atom, real_coords 
from .. import translation_matrix, affine_map, change_of_basis
import numpy as n
import spglib

class CIFParser(object):
	"""
	Collection of methods that parses CIF 1.1 files. Modified from mmLib.CIF module.
    
	Attributes
	----------
	file : str
		Absolute path to the CIF file associated with the parser's target structure.

	Methods
	-------
	space_group
		Space-group of the structure.   
    
	lattice_vectors
		Lattice vectors of the crystal structure.

	atoms
		Returns a list of atoms making up the asymmetric unit cell

	symmetry_operators
		Symmetry operators that relate the asymmetric unit cell and the unit cell.
	"""
	def __init__(self, filename):
		if not WITH_PYCIFRW:
			raise ImportError('PyCifRW must be installed to take advantage of the CIF Parser')
		self.file = CifFile(datasource = filename)

	def __enter__(self):
		return self
	
	def __exit__(self):
		try:
			self.file.close()
		except AttributeError:
			pass
	
	def space_group(self):
		"""
		Returns the Hall number for the structure space group.
		Not to be confused with Hall Symbol.

		Returns
		-------
		spg : int   
			Hall number

		Raises
		------
		IOError
			If space group could not be parsed.
		"""
		block = self.file.first_block()
		try:
			hall_symbol = block['_space_group_name_Hall'].replace(' ','')
		except KeyError:
			hall_symbol = block['_symmetry_space_group_name_Hall'].replace(' ','')      # Legacy
		except KeyError:
			raise IOError('Space group Hall symbol not contained in file.')
		
		# Iterate over all hall numbers (1 to 530) to get the corresponding space group
		for hall_number in range(1, 531):
			if spglib.get_spacegroup_type(hall_number)['hall_symbol'].replace(' ','') == hall_symbol:
				return hall_number
		
		raise IOError('Space group {} could not be found within the spglib databse.'.format(hm_symbol))

	def lattice_vectors(self):
		""" 
		Returns the lattice vectors associated to a CIF structure.
        
		Returns
		-------
		lv : list of ndarrays, shape (3,)
		"""
		block = self.file.first_block()
		
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
        
		Returns
		-------
		sym_ops : iterable
		"""
		dataset = spglib.get_symmetry_from_database(self.space_group())

		# As returned by spglib, rotations and translations act on the fractional coordinates.
		# A change of basis is required.
		COB = change_of_basis(basis1 = self.lattice_vectors(), basis2 = standard_basis)
		sym_ops = list()
		for r, t in zip(dataset['rotations'], dataset['translations']):
			# rotations and translations combined as in http://www.euclideanspace.com/maths/geometry/affine/matrix4x4/
			operator = affine_map(COB @ r)
			operator[:3,3] = COB @ t
			sym_ops.append(operator)

		return sym_ops
	
	def atoms(self):
		"""
		Asymmetric unit cell.

		Returns
		-------
		atoms : iterable of uediff.structure.Atom
		"""
		block = self.file.first_block()
		elements = block['_atom_site_type_symbol']
		xs, ys, zs = block['_atom_site_fract_x'], block['_atom_site_fract_y'], block['_atom_site_fract_y']
		coords = [n.array([float(x), float(y), float(z)]) for x,y,z in zip(xs, ys, zs)]

		return [Atom(element = element, coords = real_coords(frac_coords = coord, lattice_vectors = self.lattice_vectors()))
				for element, coord in zip(elements, coords)]