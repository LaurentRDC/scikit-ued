# -*- coding: utf-8 -*-
"""
@author: Laurent P. René de Cotret
"""
import numpy as np
from . import Transformable
from .. import is_rotation_matrix, rotation_matrix, translation_matrix

class Structure(Transformable):
	"""
	Attributes
	----------
	items

	Methods
	-------
	transform
		Transform the coordinates of children atoms or structures using a 3x3 or 4x4 transformation matrix.
    
	rotate
		Convenience method to the transform method. Rotates children atoms and structures according to a rotation
		matrix.

	translate
		Convenience method to the transform method. Translates children atoms and structures according
		to a translation vector.
	"""

	def __init__(self, items, **kwargs):
		"""
		Parameters
		----------
		items : iterable
			List or tuple of items to be contained in this object. Items must be subclasses of Transformable.
		"""
		super().__init__(**kwargs)
		self.items = list(items)

	def find_item(self, value, key):
		"""
		Search for an item by a key attrbute.

		Parameters
		----------
		value : 
			Value to find
		key : str
			Attribute name.
        
		Returns
		-------
		item or None
		"""
		return next( (item for item in self.items if getattr(item, key) == value), None)
	
	def sort(self, attribute):
		"""
		Sort items by attributes

		Parameters
		----------
		attribute : str
		"""
		self.items.sort(key = lambda x: getattr(x, attribute))
	
	def transform(self, *matrices):
		"""
		Transforms the real-space coordinates of children atoms according to a transformation matrix.
        
		Parameters
		----------
		matrix : ndarray, shape {(3,3), (4,4)}
			Transformation matrix.
		"""
		for item in self.items:
			item.transform(*matrices)
		super().transform(*matrices)
	
	def rotate(self, angle, axis):
		""" 
		Rotate all the atoms around the origin [0,0,0] according to the input 
		matrices. Shortcut to using the method 'transform'.
        
		Parameters
		----------
		angle : float
			angle in deg
		axis : array-like of length 3
			Axis about which to rotate
		"""
		self.transform(rotation_matrix(angle = angle * (np.pi/180), axis = axis))
	
	def translate(self, *vectors):
		""" 
		Translates every atom in the residue. Shortcut to using the method 'transform'.
        
		Parameters
		----------
		vectors : array-like, shape (3,)
			Coordinates added to all atoms coordinates in the residue.
		"""
		for vector in vectors:
			self.transform(translation_matrix(vector))

class AtomicStructure(Structure):
	"""
	Container object for all structures made of atoms.
    
	Attributes
	----------
	items : list
		List of Transformable objects: either AtomicStructures or Atom objects.
    
	atoms : list of uediff.structure.Atom objects
		List of Atom present in this AtomicStructure and all children AtomicStructures.
	"""

	# alias for when items are atoms
	@property
	def atoms(self):
		return self.items

	def find_atom(self, value, attribute):
		"""
		Search for an atom by a key attrbute.

		Parameters
		----------
		key : str, optional
			structure.Atom attribute name. Default is 'identification' (PDB atom identification)
        
		Returns
		-------
		item or None
		"""
		return self.find_item(value = value, key = attribute)

class CompositeStructure(AtomicStructure):
    
	@property
	def atoms(self):
		return [atom for structure in self.items for atom in structure.items]