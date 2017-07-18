# -*- coding: utf-8 -*-

import numpy as np

from .lattice import Lattice
from .. import (change_of_basis, is_rotation_matrix, transform,
				translation_matrix)
from .scattering_params import scattering_params

# Constants
m = 9.109*10**(-31)     #in kg
a0 = 0.5291             #in Angs

def real_coords(frac_coords, lattice_vectors):
	"""
	Calculates the real-space coordinates of the atom from fractional coordinates and lattice vectors.
	
	Parameters
	----------
	frac_coords : array-like, shape (3,)
		Fractional coordinates
	lattice_vectors : list of ndarrays, shape (3,)
		Lattice vectors of the crystal.
		
	Returns
	-------
	coords : ndarray, shape (3,)
	"""
	COB = change_of_basis(lattice_vectors, np.eye(3))
	return transform(COB, frac_coords)

def frac_coords(real_coords, lattice_vectors):
	"""
	Calculates and sets the real-space coordinates of the atom from fractional coordinates and lattice vectors.
	Only valid for inorganic compounds.
	
	Parameters
	----------
	real_coords : array-like, shape (3,)
		Real-space coordinates
	lattice_vectors : list of ndarrays, shape (3,)
		Lattice vectors of the crystal.
		
	Returns
	-------
	coords : ndarray, shape (3,)
	"""
	COB = change_of_basis(np.eye(3), lattice_vectors)
	return transform(COB, real_coords)

class Atom(object):
	"""
	Container object for atomic data. 
	"""

	__slots__ = ('element', 'coords', 'displacement', 'lattice', '_a', '_b', '_c', '_d')

	# TODO: PDB identification?
	def __init__(self, element, coords, lattice = Lattice(), displacement = None, **kwargs): 
		"""
		Parameters
		----------
		element : str
			Chemical element
		coords : array-like, shape (3,)
			Coordinates of the atom in Euclidiant basis. See real_coords.
		lattice : Lattice instance, optional
			Lattice in which the atom is located. Default is the trivial lattice (i.e. Euclidian basis)
		displacement : array-like or None, optional
			Atomic maximum displacement [Angs]. If None (default), set to (0,0,0).
		"""

		if element not in atomic_number:
			raise ValueError('Invalid chemical element {}'.format(element))

		displacement = (0,0,0) or displacement  # If None, -> (0,0,0)
		
		self.element = element
		self.coords = np.array(coords, dtype = np.float)
		self.displacement = np.array(displacement, dtype = np.float)
		self.lattice = lattice
		
		# Atomic potential parameters loaded on instantiation
		# These are used to compute atomic potential
		try:
			_, a1, b1, a2, b2, a3, b3, c1, d1, c2, d2, c3, d3 = scattering_params[self.atomic_number]
		except KeyError:
			raise ValueError('Scattering information for element {} is unavailable.'.format(self.element))

		self._a = np.array((a1, a2, a3)).reshape((1,3))
		self._b = np.array((b1, b2, b3)).reshape((1,3))
		self._c = np.array((c1, c2, c3)).reshape((1,3))
		self._d = np.array((d1, d2, d3)).reshape((1,3))
		
	def __repr__(self):
		x, y, z = tuple(self.coords)
		return "< Atom {} at coordinates ({:.2f}, {:.2f}, {:.2f}) >".format(self.element, x, y, z)
	
	def __sub__(self, other):
		""" Returns the distance between the two atoms. """
		return np.linalg.norm(self.coords - other.coords)
	
	def __eq__(self, other):
		return (isinstance(other, self.__class__)
				and (self.element == other.element) 
				and np.allclose(self.coords, other.coords, atol = 1e-3) 
                and np.allclose(self.real_coords, other.real_coords, atol = 1e-3)
				and np.allclose(self.displacement, other.displacement, atol = 1e-3))
	
	def __hash__(self):
		return hash( (self.element, 
					  tuple(np.round(self.coords, 3)), 
                      tuple(np.round(self.real_coords, 3)),
					  tuple(np.round(self.displacement, 3))) )

	@property
	def atomic_number(self):
		return atomic_number[self.element]
	
	@property
	def weight(self):
		return atomic_weights[self.atomic_number - 1]

	@property
	def real_coords(self):
		return real_coords(self.coords, self.lattice.lattice_vectors)
	
	def electron_form_factor(self, nG):
		"""
		Vectorized electron form factor calculation.

		Parameters
		----------
		nG : array_like
			Scattering vector norm (G = 4 pi s)
		
		Returns
		-------
		atomff : `~numpy.ndarray`
		"""
		scatt_vector_norm = nG / (2*np.pi)	# In the units of Kirkland 2010
		
		s = scatt_vector_norm.shape
		scatt_vector_norm = scatt_vector_norm.reshape((-1,1))
		scatt_vector_norm2 = np.square(scatt_vector_norm)

		sum1 = np.sum(self._a/(scatt_vector_norm2 + self._b), axis = 1)
		sum2 = np.sum(self._c * np.exp(-self._d * scatt_vector_norm2), axis = 1)
		
		return (sum1 + sum2).reshape(s)
	
	def debye_waller_factor(self, G, out = None):
		"""
		Debye-Waller factor, calculated with the average of the 
		sinusoid displacement over a full cycle.

		Parameters
		----------
		G : ndarrays, ndim 3
		
		out : ndarray or None, optional
			NumPy ufunc parameter for avoiding unnecessary copies.

		Returns
		-------
		out : ndarray
		"""
		Gx, Gy, Gz = G
		dot = self.displacement[0]*Gx + self.displacement[1]*Gy + self.displacement[2]*Gz
		return np.exp(-0.5*dot**2, out = out)   # Factor of 1/2 from average of u = sin(wt)
	
	def potential(self,r):
		"""
		Electrostatic atomic potential.

		Parameters
		----------
		r : array_like
			Radial distance from the atom [Angs].

		Returns
		-------
		out : ndarray
			Potential [V * Angs]

		References
		----------
		Kirkland 2010 Eq. C.19
		"""
		r = np.array(r, copy = False)
		
		s = r.shape
		r = r.reshape((-1,1))
		sum1 = np.sum((self._a/r) * np.exp(-2*np.pi*r*np.sqrt(self._b)), axis = 1)
		sum2 = np.sum(self._c*self._d**(-1.5) * np.exp( -(r*np.pi)**2 / self._d), axis = 1)

		e = 14.4 #[Volt-Angstrom]
		res = 2*a0*e*(np.pi**2 * sum1 + np.pi**(2.5) * sum2)
		return res.reshape(s)
	
	def transform(self, *matrices):
		"""
		Transforms the real space coordinates according to a matrix.
		
		Parameters
		----------
		matrices : ndarrays, shape {(3,3), (4,4)}
			Transformation matrices.
		"""
		for matrix in matrices:
			self.coords = transform(matrix, self.coords)

atomic_number = {"H": 1, "He":2,  "Li":3,  "Be":4,  "B": 5,  "C": 6, "N": 7, "O":8,
				 "F": 9, "Ne":10, "Na":11, "Mg":12, "Al":13, "Si":14,"P": 15,
				 "S": 16,"Cl":17, "Ar":18, "K": 19, "Ca":20, "Sc":21,"Ti":22,
				 "V": 23,"Cr":24, "Mn":25, "Fe":26, "Co":27, "Ni":28,"Cu":29,
				 "Zn":30,"Ga":31, "Ge":32, "As":33, "Se":34, "Br":35,"Kr":36,
				 "Rb":37,"Sr":38, "Y": 39, "Zr":40, "Nb":41, "Mo":42,"Tc":43,
				 "Ru":44,"Rh":45, "Pd":46, "Ag":47, "Cd":48, "In":49,"Sn":50,
				 "Sb":51,"Te":52, "I": 53, "Xe":54, "Cs":55, "Ba":56,"La":57,
				 "Hf":72,"Ta":73, "W": 74, "Re":75, "Os":76, "Ir":77,"Pt":78,
				 "Au":79,"Hg":80, "Tl":81, "Pb":82, "Bi":83, "Po":84,"At":85,
				 "Rn":86,"Fr":87, "Ra":88, "Ac":89, "Rf":104,"Ce":58,"Pr":59,
				 "Nd":60,"Pm":61, "Sm":62, "Eu":63, "Gd":64, "Tb":65,"Dy":66,
				 "Ho":67,"Er":68, "Tm":69, "Yb":70, "Lu":71, "Th":90,"Pa":91,
				 "U": 92,"Np":93, "Pu":94, "Am":95, "Cm":96, "Bk":97,"Cf":98,
				 "Es":99,"Fm":100,"Md":101,"No":102,"Lr":103,"Rf":104,"Db":105,
				 "Sg":106,"Bh":107,"Hs":108,"Mt":109,"D":110,"Rg":111,"Cn":112,
				 "Uut":113,"Uuq":114,"Uup":115,"Uuh":116,"Uus":117,"Uuo":118}

#From atomdb.atomic
atomic_weights = ( 1.00794,   4.002602,   6.941    ,   9.012182 ,  10.811    ,
			12.0107, 14.0067, 15.9994, 18.9984032,  20.1797   ,
			22.98976928,  24.3050  ,  26.9815386,  28.0855   ,  30.973762 ,
			32.065     ,  35.453   ,  39.948    ,  39.0983   ,  40.078    ,
			44.955912  ,  47.867   ,  50.9415   ,  51.9961   ,  54.938045 ,
			55.845     ,  58.933195,  58.6934   ,  63.546    ,  65.38     ,
			69.723     ,  72.64    ,  74.92160  ,  78.96     ,  79.904    ,
			83.798     ,  85.4678  ,  87.62     ,  88.90585  ,  91.224    ,
			92.90638   ,  95.96    ,  98.000    , 101.07     , 102.90550  ,
		   106.42      , 107.8682  , 112.411    , 114.818    , 118.710    ,
		   121.760     , 127.60    , 126.90447  , 131.293    , 132.9054519,
		   137.327     , 138.90547 , 140.116    , 140.90765  , 144.242    ,
		   145.000     , 150.36    , 151.964    , 157.25     , 158.92535  ,
		   162.500     , 164.93032 , 167.259    , 168.93421  , 173.054    ,
		   174.9668    , 178.49    , 180.94788  , 183.84     , 186.207    ,
		   190.23      , 192.217   , 195.084    , 196.966569 , 200.59     ,
		   204.3833    , 207.2     , 208.98040  , 209.000    , 210.000    ,
		   222.000     , 223.000   , 226.000    , 227.00     , 232.03806  ,
		   231.03588   , 238.02891)
