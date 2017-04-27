# -*- coding: utf-8 -*-

import numpy as np
from .form_factors import params
from . import Transformable
from .. import  change_of_basis, transform, translation_matrix, is_rotation_matrix
from scipy.special import k0 as bessel

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

class Atom(Transformable):
	"""
	Container object for atomic data. 
           
	Attributes
	----------
	element : str
		Chemical symbol of the atom element
	atomic_number : int
		Atomic number of the object.
	weight : float
		Standard atomic weight (in international units 'u').
	coords : ndarray, shape (3,)
		Real-space coordinates of the atom.
	identification : str
		Protein DataBank identification.

	Methods
	-------
	rotate
		Rotate coordinates.
    
	translate
		Translate coordinates.
    
	transform
		Apply 3x3 or 4x4 transformation matrices for reflections, shearing, 
		translations, rotations, projections, etc.

	debye_waller_factor
		Diffracted intensity suppression due to mean atomic displacement from
		perfect crystalline positions.
    
	potential
		Electrostatic potential computed on a three-dimensional meshgrid

	projected_potential
		Projection of the electrostatic potential onto the x-y plane.

	Special methods
	---------------
	Atom[index] is the same as Atom.coords[index]

	Atom1 - Atom2 gives the distance between atoms in Angs
	"""

	__slots__ = ('element', 'coords', 'identification', 'displacement', 
					'_a', '_b', '_c', '_d')

	# TODO: PDB identification?
	def __init__(self, element, coords, displacement = None, **kwargs): 
		"""
		Parameters
		----------
		element : str
			Chemical element
		coords : array-like, shape (3,)
			Coordinates of the atom in Euclidiant basis. See real_coords.
		displacement : array-like or None, optional
			Atomic maximum displacement [Angs]. If None (default), set to (0,0,0).
		"""
		super().__init__(**kwargs)

		if element not in atomic_number:
			raise ValueError('Invalid chemical element {}'.format(element))

		displacement = (0,0,0) or displacement  # If None, -> (0,0,0)
		
		self.element = element
		self.coords = np.array(coords, dtype = np.float)
		self.displacement = np.array(displacement, dtype = np.float)
		
		# Atomic potential parameters loaded on instantiation
		# These are used to compute atomic potential
		# TODO: add ref Kirkland 2008
		_, a1, b1, a2, b2, a3, b3, c1, d1, c2, d2, c3, d3 = params[self.atomic_number]
		self._a = np.array((a1, a2, a3)).reshape((1,3))
		self._b = np.array((b1, b2, b3)).reshape((1,3))
		self._c = np.array((c1, c2, c3)).reshape((1,3))
		self._d = np.array((d1, d2, d3)).reshape((1,3))
		
	def __repr__(self):
		return "< Atom {0} at coordinates {1} >".format(self.element, repr(self.coords))
	
	def __sub__(self, atom):
		""" Returns the distance between the two atoms. """
		return np.linalg.norm(self.coords - atom.coords)
	
	def __getitem__(self, index):
		return self.coords[index]
	
	def __setitem__(self, index, value):
		self.coords[index] = value

	@property
	def atomic_number(self):
		return atomic_number[self.element]
	
	@property
	def weight(self):
		return atomic_weights[self.atomic_number - 1]
	
	def frac_coords(self, lattice_vectors):
		"""
		Returns the fractional coordinates within a lattice.

		Parameters
		----------
		lattice_vectors : iter of 3 ndarrays, shape (3,)
			Lattice basis vectors

		Returns
		-------
		out : ndarray, shape (3,)
			Fractional coordinates
		"""
		return frac_coords(self.coords, lattice_vectors)
	
	def form_factor(self, scatt_vector_norm):
		""" 
		Vectorized atomic form factor.
        
		Parameters
		----------           
		scatt_vector_norm : array-like of numericals
			Scattering vector length |G|.
        
		Returns
		-------
		atomff : array-like of numerical
			array of the same shape as input
        
		Notes
		-----
		By convention, scattering vectors G are defined such that |G| = 4 pi s
		"""
		# TODO: vectorize better
		# TODO: compute as fourier transform of potential?
		s = np.reshape(scatt_vector_norm/(4*np.pi), newshape = (1,scatt_vector_norm.size))
		a, b = (np.expand_dims(item, axis = 1) for item in atomic_ff_dict[self.element])
		ff = np.sum(a*np.exp(-b*(s**2)), axis = 0)
		return np.reshape(ff, scatt_vector_norm.shape)
	
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
		r : ndarray
			Radial distance from the atom [Angs].

		Returns
		-------
		out : ndarray
			Potential [V * Angs]

		References
		----------
		Kirkland 2010 Eq. C.19
		"""
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
		super().transform(*matrices)

	def rotate(self, *matrices):
		"""
		Rotates the real-space coordinates according to a rotation matrix. This
		method is a shortcut of transform()
        
		Parameters
		----------
		matrices : ndarrays, shape (3,3)
			Rotation matrices
        
		See also
		--------
		transform
			More general transformations, including rotations, translations, 
			reflections, shearing.
		"""
		for matrix in matrices:
			try:
				assert(is_rotation_matrix(matrix))
			except AssertionError:
				raise ValueError('One of the matrices provided is not a rotation matrix.')
		
		self.transform(*matrices)
	
	def translate(self, *vectors):
		""" 
		Translates coordinates according to a vector.
        
		Parameters
		----------
		vectors : arrays-like, shape (3,)
		"""
		for vector in vectors:
			self.transform(translation_matrix(vector))


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

# Atomic form factors coefficients
# TODO: use Zheng et al. 2009's data that uses 8-term series of gaussian
# TODO: all atoms
a_H = np.array([0.0349,0.1201,0.1970,0.0573,0.1195]);b_H = np.array([0.5347,3.5867,12.3471,18.9525,38.6269])
a_V = np.array([0.2969,1.0774,2.1894,3.0825,1.7190]);b_V = np.array([0.1505,1.6392,7.5691,36.8741,107.8517])
a_O = np.array([0.0365,0.1729,0.5805,0.8814,0.3121]);b_O = np.array([0.0652,0.6184,2.9449,9.6298,28.2194])
a_C = np.array([0.0489,0.2091,0.7537,1.1420,0.3555]);b_C = np.array([0.1140,1.0825,5.4281,17.8811,51.1341])
a_N = np.array([0.0267,0.1328,0.5301,1.1020,0.4215]);b_N = np.array([0.0541,0.5165,2.8207,10.6297,34.3764])
a_O = np.array([0.0365,0.1729,0.5805,0.8814,0.3121]);b_O = np.array([0.0652,0.6184,2.9449,9.6298,28.2194])
a_S = np.array([0.0915,0.4312,1.0847,2.4671,1.0852]);b_S = np.array([0.0838,0.7788,4.3462,15.5846,44.6365])
a_P = np.array([0.1005,0.4615,1.0663,2.5854,1.2725]);b_P = np.array([0.0977,0.9084,4.9654,18.5471,54.3648])
a_Cu= np.array([0.4314,1.3208,1.5236,1.4671,0.8562]);b_Cu= np.array([0.2694,1.9223,7.3474,28.9892,90.6246])
a_Br= np.array([0.4798,1.1948,1.8695,2.6953,0.8203]);b_Br= np.array([0.2504, 1.5963,6.9653,19.8492,50.3233])
a_Au= np.array([0.3055,1.3956,2.9617,3.8990,2.0026]);b_Au= np.array([0.0596,0.5827,3.1035,11.9693,47.9106])
a_Pb= np.array([0.3540,1.5453,3.5975,4.3152,2.7743]);b_Pb= np.array([0.0668,0.6465,3.6968,16.2056,61.4909])
a_U = np.array([0.6410,2.2643,4.8713,5.9287,5.3935]);b_U = np.array([0.1097,1.0644,5.7907,25.0261,101.3899])

#Atomic Form Factors switch-case implementation
atomic_ff_dict = {'H':(a_H, b_H), 
                  'C': (a_C, b_C), 
                  'N': (a_N, b_N), 
                  'O': (a_O, b_O),
                  'P': (a_P, b_P), 
                  'S': (a_S, b_S), 
                  'Au':(a_Au, b_Au),
                  'Pb':(a_Pb, b_Pb),
                  'U': (a_U, b_U), 
                  'V': (a_V, b_V),
                  'Cu':(a_Cu, b_Cu),
                  'Br':(a_Br,b_Br)}