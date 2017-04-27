from copy import copy
from collections.abc import Iterable
from functools import lru_cache
from itertools import product, takewhile, count
import numpy as np
from numpy import pi
from numpy.linalg import norm
from operator import mul
from scipy.special import k0 as bessel
import spglib
from warnings import warn

from . import AtomicStructure, Lattice, CIFParser, PDBParser, real_coords
from .. import change_of_basis, transform, affine_map, change_basis_mesh, is_rotation_matrix

# Constants
m = 9.109*10**(-31)     #in kg
a0 = 0.5291             #in Angs
e = 14.4                #Volt*Angstrom

def symmetry_expansion(items, *symmetry_operators):
    """
    Expands a a list of Transformable objects according to symmetry operators.
    
    Parameters
    ----------
    items : list of uediff.structure.Transformable objects 
        Typically AtomicStructures or Atoms
    symmetry_operators : tuple of ndarrays, shape (3,3)
        Symmetry operators
    
    Returns
    -------
    expanded : generator of uediff.structure.Transformable objects
    """
    for item in items:
        for sym_op in symmetry_operators:
            new_item = copy(item)
            new_item.transform(sym_op)
            yield new_item

class Crystal(AtomicStructure, Lattice):
	"""
	This object is the basis for inorganic crystals such as VO2, 
	and protein crystals such as bR. 

	Special methods
	---------------
	__iter__
		Generator of atoms in the unit cell.

	__len__
		Number of atoms in the unit cell.
    
	Attributes
	----------
	symmetry_operators : list of ndarrays
		Symmetry operators that links the underlying AtomicStructure to the unit cell construction.
    
	unitcell : list of Atom objects
		List of atoms in the crystal unitcell. iter(Crystal) is a generator that yields
		the same atoms; this approach is preferred.
    
	atoms : list
		List of atoms in the asymmetric unit.
    
	volume : float
		Volume of the unit cell, in angstrom**3
    
	spglib_cell : tuple
		Crystal structure in spglib's `cell` format.
    
	Constructors
	------------    
	from_cif
		Return a Crystal object from a CIF file. CIF versions 1.0, 1.1 and 2.0 are supported.

	from_pdb
		Returns a Crystal object from a Protein DataBank ID number. The .pdb file will
		be downloaded, cached, and parsed appropriately.
    
	Methods
	-------    
	periodicity
		Bounding cube of the unit cell in real-space. 
        
	potential
		Atomic potential computed over euclidian space.

	projected_potential
		Atomic potential computed as the projection onto the x-y plane.
    
	scattering_vector
		Scattering vector G from Miller indices (hkl)
    
	miller_indices
		Miller indices (hkl) from scattering vector G

	structure_factor
		General static structure factor calculation. Includes Debye-Waller suppression.
    
	structure_factor_miller
		Static structure factor calculation from Miller indices, taking into account
		the diffraction conditions of the crystal geometry.
    
	bounded_reflections
		Generate a set of (hkl) reflections below a bound.
    
	intensity_normalization
		Sum of form factor squared. Useful when normalizing diffraction intensities.
     
	transform
		Apply 3x3 or 4x4 transformation (reflection, shearing, translation, rotation)
		to atomic coordinates and lattice vectors.
	"""

	def __init__(self, atoms, symmetry_operators = [np.eye(3)], **kwargs):
		"""
		Parameters
		----------
		atoms : iterable of Atoms
			Atoms in the asymmetric cell.
        lattice_vectors : list of ndarrays, shape (3,), optional
            Lattice vectors. Default is a cartesian lattice.
		symmetry_operators : list of ndarrays, shape (3,3), optional
			Symmetry operators linking the underlying AtomicStructure to the unit cell construction. Default is
			the identity transformation.
		"""
		kwargs.update({'items': atoms}) # atoms argument is an alias for AtomicStructure.items
		self.symmetry_operators = tuple(map(affine_map, symmetry_operators))
		super().__init__(**kwargs)
	
	def __iter__(self):
		yield from symmetry_expansion(self.atoms, *self.symmetry_operators)
	
	def __len__(self):
		return len(tuple(iter(self)))
	
	def __repr__(self):
		return '< Crystal object with unit cell of {} atoms>'.format(len(self))

	@classmethod
	def from_cif(cls, path):
		"""
		Returns a Crystal object created from a CIF 1.0, 1.1 or 2.0 file.

		Parameters
		----------
		path : str
			File path
		"""
		with CIFParser(filename = path) as parser:
			return Crystal(items = parser.atoms(), 
							lattice_vectors = parser.lattice_vectors(), 
							symmetry_operators = parser.symmetry_operators())

	@classmethod
	def from_pdb(cls, ID):
		"""
		Returns a Crystal object created from a Protein DataBank entry.

		Parameters
		----------
		ID : str
			Protein DataBank identification. The correct *.pdb file will be downloaded,
			cached and parsed.
		"""
		parser = PDBParser(ID = ID)
		return Crystal(atoms = parser.atoms(), 
					   lattice_vectors = parser.lattice_vectors(),
					   symmetry_operators = parser.symmetry_operators())
	
	@property
	def unitcell(self):
		return list(iter(self))
	
	@property
	def spglib_cell(self):
		""" Returns the crystal structure in spglib's `cell` format."""
		lattice = np.array(self.lattice_vectors)
		positions = np.array([atom.frac_coords(self.lattice_vectors) for atom in iter(self)])
		numbers = np.array(tuple(atom.atomic_number for atom in iter(self)))
		return (lattice, positions, numbers)
	
	def periodicity(self):
		"""
		Crystal periodicity in x, y and z direction from the lattice constants.
		This is effectively a bounding cube for the unit cell, which is itself a unit cell.

		Parameters
		----------
		lattice : Lattice

		Returns
		-------
		out : tuple
			Periodicity in x, y and z directions [angstroms]
        
		Notes
		-----
		Warning: the periodicity of the lattice depends on its orientation in real-space.
		"""
		# By definition of a lattice, moving by the projection of all Lattice
		# vectors on an axis should return you to an equivalent lattice position
		e1, e2, e3 = np.eye(3)
		per_x = sum( (abs(np.vdot(e1,a)) for a in self.lattice_vectors) )
		per_y = sum( (abs(np.vdot(e2,a)) for a in self.lattice_vectors) )
		per_z = sum( (abs(np.vdot(e3,a)) for a in self.lattice_vectors) )
		return per_x, per_y, per_z
		
	def potential(self, x, y, z):
		"""
		Scattering potential calculated on a real-space mesh.

		Parameters
		----------
		x, y, z : ndarrays
			Real space coordinates mesh. 
        
		Returns
		-------
		potential : ndarray, dtype float
			Linear superposition of atomic potential [V*Angs]
		"""
		# TODO: multicore
		potential = np.zeros_like(x, dtype = np.float)
		r = np.zeros_like(x, dtype = np.float)
		for atom in self:
			ax, ay, az = atom.coords
			r[:] = np.sqrt( (x - ax)**2 + (y - ay)**2 + (z - az)**2 )
			potential += atom.potential(r)
		
		# Due to sampling, x,y, and z might pass through the center of atoms
		# Replace np.inf by the next largest value
		potential[np.isinf(potential)] = np.max(potential[np.isfinite])
		return potential
	
	def scattering_vector(self, h, k, l):
		"""
		Returns the scattering vector G from Miller indices.
        
		Parameters
		----------
		h, k, l : int or ndarrays
			Miller indices.
        
		Returns
		-------
		G : array-like
			If `h`, `k`, `l` are integers, returns a single array of shape (3,)
			If `h`, `k`, `l` are arrays, returns three arrays Gx, Gy, Gz
		"""
		if isinstance(h, Iterable):
			return change_basis_mesh(h, k, l, basis1 = self.reciprocal_vectors, basis2 = np.eye(3))

		b1,b2,b3 = self.reciprocal_vectors
		return int(h)*b1 + int(k)*b2 + int(l)*b3
	
	def miller_indices(self, G):
		"""
		Returns the miller indices associated with a scattering vector.
        
		Parameters
		----------
		G : array-like, shape (3,)
			Scattering vector.
        
		Returns
		-------
		hkl : ndarray, shape (3,), dtype int
			Miller indices [h, k, l].
		"""
		G = np.asarray(G, dtype = np.float)
	
		# Transformation matric between reciprocal space and miller indices
		matrix_trans = np.empty(shape = (3,3), dtype = np.float)
		for i in range(len(self.reciprocal_vectors)):
			matrix_trans[:,i] = self.reciprocal_vectors[i]

		matrix_trans = np.linalg.inv(matrix_trans)
		return transform(matrix_trans, G).astype(np.int)
	
	def structure_factor_miller(self, h, k, l, normalized = False):
		"""
		Computation of the static structure factor from Miller indices.
        
		Parameters
		----------
		h, k, l : array_likes
			Miller indices. Can be given in a few different formats:
            
			``3 floats``
				returns structure factor computed for a single scattering vector
                
			``list of 3 coordinate ndarrays, shapes (L,M,N)``
				returns structure factor computed over all coordinate space
        
		normalized : bool
			If True, returns the normalized structure factor E.
			See http://www.mx.iucr.org/iucr-top/comm/cteach/pamphlets/17/node4.html
        
		Returns
		-------
		sf : ndarray, dtype complex
			Output is the same shape as h, k, or l.
        
		See also
		--------
		structure_factor
			Vectorized structure factor calculation for general scattering vectors.	
		"""
		return self.structure_factor(G = self.scattering_vector(h, k, l), normalized = normalized)
		
	def structure_factor(self, G, normalized = False):
		"""
		Computation of the static structure factor. This function is meant for 
		general scattering vectors, not Miller indices. 
        
		Parameters
		----------
		G : array-like
			Scattering vector. Can be given in a few different formats:
            
			``array-like of numericals, shape (3,)``
				returns structure factor computed for a single scattering vector
                
			``list of 3 coordinate ndarrays, shapes (L,M,N)``
				returns structure factor computed over all coordinate space
            
			WARNING: Scattering vector is not equivalent to the Miller indices.
        
		normalized : bool
			If True, returns the normalized structure factor E.
			See http://www.mx.iucr.org/iucr-top/comm/cteach/pamphlets/17/node4.html
        
		Returns
		-------
		sf : ndarray, dtype complex
			Output is the same shape as input G[0]. Takes into account
			the Debye-Waller effect.
        
		See also
		--------
		structure_factor_miller 
			For structure factors calculated from Miller indices.
		        
		Notes
		-----
		By convention, scattering vectors G are defined such that |G| = 4 pi s
		"""
		# Distribute input
		# This works whether G is a list of 3 numbers, a ndarray shape(3,) or 
		# a list of meshgrid arrays.
		Gx, Gy, Gz = G
		nG = np.sqrt(Gx**2 + Gy**2 + Gz**2)
		
		# Separating the structure factor into sine and cosine parts avoids adding
		# complex arrays together. About 3x speedup vs. using complex exponentials
		SFsin, SFcos = np.zeros(shape = nG.shape, dtype = np.float), np.zeros(shape = nG.shape, dtype = np.float)

		# Pre-allocation
		normalization = 1.0
		atomff_dict = self._atomic_ff_dict(nG)  #Precalculation of form factors
		dwf = np.empty_like(SFsin)

		for atom in self: #TODO: implement in parallel_sum?
			x, y, z = atom.coords
			arg = x*Gx + y*Gy + z*Gz
			atom.debye_waller_factor((Gx, Gy, Gz), out = dwf)
			atomff = atomff_dict[atom.element]
			SFsin += atomff * dwf * np.sin(arg)
			SFcos += atomff * dwf * np.cos(arg)
		
		if normalized:
			normalization = 1/np.sqrt(self.intensity_normalization(nG))
		
		return normalization*(SFcos + 1j*SFsin)
	
	def bounded_reflections(self, nG):
		"""
		Returns iterable of reflections (hkl) with |G| < nG
        
		Parameters
		----------
		nG : float
			Maximal scattering vector norm. By our convention, |G| = 4 pi s.
        
		Returns
		-------
		h, k, l : ndarrays, shapes (N,), dtype int
		"""
		if nG < 0:
			raise ValueError('Bound {} is negative.'.format(nG))
		
		# Determine the maximum index such that (i00) family is still within data limits
		#TODO: cache results based on max_index?
		bounded = lambda i : any([norm(self.scattering_vector(i,0,0)) <= nG, 
									norm(self.scattering_vector(0,i,0)) <= nG, 
									norm(self.scattering_vector(0,0,i)) <= nG])
		max_index = max(takewhile(bounded, count(0)))
		extent = range(-max_index, max_index + 1)
		h, k, l = np.split(np.array(list(product(extent, extent, extent)), dtype = np.int), 3, axis = -1)
		h, k, l = h.ravel(), k.ravel(), l.ravel()

		# we only have an upper bound on possible reflections
		# Let's filter down
		Gx, Gy, Gz = self.scattering_vector(h, k, l)
		norm_G = np.sqrt(Gx**2 + Gy**2 + Gz**2)
		in_bound = norm_G <= nG
		return h.compress(in_bound), k.compress(in_bound), l.compress(in_bound)
	
	def intensity_normalization(self, nG):
		""" 
		Vectorized sum of form factor squared.
        
		Parameters
		----------           
		scatt_vector_norm : array-like of numericals
			Scattering vector length |G|.
        
		Returns
		-------
		total : array-like of numerical
			array of the same shape as input.
        
		Notes
		-----
		By convention, scattering vectors G are defined such that |G| = 4 pi s
		"""
		return sum(atom.form_factor(nG)**2 for atom in self)
	
	def transform(self, *matrices):
		"""
		Transforms the real space coordinates according to a matrix.
        
		Parameters
		----------
		matrices : ndarrays, shape {(3,3), (4,4)}
			Transformation matrices.
		"""
		# Only rotation matrices should affect the symmetry operations
		for matrix in matrices:
			matrix = affine_map(matrix)  # Make sure it is 4x4
			if is_rotation_matrix(matrix):
				matrix[:3, 3] = 0  # remove translations
			self.symmetry_operators = tuple(transform(matrix, sym_op) for sym_op in self.symmetry_operators)
		
		super().transform(*matrices)

	def _atomic_ff_dict(self, scatt_vector_norm):
		""" 
		Returns a dictionary containing the vectorized atomic form factor
		of each variety of atom in a crystal unit cell. Using this function 
		avoids recalculating atomic form factors over and over.
        
		Parameters
		----------            
		scatt_vector_norm : array-like of numericals
			Scattering vector length |G|
        
		Returns
		-------
		atomff_fict : dict
			Dictionnary with keys as atomic elements (e.g. 'V', 'H', ...) and
			values as the atomic form factor at a certain scattering vector
			norm.
		"""
		atomff_dict = dict()
		for atom in self.atoms: # No need to check beyond the irreducible group of atoms
			if atom.element not in atomff_dict:
				atomff_dict[atom.element] = atom.form_factor(scatt_vector_norm)
		return atomff_dict