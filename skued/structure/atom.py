# -*- coding: utf-8 -*-

from functools import lru_cache

import numpy as np

from .. import change_of_basis, transform
from .atom_data import (ELEM_TO_MAGMOM, ELEM_TO_MASS, ELEM_TO_NAME,
                        ELEM_TO_NUM, NUM_TO_ELEM)
from .lattice import Lattice


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
    COB = change_of_basis(np.array(lattice_vectors), np.eye(3))
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
    COB = change_of_basis(np.eye(3), np.array(lattice_vectors))
    return np.mod(transform(COB, real_coords), 1)

# TODO: store atomic data as class attributes?
class Atom(object):
    """
    Container object for atomic data. 

    Parameters
    ----------
    element : str or int
        Chemical element symbol or atomic number.
    coords : array-like, shape (3,)
        Coordinates of the atom in fractional form.
    displacement : array-like or None, optional
        Atomic maximum displacement [Angs].
    magmom : float, optional
        Magnetic moment. If None (default), the ground-state magnetic moment is used.
    """
    __slots__ = ('element', 'coords', 'displacement', 'magmom', 
                 '_a', '_b', '_c', '_d')

    def __init__(self, element, coords, displacement = (0,0,0), magmom = None, **kwargs): 
        if isinstance(element, int):
            element = NUM_TO_ELEM[element]
        elif element not in ELEM_TO_NUM:
            raise ValueError('Invalid chemical element {}'.format(element))
        
        if magmom is None:
            magmom = ELEM_TO_MAGMOM[element]
        
        self.element = element
        self.coords = np.array(coords, dtype = np.float)
        self.displacement = np.array(displacement, dtype = np.float)
        self.magmom = magmom
        
    def __repr__(self):
        return "< Atom {:<2} @ ({:.2f}, {:.2f}, {:.2f}) >".format(self.element, *tuple(self.coords))
    
    # TODO: add `distance_from` function for atoms on a lattice
    def __sub__(self, other):
        return np.linalg.norm(self.coords - other.coords)
    
    def __eq__(self, other):
        return (isinstance(other, self.__class__)
                and (self.element == other.element) 
                and (self.magmom == other.magmom)
                and np.allclose(self.coords, other.coords, atol = 1e-3) 
                and np.allclose(self.displacement, other.displacement, atol = 1e-3))
    
    def __hash__(self):
        return hash( (self.element, 
                      self.magmom,
                      tuple(np.round(self.coords, 3)), 
                      tuple(np.round(self.displacement, 3))) )
    
    @classmethod
    def from_ase(cls, atom):
        """ 
        Returns an Atom instance from an ASE atom 
        
        Parameters
        ----------
        atom : ase.Atom
        """
        lattice = np.eye(3)
        if atom.atoms is not None:
            lattice = np.array(atom.atoms.cell)

        return cls(element = atom.symbol, 
                   coords = frac_coords(atom.position, lattice), 
                   magmom = atom.magmom)

    @property
    def atomic_number(self):
        return ELEM_TO_NUM[self.element]
    
    @property
    def mass(self):
        return ELEM_TO_MASS[self.element]

    def ase_atom(self, lattice = None, **kwargs):
        """
        Returns an ``ase.Atom`` object. 
        
        Parameters
        ----------
        lattice : skued.Lattice, optional
            Lattice on which the atoms sit. Default is free space.
        kwargs
            Keyword arguments are passed to the ``ase.Atom`` constructor.
        
        Returns
        -------
        atom: ase.Atom

        Raises
        ------
        Import : If ASE is not installed
        """
        from ase import Atom

        if lattice is None:
            lattice = Lattice(np.eye(3))

        return Atom(symbol = self.element, 
                    position = self.xyz(lattice), 
                    magmom = self.magmom,
                    mass = self.mass, **kwargs)

    @lru_cache()
    def xyz(self, lattice):
        """ 
        Real-space position of the atom

        Parameters
        ----------
        lattice : Lattice or iterable
            Lattice or Crystal instance in which the atom is located.
                    
        Returns
        -------
        pos : `~numpy.ndarray`, shape (3,)
            Atomic position
        """
        return real_coords(self.coords, lattice.lattice_vectors)
    
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
