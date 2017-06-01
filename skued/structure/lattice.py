# -*- coding: utf-8 -*-
from math import sin, cos, tan, sqrt, radians
import numpy as np
from numpy.linalg import norm
from . import Transformable
from .. import transform

e1, e2, e3 = np.eye(3) # Euclidian basis

# TODO: Introduce conventions on ordering a, b, c and angles
#       based on http://atztogo.github.io/spglib/definition.html
def lattice_vectors_from_parameters(a, b, c, alpha, beta, gamma):
    """ 
    Returns the lattice vectors from three lengths and angles 
    
    Parameters
    ----------
    a, b, c : float
        Lengths of lattice vectors [Angstroms]
    alpha, beta, gamma : float
        Angles of lattice vectors [degrees]
    
    Returns
    -------
    a1, a2 a3 : `~numpy.ndarray`, shape (3,)
        Lattice vectors
    """
    alpha, beta, gamma = map(radians, (alpha, beta, gamma))

    a1 = a*e1
    a2 = b * (cos(gamma)*e1 + sin(gamma)*e2)

    # Determine a3 = c1 *e1 + c2 * e2 + c3 * e3
    c1 = cos(beta)
    c2 = cos(alpha)/sin(gamma) - cos(beta)/tan(gamma)
    try:
        c3 = sqrt(1 - c1**2 - c2**2)    #
    except ValueError:
        raise ValueError('Invalid lattice parameters')

    return a1, a2, c*(c1*e1 + c2*e2 + c3*e3)

class Lattice(Transformable):
    """
    Abstract class handling Lattice information and manipulations.
    """

    def __init__(self, lattice_vectors = np.eye(3), **kwargs):
        """
        Parameters
        ----------
        lattice_vectors: iterable of `~numpy.ndarray`, shape (3,), optional
            Lattice vectors. Default is a cartesian lattice.
        """
        super().__init__(**kwargs)
        a1, a2, a3 = lattice_vectors
        self.a1 = np.asarray(a1, dtype = np.float) 
        self.a2 = np.asarray(a2, dtype = np.float) 
        self.a3 = np.asarray(a3, dtype = np.float)
    
    def __repr__(self):
        return '< Lattice object. a1 : {} \n, a2 : {} \n, a3 : {}>'.format(self.a1, self.a2, self.a3)
    
    @classmethod
    def from_parameters(cls, a, b, c, alpha, beta, gamma):
        """ 
        Create a lattice instance from three lengths and angles.

        Parameters
        ----------
        a, b, c : floats
            Lattice vectors lengths [Angs]
        alpha, beta, gamma : floats
            Angles between lattice vectors [deg]
        """
        return cls(lattice_vectors = lattice_vectors_from_parameters(a, b, c, alpha, beta, gamma))
    
    @property
    def lattice_parameters(self):
        """ Lattice parameters as three lengths and three angles. """
        a, b, c = norm(self.a1), norm(self.a2), norm(self.a3)
        alpha = np.arccos(np.vdot(self.a2, self.a3)/(b*c))
        beta = np.arccos(np.vdot(self.a1, self.a3)/(a*c))
        gamma = np.arccos(np.vdot(self.a1, self.a2)/(a*b))
        return a, b, c, np.rad2deg(alpha), np.rad2deg(beta), np.rad2deg(gamma)
    
    @property
    def volume(self):
        """ Lattice cell volume in angstroms cubed """
        return np.dot(self.a1, np.cross(self.a2, self.a3))
    
    @property
    def lattice_vectors(self):
        """ Iterable of lattice vectors """
        return self.a1, self.a2, self.a3
    
    @lattice_vectors.setter
    def lattice_vectors(self, vectors):
        self.a1, self.a2, self.a3 = vectors
    
    @property
    def reciprocal(self):
        """ Reciprocal lattice """
        return Lattice(lattice_vectors = self.reciprocal_vectors)

    @property
    def reciprocal_vectors(self):
        cell_volume = self.volume
        b1 = 2*np.pi*np.cross(self.a2, self.a3)/cell_volume
        b2 = 2*np.pi*np.cross(self.a3, self.a1)/cell_volume
        b3 = 2*np.pi*np.cross(self.a1, self.a2)/cell_volume
        return b1, b2, b3

    def transform(self, *matrices):
        """
        Transforms the real space coordinates according to a matrix.
        
        Parameters
        ----------
        matrices : ndarrays, shape {(3,3), (4,4)}
            Transformation matrices.
        """
        # Transform lattice vectors 
        for matrix in matrices:
            self.a1 = transform(matrix, self.a1)
            self.a2 = transform(matrix, self.a2)
            self.a3 = transform(matrix, self.a3)
        
        # Transform items
        super().transform(*matrices)
