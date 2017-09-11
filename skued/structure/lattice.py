# -*- coding: utf-8 -*-
from math import sin, cos, tan, sqrt, radians
import numpy as np
from numpy.linalg import norm
from .. import transform, change_basis_mesh

e1, e2, e3 = np.eye(3) # Euclidian basis

# TODO: Introduce conventions on ordering a, b, c and angles
#       based on http://atztogo.github.io/spglib/definition.html
def lattice_vectors_from_parameters(a, b, c, alpha, beta, gamma):
    """ 
    Returns the lattice vectors from three lengths and angles.
    
    Parameters
    ----------
    a, b, c : float
        Lengths of lattice vectors [Angstroms]
    alpha, beta, gamma : float
        Angles of lattice vectors [degrees]. 
    
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

class Lattice(object):
    """
    Container class for lattice information and manipulations.

    Instances can also be create from the standard 'three lengths and angles'
    parameters via ``Lattice.from_parameters``:

    Parameters
    ----------
    lattice_vectors: iterable of `~numpy.ndarray`, shape (3,), optional
        Lattice vectors. Default is a cartesian lattice.
    """
    def __init__(self, lattice_vectors, **kwargs):
        a1, a2, a3 = lattice_vectors
        self.a1 = np.asarray(a1, dtype = np.float) 
        self.a2 = np.asarray(a2, dtype = np.float) 
        self.a3 = np.asarray(a3, dtype = np.float)
    
    def __repr__(self):
        return '< Lattice object. a1 : {} \n, a2 : {} \n, a3 : {}>'.format(self.a1, self.a2, self.a3)

    def __eq__(self, other):
        return np.allclose(self.lattice_vectors, other.lattice_vectors)

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
        """ Lattice parameters as three lengths [:math:`\AA`] and three angles [degrees]. """
        a, b, c = norm(self.a1), norm(self.a2), norm(self.a3)
        alpha = np.arccos(np.vdot(self.a2, self.a3)/(b*c))
        beta = np.arccos(np.vdot(self.a1, self.a3)/(a*c))
        gamma = np.arccos(np.vdot(self.a1, self.a2)/(a*b))
        return a, b, c, np.rad2deg(alpha), np.rad2deg(beta), np.rad2deg(gamma)
    
    @property
    def volume(self):
        """ Lattice cell volume [:math:`\^3`] """
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
        """
        Reciprocal lattice vectors, defined as:

        .. math::

            b_i = 2 \pi \\frac{a_j \\times a_k}{v}
        
        For :math:`v` the unit cell volume.
        """
        cell_volume = self.volume
        b1 = 2*np.pi*np.cross(self.a2, self.a3)/cell_volume
        b2 = 2*np.pi*np.cross(self.a3, self.a1)/cell_volume
        b3 = 2*np.pi*np.cross(self.a1, self.a2)/cell_volume
        return b1, b2, b3

    @property
    def periodicity(self):
        """ Crystal periodicity in x, y and z direction from the lattice constants.
        This is effectively a bounding cube for the unit cell, which is itself a unit cell. """
        e1, e2, e3 = np.eye(3)
        per_x = sum( (abs(np.vdot(e1,a)) for a in self.lattice_vectors) )
        per_y = sum( (abs(np.vdot(e2,a)) for a in self.lattice_vectors) )
        per_z = sum( (abs(np.vdot(e3,a)) for a in self.lattice_vectors) )
        return per_x, per_y, per_z

    def scattering_vector(self, h, k, l):
        """
        Scattering vector from Miller indices.

        Parameters
        ----------
        h, k, l : array_like
            Miller indices. 

        Returns
        -------
        Gx, Gy, Gz : `~numpy.ndarray`
            Components of the scattering vectors, of the same shape 
            as ``h``, ``k``, and ``l``.
        """
        h, k, l = np.atleast_1d(h, k, l)
        return change_basis_mesh(h, k, l, basis1 = self.reciprocal_vectors, basis2 = np.eye(3))

    def miller_indices(self, Gx, Gy, Gz):
        """
        Miller indices from scattering vector components.

        Parameters
        ----------
        Gx, Gy, Gz : `~numpy.ndarray`
            Scattering vector components, in :math:`A^{-1}`.
        
        Returns
        -------
        h, k, l : `~numpy.ndarray`
            Miller indices.
        """
        Gx, Gy, Gz = np.atleast_1d(Gx, Gy, Gz)
        return change_basis_mesh(Gx, Gy, Gz, basis1 = np.eye(3), basis2 = self.reciprocal_vectors)

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