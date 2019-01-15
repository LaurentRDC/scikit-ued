# -*- coding: utf-8 -*-
"""
Linear algebra operations and helpers.

Inspired by Christoph Gohlke's transformation.py <http://www.lfd.uci.edu/~gohlke/>
"""

import math

import numpy as np

# standard basis
e1, e2, e3 = np.eye(3)


def affine_map(array):
    """
    Extends 3x3 transform matrices to 4x4, i.e. general affine transforms.
    
    Parameters
    ----------
    array : ndarray, shape {(3,3), (4,4)}
        Transformation matrix. If shape = (4,4), returned intact.
    
    Returns
    -------
    extended : ndarray, shape (4,4)
        Extended array

	Raises
	------
	ValueError : If the transformation matrix is neither 3x3 or 4x4
    """
    if array.shape == (4, 4):  # Already the right shape
        return array
    elif array.shape == (3, 3):
        extended_matrix = np.zeros(shape=(4, 4), dtype=array.dtype)
        extended_matrix[-1, -1] = 1
        extended_matrix[:3, :3] = array
        return extended_matrix
    else:
        raise ValueError(
            "Array shape not 3x3 or 4x4, and thus is not a transformation matrix."
        )


def transform(matrix, array):
    """
	Applies a matrix transform on an array.

    Parameters
    ----------
    matrix : ndarray, shape {(3,3), (4,4)}
        Transformation matrix.
    array : ndarray, shape {(3,), (3,3), (4,4)}
        Array to be transformed. Either a 1x3 vector, or a transformation
        matrix in 3x3 or 4x4 shape.
    
    Returns
    -------
    transformed : ndarray
        Transformed array, either a 1D vector or a 4x4 transformation matrix

	Raises
	------
	ValueError : If the transformation matrix is neither 3x3 or 4x4
    """
    if matrix.shape not in [(3, 3), (4, 4)]:
        raise ValueError(
            "Input matrix is neither a 3x3 or 4x4 matrix, but \
						  rather of shape {}.".format(
                matrix.shape
            )
        )

    matrix = affine_map(matrix)
    # Case of a vector (e.g. position vector):
    if array.ndim == 1:
        extended_vector = np.array([0, 0, 0, 1], dtype=array.dtype)
        extended_vector[:3] = array
        return np.dot(matrix, extended_vector)[:3]
    else:
        array = affine_map(array)
    return np.dot(matrix, array)


def translation_matrix(direction):
    """
	Return matrix to translate by direction vector.

	Parameters
	----------
	direction : array_like, shape (3,)

	Returns
	-------
	translation : `~numpy.ndarray`, shape (4,4)
		4x4 translation matrix.
    """
    matrix = np.eye(4)
    matrix[:3, 3] = np.asarray(direction)[:3]
    return matrix


def change_of_basis(basis1, basis2=(e1, e2, e3)):
    """
	Returns the matrix that goes from one basis to the other.
    
	Parameters
	----------
	basis1 : list of array_like, shape (3,)
		First basis
	basis2 : list of array_like, shape (3,), optional
		Second basis. By default, this is the standard basis

	Returns
	-------
	cob : `~numpy.ndarray`, shape (3,3)
		Change-of-basis matrix that, applied to `basis`, will
		return `basis2`.
	"""
    # Calculate the transform that goes from basis 1 to standard basis
    basis1 = [np.asarray(vector).reshape(3, 1) for vector in basis1]
    basis1_to_standard = np.hstack(tuple(basis1))

    # Calculate the transform that goes from standard basis to basis2
    basis2 = [np.asarray(vector).reshape(3, 1) for vector in basis2]
    standard_to_basis2 = np.linalg.inv(np.hstack(tuple(basis2)))

    return np.dot(standard_to_basis2, basis1_to_standard)


def is_basis(basis):
    """ 
	Returns true if the set of vectors forms a basis. This is done by checking
	whether basis vectors are independent via an eigenvalue calculation.
    
	Parameters
	----------
	basis : list of array-like, shape (3,)

	Returns
	-------
	out : bool
		Whether or not the basis is valid.
	"""
    return 0 not in np.linalg.eigvals(np.asarray(basis))


def is_rotation_matrix(matrix):
    """
    Checks whether a matrix is orthogonal with unit determinant (1 or -1), properties
    of rotation matrices.
    
    Parameters
    ----------
    matrix : ndarray, shape {(3,3), (4,4)}
        Rotation matrix candidate. If (4,4) matrix is provided,
        only the top-left block matrix of (3,) is checked
       
    Returns
    -------
    result : bool
        If True, input could be a rotation matrix.
    """
    # TODO: is this necessary? should a composite transformation
    # 		of translation and rotation return True?
    # if matrix.shape == (4,4):
    #    matrix = matrix[:3,:3]

    is_orthogonal = np.allclose(np.linalg.inv(matrix), np.transpose(matrix))
    unit_determinant = np.allclose(abs(np.linalg.det(matrix)), 1)
    return is_orthogonal and unit_determinant


def rotation_matrix(angle, axis=(0, 0, 1)):
    """ 
	Return matrix to rotate about axis defined by direction around the origin [0,0,0].
	To combine rotation and translations, see http://www.euclideanspace.com/maths/geometry/affine/matrix4x4/index.htm
    
	Parameters
	----------
	angle : float
		Rotation angle [rad]
	axis : array-like of length 3
		Axis about which to rotate
    
	Returns
	-------
	matrix : `~numpy.ndarray`, shape (3,3)
		Rotation matrix.

	See also
	--------
	translation_rotation_matrix
	"""
    sina, cosa = math.sin(angle), math.cos(angle)

    # Make sure direction is a numpy vector of unit length
    direction = np.asarray(axis)
    direction = direction / np.linalg.norm(direction)

    # rotation matrix around unit vector
    R = np.diag([cosa, cosa, cosa])
    R += np.outer(direction, direction) * (1.0 - cosa)
    direction *= sina
    R += np.array(
        [
            [0.0, -direction[2], direction[1]],
            [direction[2], 0.0, -direction[0]],
            [-direction[1], direction[0], 0.0],
        ]
    )

    return R


def translation_rotation_matrix(angle, axis, translation):
    """
	Returns a 4x4 matrix that includes a rotation and a translation.

	Parameters
	----------
	angle : float
		Rotation angle [rad]
	axis : array-like of length 3
		Axis about which to rotate
	translation : array_like, shape (3,)
		Translation vector

	Returns
	-------
	matrix : `~numpy.ndarray`, shape (4,4)
		Affine transform matrix.
	"""
    rmat = affine_map(rotation_matrix(angle=angle, axis=axis))
    rmat[:3, 3] = np.asarray(translation)
    return rmat


def change_basis_mesh(xx, yy, zz, basis1, basis2):
    """
    Changes the basis of meshgrid arrays.

    Parameters
    ----------
    xx, yy, zz : ndarrays
        Arrays of equal shape, such as produced by numpy.meshgrid.
    basis1 : list of ndarrays, shape(3,)
        Basis of the mesh
    basis2 : list of ndarrays, shape(3,)
        Basis in which to express the mesh
    
    Returns
    -------
    XX, YY, ZZ : `~numpy.ndarray`
    """
    # Build coordinate array row-wise
    changed = np.empty(shape=(3, xx.size), dtype=np.float)
    linearized = np.empty(shape=(3, xx.size), dtype=np.float)
    linearized[0, :] = xx.ravel()
    linearized[1, :] = yy.ravel()
    linearized[2, :] = zz.ravel()

    # Change the basis at each row
    COB = change_of_basis(basis1, basis2)
    np.dot(COB, linearized, out=changed)
    return (
        changed[0, :].reshape(xx.shape),
        changed[1, :].reshape(yy.shape),
        changed[2, :].reshape(zz.shape),
    )


def minimum_image_distance(xx, yy, zz, lattice):
    """
	Returns a periodic array according to the minimum image convention.

	Parameters
	----------
	xx, yy, zz : ndarrays
		Arrays of equal shape, such as produced by numpy.meshgrid.
	lattice : list of ndarrays, shape(3,)
		Basis of the mesh
    
	Returns
	-------
	r : `~numpy.ndarray`
		Minimum image distance over the lattice
	"""
    COB = change_of_basis(np.eye(3), lattice)
    linearized = np.empty(shape=(3, xx.size), dtype=np.float)  # In the standard basis
    ulinearized = np.empty_like(linearized)  # In the unitcell basis

    linearized[0, :] = xx.ravel()
    linearized[1, :] = yy.ravel()
    linearized[2, :] = zz.ravel()

    # Go to unitcell basis, where the cell is cubic of side length 1
    np.dot(COB, linearized, out=ulinearized)
    ulinearized -= np.rint(ulinearized)
    np.dot(np.linalg.inv(COB), ulinearized, out=linearized)

    return np.reshape(np.linalg.norm(linearized, axis=0), xx.shape)
