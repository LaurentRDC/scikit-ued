"""
Linear algebra operations and helpers.

Inspired by Christoph Gohlke's transformation.py <http://www.lfd.uci.edu/~gohlke/>
"""

import numpy as np

#standard basis
e1,e2,e3 = np.eye(3)

def extend_matrix(array):
    """
    Extends 3x3 transform matrices to 4x4.
    
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
	ValueError
		If the transformation matrix is neither 3x3 or 4x4
    """
    if array.shape == (4,4):        #Already the right shape
        return array
    elif array.shape == (3,3):
        extended_matrix = np.zeros(shape = (4,4), dtype = array.dtype)
        extended_matrix[-1, -1] = 1
        extended_matrix[:3,:3] = array
        return extended_matrix
    else:
        raise ValueError('Array shape not 3x3 or 4x4, and thus is not a transformation matrix.')

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
	ValueError
		If the transformation matrix is neither 3x3 or 4x4
    """
	if matrix.shape not in [(3,3), (4,4)]:
		raise ValueError('Input matrix is neither a 3x3 or 4x4 matrix, but \
						  rather of shape {}.'.format(matrix.shape))

	# Case of a vector (e.g. position vector):
	if array.ndim == 1:

		if matrix.shape == (3,3):
			matrix = extend_matrix(matrix)
		extended_vector = np.array([0,0,0,1], dtype = array.dtype)
		extended_vector[:3] = array 
		return np.dot(matrix, extended_vector)[:3]

	# Transformation matrix is either 3x3 or 4x4. Extend everything to 4x4
	matrix, array = extend_matrix(matrix), extend_matrix(array)
	return np.dot(matrix, array)

def change_of_basis(basis1, basis2 = [e1,e2,e3]):
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
	cob : ndarray, shape (3,3)
		Change-of-basis matrix that, applied to `basis`, will
		return `basis2`.
	"""
	# Calculate the transform that goes from basis 1 to standard basis 
	basis1 = [np.asarray(vector).reshape(3,1) for vector in basis1]
	basis1_to_standard = np.hstack(tuple(basis1))  

	# Calculate the transform that goes from standard basis to basis2
	basis2 = [np.asarray(vector).reshape(3,1) for vector in basis2]
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
	"""
	return (0 not in np.linalg.eigvals(np.asarray(basis)))

def is_rotation_matrix(matrix):
    """
    Checks whether a matrix is orthogonal with unit determinant, properties
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
	#		of translation and rotation return True?
    #if matrix.shape == (4,4):
    #    matrix = matrix[:3,:3]
        
    is_orthogonal = np.allclose(np.linalg.inv(matrix), np.transpose(matrix))
    unit_determinant = np.allclose(np.linalg.det(matrix), 1)
    return is_orthogonal and unit_determinant

def rotation_matrix(angle, axis = [0,0,1]):
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
	matrix : ndarray, shape (3,3)
		Rotation matrix.

	See also
	--------
	translation_rotation_matrix
	"""
	sina, cosa = np.sin(angle), np.cos(angle)
	
	# Make sure direction is a numpy vector of unit length
	direction = np.asarray(axis)
	direction = direction/np.linalg.norm(direction)
	
	# rotation matrix around unit vector
	R = np.diag([cosa, cosa, cosa])
	R += np.outer(direction, direction) * (1.0 - cosa)
	direction *= sina
	R += np.array([[ 0.0,         -direction[2],  direction[1]],
					[ direction[2], 0.0,          -direction[0]],
					[-direction[1], direction[0],  0.0]])
	
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
	matrix : ndarray, shape (4,4)
		Affine transform matrix.
	"""
	rmat = extend_matrix(rotation_matrix(angle = angle, axis = axis))
	rmat[:3, 3] = np.asarray(translation)
	return rmat
	