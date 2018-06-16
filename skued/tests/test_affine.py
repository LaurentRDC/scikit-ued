# -*- coding: utf-8 -*-
from itertools import product
import numpy as np
from .. import affine as tr
import unittest

np.random.seed(23)

class TestAffineMap(unittest.TestCase):

	def test_exception_raised(self):
		""" Tests that affine_map() raises an exceptions for incorrect input """
		with self.assertRaises(ValueError):
			
			m = tr.affine_map(np.eye(5)) # Matrix too large
			m = tr.affine_map(np.eye(2)) # Matrix too small

	def test_output_shape(self):
		""" Tests that affine_map() always returns a 4x4 array """
		in3, in4 = np.eye(3), np.eye(4)

		self.assertSequenceEqual(tr.affine_map(in3).shape, (4,4))
		self.assertSequenceEqual(tr.affine_map(in4).shape, (4,4))

	def test_fill(self):
		""" Tests that affine_map() fills with zeros (and one on diagonal) for
		an input of shape (3,3) """
		arr = np.random.random(size = (3,3))
		extended = tr.affine_map(arr)
		self.assertTrue(np.allclose(extended[:, -1], [0,0,0,1]))
		self.assertTrue(np.allclose(extended[-1, :], [0,0,0,1]))

class TestTranslationMatrix(unittest.TestCase):
	
	def test_random(self):
		""" Tests that translation_matrix() has the same effect on a point
		as directly translating."""
		pt = np.random.random((3,))
		translation = np.random.random((3,))

		mat = tr.translation_matrix(translation)
		transformed = tr.transform(mat, pt)

		self.assertTrue(np.allclose(pt + translation, transformed))

class TestTransform(unittest.TestCase):

	def test_transform_vector(self):
		""" Tests that transform() on vector returns a vector """
		vec = np.random.random((3,))

		# Test for vector of length 3
		mat3, mat4 = np.eye(3), np.eye(4)
		for mat in [mat3, mat4]:
			self.assertSequenceEqual(vec.shape, tr.transform(mat, vec).shape)

	def test_transform_matrix(self):
		""" Tests that transform() on matrix returns a 4x4 matrix, always """
		arr3, arr4 = np.random.random((3,3)), np.random.random((4,4))
		mat3, mat4 = np.eye(3), np.eye(4)

		for mat, arr in product([mat3, mat4], [arr3, arr4]):
			self.assertSequenceEqual(tr.transform(mat, arr).shape, (4,4))

class TestChangeOfBasis(unittest.TestCase):

	def test_trivial(self):
		""" Test that change_of_basis() returns the identity operator
		for a trivial change of basis """
		cob = tr.change_of_basis(np.eye(3), np.eye(3))
		self.assertTrue(np.allclose(cob, np.eye(3)))

	def test_change_of_basis(self):
		""" Test that change_of_basis() returns a correct change-of-basis matrix"""

		b1 = np.array([[0,0,1],
					   [1,0,0], 
					   [0,1,0]])

		# Generate a new basis as a rotation of pi/3 around z-axis
		b2 = np.dot(tr.rotation_matrix(np.pi/3, axis = [0,0,1]), b1)

		cob = tr.change_of_basis(b1, b2)
		self.assertTrue(np.allclose(np.dot(cob, b1), b2))

class TestIsBasis(unittest.TestCase):

	def test_trivial(self):
		""" Test that is_basis() correctly identifies that a 
		basis of zeros is not a basis and that the standard
	    basis is a basis"""
		self.assertFalse(tr.is_basis(np.zeros((3,3))))
		self.assertTrue(tr.is_basis(np.eye(3)))

class TestIsRotationMatrix(unittest.TestCase):
	
	def test_trivial(self):
		""" test that the identity matrix is a rotation matrix """
		self.assertTrue(tr.is_rotation_matrix(np.eye(3)))
		self.assertTrue(tr.is_rotation_matrix(np.eye(4)))

	def test_from_rotation_matrix(self):
		""" test that the rotated identity operator is
		a rotation matrix"""
		self.assertTrue(tr.is_rotation_matrix(tr.rotation_matrix(np.pi/3, axis = [0,0,1])))

class TestRotationMatrix(unittest.TestCase):
	
	def test_random(self):
		""" Test that rotation_matrix() returns a valid rotation matrix
		for random axes and angles """
		axis = np.random.random((3,))
		angle = np.random.random()
		
		mat = tr.rotation_matrix(angle, axis)
		self.assertTrue(tr.is_rotation_matrix(mat))

class TestTranslationRotationMatrix(unittest.TestCase):
	
	def test_trivial(self):
		""" Test that a translation_rotation_matrix() reduces to rotation_matrix()
		for zero translation """

		axis = np.random.random((3,))
		angle = np.random.random()
		
		mat = tr.translation_rotation_matrix(angle, axis, translation = [0,0,0])
		self.assertTrue(tr.is_rotation_matrix(mat))

	def test_random(self):
		""" Test that translation_rotation_matrix() produces a matrix that correctly
		transforms a random point """
		pt = np.random.random((3,))

		axis = np.random.random((3,))
		angle = np.random.random()
		translation = np.random.random((3,))

		trmat = tr.translation_rotation_matrix(angle, axis, translation = translation)
		v1 = tr.transform(trmat, pt) #translated rotated point

		# Transform the point once operator at a time
		v2 = tr.transform(tr.rotation_matrix(angle, axis), pt)
		v2 += translation
		
		self.assertTrue(np.allclose(v1, v2))

class TestChangeBasisMesh(unittest.TestCase):
    """
    Tests related to the change_basis_mesh function
    """
        
    def test_trivial_basis_change(self):
        """
        Test the change of basis from standard basis to standard basis.
        """
        extent = np.linspace(0, 10, 10, dtype = np.int)
        xx, yy, zz = np.meshgrid(extent, extent, extent)

        XX, YY, ZZ = tr.change_basis_mesh(xx = xx, yy = yy, zz = zz, basis1 = np.eye(3), basis2 = np.eye(3))
        self.assertTrue(np.allclose(xx, XX))
        self.assertTrue(np.allclose(yy, YY))
        self.assertTrue(np.allclose(zz, ZZ))
    
    def test_coordinate_swap(self):
        """
        Tests the change of basis from (x, y, z) to (x, z, y)
        """
        extent = np.linspace(0, 10, 10, dtype = np.int)
        xx, yy, zz = np.meshgrid(extent, extent, extent)

        e1, e2, e3 = np.eye(3)
        swapped_basis = [e1, e3, e2]

        XX, YY, ZZ = tr.change_basis_mesh(xx = xx, yy = yy, zz = zz, basis1 = np.eye(3), basis2 = swapped_basis)
        self.assertTrue(np.allclose(xx, XX))
        self.assertTrue(np.allclose(yy, ZZ))
        self.assertTrue(np.allclose(zz, YY))
    
    def test_scaling(self):
        extent = np.linspace(0, 10, 10, dtype = np.int)
        xx, yy, zz = np.meshgrid(extent, extent, extent)

        e1, e2, e3 = np.eye(3)
        
        scaled_basis = [0.5*e1, 0.5*e2, 0.5*e3]

        XX, YY, ZZ = tr.change_basis_mesh(xx = xx, yy = yy, zz = zz, basis1 = np.eye(3), basis2 = scaled_basis)
        self.assertTrue(np.allclose(2*xx, XX))
        self.assertTrue(np.allclose(2*yy, YY))
        self.assertTrue(np.allclose(2*zz, ZZ))



if __name__ == '__main__':
	unittest.main()