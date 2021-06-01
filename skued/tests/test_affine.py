# -*- coding: utf-8 -*-

from itertools import product

import numpy as np

from skued import affine as tr
import pytest

np.random.seed(23)


def test_affine_map_exception_raised():
    """Tests that affine_map() raises an exceptions for incorrect input"""
    with pytest.raises(ValueError):

        m = tr.affine_map(np.eye(5))  # Matrix too large
        m = tr.affine_map(np.eye(2))  # Matrix too small


def test_affine_map_output_shape():
    """Tests that affine_map() always returns a 4x4 array"""
    in3, in4 = np.eye(3), np.eye(4)

    assert tr.affine_map(in3).shape == (4, 4)
    assert tr.affine_map(in4).shape == (4, 4)


def test_affine_map_fill():
    """Tests that affine_map() fills with zeros (and one on diagonal) for
    an input of shape (3,3)"""
    arr = np.random.random(size=(3, 3))
    extended = tr.affine_map(arr)
    assert np.allclose(extended[:, -1], [0, 0, 0, 1])
    assert np.allclose(extended[-1, :], [0, 0, 0, 1])


def test_translation_random():
    """Tests that translation_matrix() has the same effect on a point
    as directly translating."""
    pt = np.random.random((3,))
    translation = np.random.random((3,))

    mat = tr.translation_matrix(translation)
    transformed = tr.transform(mat, pt)

    assert np.allclose(pt + translation, transformed)


def test_transform_vector():
    """Tests that transform() on vector returns a vector"""
    vec = np.random.random((3,))

    # Test for vector of length 3
    mat3, mat4 = np.eye(3), np.eye(4)
    for mat in [mat3, mat4]:
        assert vec.shape == tr.transform(mat, vec).shape


def test_transform_matrix():
    """Tests that transform() on matrix returns a 4x4 matrix, always"""
    arr3, arr4 = np.random.random((3, 3)), np.random.random((4, 4))
    mat3, mat4 = np.eye(3), np.eye(4)

    for mat, arr in product([mat3, mat4], [arr3, arr4]):
        assert tr.transform(mat, arr).shape == (4, 4)


def test_change_of_basis_trivial():
    """Test that change_of_basis() returns the identity operator
    for a trivial change of basis"""
    cob = tr.change_of_basis(np.eye(3), np.eye(3))
    assert np.allclose(cob, np.eye(3))


def test_change_of_basis():
    """Test that change_of_basis() returns a correct change-of-basis matrix"""

    b1 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])

    # Generate a new basis as a rotation of pi/3 around z-axis
    b2 = np.dot(tr.rotation_matrix(np.pi / 3, axis=[0, 0, 1]), b1)

    cob = tr.change_of_basis(b1, b2)
    assert np.allclose(np.dot(cob, b1), b2)


def test_is_basis():
    """Test that is_basis() correctly identifies that a
        basis of zeros is not a basis and that the standard
    basis is a basis"""
    assert not tr.is_basis(np.zeros((3, 3)))
    assert tr.is_basis(np.eye(3))


def test_is_rotation_matrix_trivial():
    """test that the identity matrix is a rotation matrix"""
    assert tr.is_rotation_matrix(np.eye(3))
    assert tr.is_rotation_matrix(np.eye(4))


def test_from_rotation_matrix():
    """test that the rotated identity operator is
    a rotation matrix"""
    assert tr.is_rotation_matrix(tr.rotation_matrix(np.pi / 3, axis=[0, 0, 1]))


def test_rotation_matrix_random():
    """Test that rotation_matrix() returns a valid rotation matrix
    for random axes and angles"""
    axis = np.random.random((3,))
    angle = np.random.random()

    mat = tr.rotation_matrix(angle, axis)
    assert tr.is_rotation_matrix(mat)


def test_translation_rotation_matrix_trivial():
    """Test that a translation_rotation_matrix() reduces to rotation_matrix()
    for zero translation"""

    axis = np.random.random((3,))
    angle = np.random.random()

    mat = tr.translation_rotation_matrix(angle, axis, translation=[0, 0, 0])
    assert tr.is_rotation_matrix(mat)


def test_translation_rotation_matrix_random():
    """Test that translation_rotation_matrix() produces a matrix that correctly
    transforms a random point"""
    pt = np.random.random((3,))

    axis = np.random.random((3,))
    angle = np.random.random()
    translation = np.random.random((3,))

    trmat = tr.translation_rotation_matrix(angle, axis, translation=translation)
    v1 = tr.transform(trmat, pt)  # translated rotated point

    # Transform the point once operator at a time
    v2 = tr.transform(tr.rotation_matrix(angle, axis), pt)
    v2 += translation

    assert np.allclose(v1, v2)


def test_change_basis_mesh_trivial_basis_change():
    """
    Test the change of basis from standard basis to standard basis.
    """
    extent = np.linspace(0, 10, 10, dtype=int)
    xx, yy, zz = np.meshgrid(extent, extent, extent)

    XX, YY, ZZ = tr.change_basis_mesh(
        xx=xx, yy=yy, zz=zz, basis1=np.eye(3), basis2=np.eye(3)
    )
    assert np.allclose(xx, XX)
    assert np.allclose(yy, YY)
    assert np.allclose(zz, ZZ)


def test_change_basis_mesh_coordinate_swap():
    """
    Tests the change of basis from (x, y, z) to (x, z, y)
    """
    extent = np.linspace(0, 10, 10, dtype=int)
    xx, yy, zz = np.meshgrid(extent, extent, extent)

    e1, e2, e3 = np.eye(3)
    swapped_basis = [e1, e3, e2]

    XX, YY, ZZ = tr.change_basis_mesh(
        xx=xx, yy=yy, zz=zz, basis1=np.eye(3), basis2=swapped_basis
    )
    assert np.allclose(xx, XX)
    assert np.allclose(yy, ZZ)
    assert np.allclose(zz, YY)


def test_change_basis_mesh_scaling():
    extent = np.linspace(0, 10, 10, dtype=int)
    xx, yy, zz = np.meshgrid(extent, extent, extent)

    e1, e2, e3 = np.eye(3)

    scaled_basis = [0.5 * e1, 0.5 * e2, 0.5 * e3]

    XX, YY, ZZ = tr.change_basis_mesh(
        xx=xx, yy=yy, zz=zz, basis1=np.eye(3), basis2=scaled_basis
    )
    assert np.allclose(2 * xx, XX)
    assert np.allclose(2 * yy, YY)
    assert np.allclose(2 * zz, ZZ)
