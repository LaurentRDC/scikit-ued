# -*- coding: utf-8 -*-

from copy import deepcopy

import numpy as np

from crystals import Crystal
from skued import electrostatic, pelectrostatic


def test_return_shape():
    """Test that the return shape of pelectrostatic is the same as input arrays"""
    crystal = Crystal.from_database("C")
    xx, yy, zz = np.meshgrid(
        np.linspace(-10, 10, 16), np.linspace(-10, 10, 16), np.linspace(-10, 10, 16)
    )
    potential = electrostatic(crystal, xx, yy, zz)

    assert xx.shape == potential.shape


def test_side_effects():
    """Test that mesh arrays are not written to in pelectrostatic"""
    xx, yy, zz = np.meshgrid(
        np.linspace(-10, 10, 16), np.linspace(-10, 10, 16), np.linspace(-10, 10, 16)
    )

    xx.setflags(write=False)
    yy.setflags(write=False)
    zz.setflags(write=False)

    crystal = Crystal.from_database("C")
    potential = electrostatic(crystal, xx, yy, zz)


def test_return_shape():
    """Test that the return shape of pelectrostatic is the same as input arrays"""
    crystal = Crystal.from_database("C")
    xx, yy = np.meshgrid(np.linspace(-10, 10, 32), np.linspace(-10, 10, 32))
    potential = pelectrostatic(crystal, xx, yy)

    assert xx.shape == potential.shape


def test_side_effects():
    """Test that mesh arrays are not written to in pelectrostatic"""
    crystal = Crystal.from_database("C")
    xx, yy = np.meshgrid(np.linspace(-10, 10, 32), np.linspace(-10, 10, 32))
    xx.setflags(write=False)
    yy.setflags(write=False)
    potential = pelectrostatic(crystal, xx, yy)


def test_trivial():
    """Test that the projected electrostatic potential from an empty slice of crystal is zero"""
    crystal = Crystal.from_database("C")
    xx, yy = np.meshgrid(np.linspace(-10, 10, 32), np.linspace(-10, 10, 32))
    potential = pelectrostatic(crystal, xx, yy, bounds=(0, 1))
    assert np.allclose(potential, 0)
