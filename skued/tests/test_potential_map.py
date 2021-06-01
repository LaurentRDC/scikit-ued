# -*- coding: utf-8 -*-
from itertools import combinations_with_replacement

import numpy as np

from crystals import Crystal
from skued import (
    plane_mesh,
    potential_map,
    potential_synthesis,
    powdersim,
    structure_factor,
)
from skued.utils import suppress_warnings
import pytest


def test_potential_map_shape():
    """Test that potential_map returns a map with the same shape as the mesh"""
    crystal = Crystal.from_database("vo2-rutile")
    q = np.linspace(1, 5, 256)
    I = powdersim(crystal, q)

    aR1, aR2, aR3 = crystal.lattice_vectors
    extent = np.arange(0, 5, 0.1)

    with suppress_warnings():
        plane = plane_mesh(aR3, aR1 + aR2, x1=extent)

    potmap = potential_map(q, I, crystal, plane)

    xx, yy, zz = plane
    for arr in plane:
        assert potmap.shape == arr.shape


def test_potential_map_positive_intensity():
    """Test that potential_map raises an error if diffraction intensity is not positive"""
    crystal = Crystal.from_database("vo2-rutile")
    q = np.linspace(1, 5, 256)
    I = powdersim(crystal, q)

    aR1, aR2, aR3 = crystal.lattice_vectors
    extent = np.arange(0, 5, 0.1)
    with suppress_warnings():
        plane = plane_mesh(aR3, aR1 + aR2, x1=extent)

    I[0] = -1
    with pytest.raises(ValueError):
        potmap = potential_map(q, I, crystal, plane)


def test_potential_map_trivial():
    """Test that potential_map calculated from zero intensity is zero everywhere"""
    crystal = Crystal.from_database("vo2-rutile")
    q = np.linspace(1, 5, 256)
    I = powdersim(crystal, q)

    aR1, aR2, aR3 = crystal.lattice_vectors
    extent = np.arange(0, 10, 0.1)

    with suppress_warnings():
        plane = plane_mesh(aR3, aR1 + aR2, x1=extent)

    potmap = potential_map(q, np.zeros_like(I), crystal, plane)

    assert np.allclose(potmap, 0)


def test_potential_synthesis_trivial():
    """Test that potential_synthesis calculated from zero intensity is zero everywhere"""
    crystal = Crystal.from_database("C")
    reflections = list(combinations_with_replacement(range(-3, 4), 3))
    intensities = [
        np.abs(structure_factor(crystal, *reflection)) ** 2
        for reflection in reflections
    ]

    aR1, aR2, aR3 = crystal.lattice_vectors
    extent = np.arange(0, 10, 0.1)

    with suppress_warnings():
        plane = plane_mesh(aR3, aR1 + aR2, x1=extent)

    potmap = potential_synthesis(
        reflections, np.zeros_like(intensities), crystal, plane
    )

    assert np.allclose(potmap, 0)


def test_potential_synthesis_positive_intensity():
    """Test that potential_synthesis raises an error if diffraction intensity is not positive"""
    crystal = Crystal.from_database("C")
    reflections = list(combinations_with_replacement(range(-3, 4), 3))
    intensities = [
        np.abs(structure_factor(crystal, *reflection)) ** 2
        for reflection in reflections
    ]

    aR1, aR2, aR3 = crystal.lattice_vectors
    extent = np.arange(0, 5, 0.1)
    with suppress_warnings():
        plane = plane_mesh(aR3, aR1 + aR2, x1=extent)

    intensities[0] = intensities[0] * -1
    with pytest.raises(ValueError):
        potmap = potential_synthesis(reflections, intensities, crystal, plane)


def test_potential_synthesis_shape():
    """Test that potential_synthesis returns a map with the same shape as the mesh"""
    crystal = Crystal.from_database("C")
    reflections = list(combinations_with_replacement(range(-3, 4), 3))
    intensities = [
        np.abs(structure_factor(crystal, *reflection)) ** 2
        for reflection in reflections
    ]

    aR1, aR2, aR3 = crystal.lattice_vectors
    extent = np.arange(0, 5, 0.1)

    with suppress_warnings():
        plane = plane_mesh(aR3, aR1 + aR2, x1=extent)

    potmap = potential_synthesis(reflections, intensities, crystal, plane)

    xx, yy, zz = plane
    for arr in plane:
        assert potmap.shape == arr.shape
