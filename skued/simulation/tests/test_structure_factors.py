# -*- coding: utf-8 -*-

from random import randint

import numpy as np

from crystals import Atom, Crystal
from skued import affe, structure_factor


def test_side_effects():
    nG = np.random.random(size=(16, 32))
    nG.setflags(write=False)  # if nG is written to, Exception is raised
    affe(Atom("He", coords=[0, 0, 0]), nG)


def test_out_shape():
    nG = np.random.random(size=(16, 32))
    eff = affe(Atom("He", coords=[0, 0, 0]), nG)
    assert eff.shape == nG.shape


def test_int():
    """Test that affe(int, ...) also works"""
    atomic_number = randint(1, 103)
    nG = np.random.random(size=(16, 32))

    from_int = affe(atomic_number, nG)
    from_atom = affe(Atom(atomic_number, [0, 0, 0]), nG)

    assert np.allclose(from_int, from_atom)


def test_str():
    """Test that affe(str, ...) also works"""
    # Try with Chlorine (Z = 17)
    atomic_number = 17
    nG = np.random.random(size=(16, 32))

    from_int = affe(atomic_number, nG)
    from_str = affe("Cl", nG)

    assert np.allclose(from_int, from_str)


def test_shape_and_dtype():
    """Test that output of structure_factor is same shape as input,
    and that the dtype is complex"""
    crystal = Crystal.from_database(next(iter(Crystal.builtins)))
    h, k, l = np.meshgrid([1, 2, 3], [1, 2, 3], [1, 2, 3])
    sf = structure_factor(crystal, h, k, l)

    assert sf.shape == h.shape
    assert sf.dtype == complex
