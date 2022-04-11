# -*- coding: utf-8 -*-
"""
Electrostatic potential simulation
==================================
"""
from math import sqrt

import numpy as np
from numpy import pi
from scipy.special import k0 as bessel

from .. import minimum_image_distance
from .scattering_params import scattering_params

m = 9.109 * 10 ** (-31)  # in kg
a0 = 0.5291  # in Angs
e = 14.4  # Volt*Angstrom


def _electrostatic_atom(atom, r):
    try:
        _, a1, b1, a2, b2, a3, b3, c1, d1, c2, d2, c3, d3 = scattering_params[
            atom.atomic_number
        ]
    except KeyError:
        raise ValueError(
            f"Scattering information for element {atom.element} is unavailable."
        )

    sum1 = np.zeros_like(r, dtype=float)
    for a, b in zip((a1, a2, a3), (b1, b2, b3)):
        sum1 += a / r * np.exp(-2 * pi * r * np.sqrt(b))

    sum2 = np.zeros_like(r, float)
    for c, d in zip((c1, c2, c3), (d1, d2, d3)):
        sum2 += c * (d ** (-1.5)) * np.exp(-((r * pi) ** 2) / d)

    e = 14.4  # [Volt-Angstrom]
    a0 = 0.5291  # [Angs]
    return 2 * a0 * e * (pi**2 * sum1 + pi ** (2.5) * sum2)


def electrostatic(crystal, x, y, z):
    """
    Electrostatic potential from a crystal calculated on a real-space mesh,
    assuming an infinite crystal.

    Parameters
    ----------
    crystal : crystals.Crystal

    x, y, z : `~numpy.ndarray`
        Real space coordinates mesh.

    Returns
    -------
    potential : `~numpy.ndarray`, dtype float
        Linear superposition of atomic potential [V*Angs]

    See also
    --------
    pelectrostatic: projected electrostatic potential of an infinite crystal.
    """
    # TODO: split xx and yy into smalled non-repeating unit
    # TODO: multicore

    potential = np.zeros_like(x, dtype=float)
    r = np.zeros_like(x, dtype=float)
    for atom in crystal:
        ax, ay, az = atom.coords_cartesian
        r[:] = minimum_image_distance(
            x - ax, y - ay, z - az, lattice=crystal.lattice_vectors
        )
        potential += _electrostatic_atom(atom, r)

    # Due to sampling, x,y, and z might pass through the center of atoms
    # Replace np.inf by the next largest value
    m = potential[np.isfinite(potential)].max()
    potential[np.isinf(potential)] = m
    return potential


def _pelectrostatic_atom(atom, r):
    try:
        _, a1, b1, a2, b2, a3, b3, c1, d1, c2, d2, c3, d3 = scattering_params[
            atom.atomic_number
        ]
    except KeyError:
        raise ValueError(
            f"Scattering information for element {atom.element} is unavailable."
        )

    potential = np.zeros_like(r, dtype=float)
    for a, b, c, d in zip((a1, a2, a3), (b1, b2, b3), (c1, c2, c3), (d1, d2, d3)):
        potential += 2 * a * bessel(2 * pi * r * sqrt(b)) + (c / d) * np.exp(
            -((r * pi) ** 2) / d
        )

    return 2 * a0 * e * (pi**2) * potential


def pelectrostatic(crystal, x, y, bounds=None):
    """
    Projected electrostatic potential from a crystal calculated on a real-space mesh,
    assuming an infinite crystal in x and y. Projection axis is defined as the z-axis.
    To project the potential along a different axis, the crystal can be rotated with ``Crystal.transform``.

    Parameters
    ----------
    crystal : crystals.Crystal

    x, y:  `~numpy.ndarray`
        Real-space coordinates.
    bounds : iterable or None, optional
        Bounds of atom inclusion. Atoms with real-space z-position outside [ min(bounds), max(bounds) )
        are not counted in the computation.

    Returns
    -------
    potential : `~numpy.ndarray`, dtype float
        Linear superposition of electrostatic potential [V*Angs]

    See also
    --------
    electrostatic: three-dimensional electrostatic potential of an infinite crystal.
    """
    # TODO: split xx and yy into smalled non-repeating unit
    #       np.unique(np.mod(xx, per_x))
    # TODO: multicore

    if bounds:
        min_z, max_z = min(bounds), max(bounds)
        atoms = (
            atom for atom in iter(crystal) if min_z <= atom.coords_cartesian[2] < max_z
        )
    else:
        atoms = iter(crystal)

    potential = np.zeros_like(x, dtype=float)
    z = np.zeros_like(x)
    lattice = np.array(crystal.lattice_vectors)
    for atom in atoms:
        xa, ya, _ = atom.coords_cartesian
        r = minimum_image_distance(x - xa, y - ya, z, lattice=lattice)
        potential += _pelectrostatic_atom(atom, r)

    # Due to sampling, x,y, and z might pass through the center of atoms
    # Replace n.inf by the next largest value
    potential[np.isinf(potential)] = np.nan
    potential[np.isnan(potential)] = np.nanmax(potential)
    return potential
