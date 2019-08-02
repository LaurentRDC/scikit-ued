# -*- coding: utf-8 -*-
"""
Modelling phonons
"""

import numpy as np

from crystals.affine import change_of_basis


def rowdot(arr, brr):
    """ Row-wise dot product """
    # This is much, much faster than np.inner for some reason.
    return np.einsum("ij,ij->i", arr, brr).reshape((-1, 1))


def mapply(matrix, table):
    """ Apply a matrix transformation to a table where every row is considered a vector. """
    return np.transpose(matrix @ table.T)


def unique_by_row_idx(arr):
    """ Return the indices of the unique rows in arr """
    # np.unique can be rather slow when checking for uniqueness in rows
    # The optimization below results in ~3X faster performance
    #
    # A discussion of why is located here:
    #   https://github.com/numpy/numpy/issues/11136
    arr = np.ascontiguousarray(arr)
    arr_row_view = arr.view("|S%d" % (arr.itemsize * arr.shape[1]))
    _, unique_row_indices = np.unique(arr_row_view, return_index=True)
    return unique_row_indices


def unique_by_rows(*arrs):
    """ Filter arrays by the unique rows of the first array """
    unique_row_indices = unique_by_row_idx(arrs[0])
    return tuple(arr[unique_row_indices] for arr in arrs)


def roughly_unique_by_rows(*arrs, decimals, axis=0):
    """ Apply uniqueness rules on an array, based on lower precision. """
    rough = np.copy(arrs[0])
    np.around(rough, decimals=decimals, out=rough)
    return unique_by_rows(rough, *arrs[1:])


def tile_over_rows(*arrs):
    """ Tile arrays over rows time until all arrays have the same number of rows as the first array"""
    nrows = arrs[0].shape[0]

    arrays = [arrs[0]]
    for array in arrs[1:]:
        missing_reps = int(nrows / array.shape[0])
        reps_tuple = [1] * array.ndim
        reps_tuple[0] = missing_reps
        newarr = np.tile(array, reps=tuple(reps_tuple))
        arrays.append(newarr)

    return tuple(arrays)


class PhononMode(object):
    """
    Modelling of a phonon mode, in reciprocal-space.

    Parameters
    ----------
    name : str
        Mode name, e.g. "LA".
    q_points : ndarray, shape (N, 3)
        Table of reciprocal space vectors where the mode is defined.
    frequencies : ndarray, shape (N, 1)
        Mode frequencies at every point in ``k_points`` [Hz]
    polarizations : ndarray, shape (N, 3, natoms), dtype complex
        Complex mode polarization PER ATOM.
    crystal : crystals.Crystal instance
        Associated crystal
    hkls : ndarray, shape (N, 3), dtype int, optional
        Nearest Bragg-peak associated with each row in k_points. 
        Default is the (000) reflection only.
    """

    def __init__(self, name, q_points, frequencies, polarizations, crystal, hkls=None):
        if hkls is None:
            hkls = np.zeros_like(q_points)

        self.name = name
        self.q_points = q_points
        self.frequencies = frequencies
        self.polarizations = polarizations
        self.crystal = crystal
        self.hkls = hkls

    def k_points(self):
        """ Determine the unique k-points in this mode. """
        from_miller = change_of_basis(
            np.array(self.crystal.reciprocal_vectors), np.eye(3)
        )
        bragg = mapply(from_miller, self.hkls)
        return self.q_points - bragg


def symmetrize(mode):
    """
    Extend mode information by symmetrization. 

    The input is assumed to be a Mode representing information in 
    a single Brillouin zone (000), unsymmetrized.
    """
    assert np.allclose(mode.hkls, np.zeros_like(mode.q_points))

    k_points_frac = mode.q_points  # called k_points because only one brillouin zone
    polarizations = mode.polarizations
    frequencies = mode.frequencies
    cryst = mode.crystal

    # Conversion matrices
    to_fractional = change_of_basis(np.eye(3), np.array(cryst.reciprocal_vectors))
    from_fractional = np.linalg.inv(to_fractional)

    # Apply symmetry information and tile arrays to same shape
    k_points_frac, polarizations = apply_symops(
        k_points_frac, polarizations, crystal=cryst
    )
    k_points_frac, polarizations, frequencies = tile_over_rows(
        k_points_frac, polarizations, frequencies
    )

    # Change of basis to inverse angstroms
    # Called k_points because still inside a single Brillouin zone.
    k_points = mapply(from_fractional, k_points_frac)

    return Mode(
        name=mode.name,
        q_points=k_points,
        frequencies=frequencies,
        polarizations=polarizations,
        crystal=cryst,
        hkls=mode.hkls,
    )


def extend_bragg(mode, reflections):
    """
    Expand mode information so that it covers the entire detector range. 
    
    The input is assumed to be a Mode representing information in a single Brillouin zone (000).
    Symmetry operations are applied, as well as some filtering.
    """
    q_points = mode.q_points
    polarizations = mode.polarizations
    frequencies = mode.frequencies

    over_reflections, hkls = list(), list()
    astar, bstar, cstar = mode.crystal.reciprocal_vectors
    for (h, k, l) in reflections:
        H = h * astar + k * bstar + l * cstar
        over_reflections.append(q_points + H[None, :])

        # Quick way to copy a chunk of h, k, l row-wise
        hkls.append(
            np.zeros_like(q_points) + np.array([h, k, l], dtype=np.int)[None, :]
        )

    # Important to distinguish
    q_points = np.vstack(over_reflections)
    hkls = np.vstack(hkls)
    q_points, polarizations, frequencies, hkls = tile_over_rows(
        q_points, polarizations, frequencies, hkls
    )

    return Mode(
        name=mode.name,
        q_points=q_points,
        frequencies=frequencies,
        polarizations=polarizations,
        crystal=mode.crystal,
        hkls=hkls,
    )


def decimate(mode, decimals=2):
    """
    Decimate the information contained in modes based on similarity
    i.e. q-points that are roughly the same will be purged.
    """
    # Filter qs to a slightly lower precision
    # Because of QE rounding error + symmetry operations,
    # lots of duplicated points...
    q_points, polarizations, frequencies, hkls = roughly_unique_by_rows(
        mode.q_points,
        mode.polarizations,
        mode.frequencies,
        mode.hkls,
        decimals=decimals,
    )

    return Mode(
        name=mode.name,
        q_points=q_points,
        polarizations=polarizations,
        frequencies=frequencies,
        hkls=hkls,
        crystal=mode.crystal,
    )
