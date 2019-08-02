# -*- coding: utf-8 -*-
"""
Utilities to handle density-functional perturbation theory results.
"""

from itertools import product

import numpy as np


def tensor_distance(a, b):
    """
    Calculate the row-wise distance between tensors

    Parameters
    ----------
    a, b : ndarray
        Tensors representing phonon polarization vectors.
    """
    nmodes = a.shape[0]

    distance = np.zeros((nmodes, nmodes), dtype=np.float)
    for i, j in product(range(nmodes), repeat=2):
        distance[i, j] = np.linalg.norm(a[i] - b[j])

    return distance


def estimate_band_connection(prev_eigvecs, eigvecs, prev_freqs, freqs, prev_band_order):
    """
    A function to order the phonon eigenvectors 
    """
    pols_distance = tensor_distance(a=prev_eigvecs, b=eigvecs)
    freq_distance = np.abs(prev_freqs[:, None] - freqs[None, :]) ** 2

    nmodes = len(prev_band_order)
    connection_order = []

    # Loop over the overlap of the eigenvectors
    for mode_index in range(nmodes):
        best_match_index = np.nanargmin(pols_distance[mode_index, :])
        best_match_freq_index = np.nanargmin(freq_distance[mode_index, :])

        # This mode has been connected, and therefore should not be considered
        # We cancel the column associated to it in so that it will never
        # be picked again.
        pols_distance[:, best_match_index] = np.nan
        freq_distance[:, best_match_index] = np.nan

        connection_order.append(best_match_index)

    return [connection_order[x] for x in prev_band_order]


# TODO: add test
def band_order(frequencies, polarizations, crystal):
    """
    Determine the band order based on q-points and frequencies

    Parameters
    ----------
    frequencies : ndarray, shape (N, M)
        Frequencies [THz] for each of the `M` phonon modes.
    polarizations : ndarray, shape (N, M, natoms, 3)
        Complex polarization vectors of the `M` phonon modes and `natoms` atoms.
    
    Returns
    -------
    frequencies : ndarray, shape (N, M)
        Ordered frequencies [THz] for each of the `M` phonon modes. Each column
        should represent a single band.
    polarizations : ndarray, shape (N, M, natoms)
        Ordered polarization vectors of the `M` phonon modes and `natoms` atoms.
    """
    pols = np.copy(polarizations)
    freqs = np.copy(frequencies)

    # Scale polarization vectors by atomic mass
    # Displacement of a large atom should count more than the displacement
    # of a smaller one
    # TODO: this does not seem necessary...
    atoms = sorted(crystal, key=lambda a: a.tag)
    sqrt_atm_masses = np.sqrt([atm.mass for atm in atoms])
    pols[:, :] *= sqrt_atm_masses[:, None]

    # Normalize polarization vectors, just in case
    pols /= np.linalg.norm(pols, axis=-1, keepdims=True)

    nqpoints, nphonons = frequencies.shape
    order = np.zeros((nqpoints, nphonons), dtype=np.int)
    order[0, :] = np.array(range(nphonons))

    for q_index in range(1, nqpoints):
        old_pols, new_pols = pols[q_index - 1], pols[q_index]
        old_freqs, new_freqs = freqs[q_index - 1], freqs[q_index]
        order[q_index] = estimate_band_connection(
            prev_eigvecs=old_pols,
            eigvecs=new_pols,
            prev_freqs=old_freqs,
            freqs=new_freqs,
            prev_band_order=order[q_index - 1],
        )

    for q_index in range(1, nqpoints):
        polsq = pols[q_index, :]
        freqsq = freqs[q_index, :]
        pols[q_index] = polsq[order[q_index]]
        freqs[q_index] = freqsq[order[q_index]]

    return freqs, pols


def apply_symops(kpoints, polarizations, crystal, symprec=1e-1):
    """
    Apply symmetry operations to polarizations vectors and q-points

    kpoints : ndarray, shape (N,3)
        Scattering vector within one Brillouin zone
    polarizations : ndarray, shape (N, natoms, 3)
        Complex polarization vectors. Every row is associated with the corresponding row
        in `kpoints`.
    crystal: crystals.Crystal
        Crystal object with the appropriate symmetry.
    symprec : float, optional
        Symmetry-determination precision.
    
    Returns
    -------
    """
    # Change of basis matrices allow to express
    # transformations in other bases
    to_reciprocal = change_of_basis(
        np.array(crystal.lattice_vectors), np.array(crystal.reciprocal_vectors)
    )
    from_reciprocal = np.linalg.inv(to_reciprocal)
    reciprocal = lambda m: to_reciprocal @ m @ from_reciprocal

    to_frac = change_of_basis(np.eye(3), np.array(crystal.lattice_vectors))
    from_frac = np.linalg.inv(to_frac)
    cartesian = lambda m: from_frac @ m @ to_frac

    transformed_k, transformed_p = list(), list()

    transformations = crystal.symmetry_operations(symprec=symprec)
    for rotation, _ in transformations:
        transformed_k.append(mapply(reciprocal(rotation), kpoints))

        # Transforming polarizations is a little more complex
        # because of the extra dimension.
        newpols = np.copy(polarizations)
        for atm_index in range(len(crystal)):
            newpols[:, atm_index, :] = mapply(
                cartesian(rotation), polarizations[:, atm_index, :]
            )

        transformed_p.append(newpols)

    return np.vstack(transformed_k), np.vstack(transformed_p)
