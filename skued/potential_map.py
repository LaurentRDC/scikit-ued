# -*- coding: utf-8 -*-
"""
2D/3D Reconstruction of scattering potential
=========================================
"""

import numpy as np

from crystals.affine import change_basis_mesh

from .simulation import powdersim, structure_factor


# TODO: add tutorial
# TODO: use potential_synthesis inside potential_map
def potential_map(q, I, crystal, mesh):
    """
    Compute the electrostatic potential from powder diffraction data.

    Parameters
    ----------
    q : ndarray, shape (N,)
        Scattering vector norm (:math:`Å^{-1}`).
    I : ndarray, shape (N,)
        Experimental diffracted intensity.
    crystal : crystals.Crystal
        Crystal that gave rise to diffraction pattern `I`.
    mesh : 3-tuple ndarrays, ndim 2 or ndim 3
        Real-space mesh over which to calculate the scattering map.
        Format should be similar to the output of numpy.meshgrid.

    Returns
    -------
    out : ndarray, ndim 2 or ndim 3
        Electrostatic potential computed over the mesh.

    Raises
    ------
    ValueError: if intensity data is not strictly positive.

    Notes
    -----
    To compute the scattering map from a difference of intensities, note that scattering
    maps are linear in *structure factor* norm. Thus, to get the map of difference data :code:`I1 - I2`:

    .. math::

        I = (\sqrt{I_1} - \sqrt{I_2})^2

    References
    ----------
    .. [#] Otto et al., How optical excitation controls the structure and properties of vanadium dioxide.
           PNAS, vol. 116 issue 2, pp. 450-455 (2018). :DOI:`10.1073/pnas.1808414115`
    """
    if np.any(I < 0):
        raise ValueError("Diffracted intensity cannot physically be negative.")

    # We want to support 2D and 3D meshes, therefore
    # expand mesh until 4D (last dim is for loop over reflections)
    # Extra dimensions will be squeezed out later
    xx, yy, zz = mesh
    while xx.ndim < 4:
        xx, yy, zz = (
            np.expand_dims(xx, xx.ndim),
            np.expand_dims(yy, yy.ndim),
            np.expand_dims(zz, zz.ndim),
        )

    # Prepare reflections
    # G is reshaped so that it is perpendicular to xx, yy, zz to enables broadcasting
    reflections = np.vstack(tuple(crystal.bounded_reflections(q.max())))
    hs, ks, ls = np.hsplit(reflections, 3)
    SF = structure_factor(crystal, hs, ks, ls)

    # Extract structure factor with correction factors
    # Diffracted intensities add up linearly (NOT structure factors)
    qx, qy, qz = change_basis_mesh(
        hs, ks, ls, basis1=crystal.reciprocal_vectors, basis2=np.eye(3)
    )
    qx, qy, qz = (
        qx.reshape((1, 1, 1, -1)),
        qy.reshape((1, 1, 1, -1)),
        qz.reshape((1, 1, 1, -1)),
    )
    SF = SF.reshape((1, 1, 1, -1))
    q_theo = np.squeeze(np.sqrt(qx**2 + qy**2 + qz**2))
    theo_I = powdersim(crystal, q_theo)
    peak_mult_corr = np.abs(SF) ** 2 / theo_I
    exp_SF = np.sqrt(np.interp(q_theo, q, I) * peak_mult_corr)

    # Squeeze out extra dimensions (e.g. if mesh was 2D)
    potential_map = np.sum(
        exp_SF
        * np.real(np.exp(1j * np.angle(SF)))
        * np.cos(xx * qx + yy * qy + zz * qz),
        axis=3,
    )
    return np.squeeze(potential_map)


def potential_synthesis(reflections, intensities, crystal, mesh):
    """
    Synthesize the electrostatic potential from a list of experimental
    reflections and associated diffracted intensities. Diffraction phases are
    taken from a known structure

    Parameters
    ----------
    reflections : iterable of tuples
        Iterable of Miller indices as tuples (e.g. `[(0,1,0), (0, -1, 2)]`)
    intensities : Iterable of floats
        Experimental diffracted intensity for corresponding reflections.
    crystal : crystals.Crystal
        Crystal that gave rise to the diffracted intensities.
    mesh : 3-tuple ndarrays, ndim 2 or ndim 3
        Real-space mesh over which to calculate the scattering map.
        Format should be similar to the output of numpy.meshgrid.

    Returns
    -------
    out : ndarray, ndim 2 or ndim 3
        Electrostatic potential computed over the mesh.
    """
    assert len(intensities) == len(reflections)

    intensities = np.array(intensities)
    if np.any(intensities < 0):
        raise ValueError("Diffracted intensity cannot physically be negative.")

    # We want to support 2D and 3D meshes, therefore
    # expand mesh until 4D (last dim is for loop over reflections)
    # Extra dimensions will be squeezed out later
    xx, yy, zz = mesh
    while xx.ndim < 4:
        xx, yy, zz = (
            np.expand_dims(xx, xx.ndim),
            np.expand_dims(yy, yy.ndim),
            np.expand_dims(zz, zz.ndim),
        )

    # Reconstruct the structure factors from experimental data
    # We need to compute the theoretical phases from the crystal structure
    # To do this, we need to change 'reflections' into three iterables:
    # h, k, and l arrays
    hs, ks, ls = np.hsplit(np.array(reflections), 3)
    theoretical_SF = structure_factor(crystal, hs, ks, ls)
    phases = np.angle(theoretical_SF)
    experimental_SF = np.sqrt(intensities) * np.exp(1j * phases)
    experimental_SF = experimental_SF.reshape((1, 1, 1, -1))

    qx, qy, qz = change_basis_mesh(
        hs, ks, ls, basis1=crystal.reciprocal_vectors, basis2=np.eye(3)
    )
    qx, qy, qz = (
        qx.reshape((1, 1, 1, -1)),
        qy.reshape((1, 1, 1, -1)),
        qz.reshape((1, 1, 1, -1)),
    )
    p = np.sum(experimental_SF * np.cos(xx * qx + yy * qy + zz * qz), axis=3)
    return np.squeeze(np.real(p))
