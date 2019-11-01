# -*- coding: utf-8 -*-
"""
Calculation of one-phonon structure factors
"""

from .phonons import PhononMode
import numpy as np


def phonon_amplitude(frequencies, temperature):
    """
    Phonon amplitude within the Debye-Waller factor

    Parameters
    ----------
    frequencies : ndarray, shape (N,)
        Frequencies of a single phonon mode [Hz].
    temperature : float
        Mode temperature [K].
    
    Returns
    -------
    amplitude : ndarray, shape (N,)
        Amplitude SQUARED

    References
    ----------
    Xu and Chiang (2005) Eq. 23
    """
    # Factor of 1/N is calculated in the parent caller
    return (HBAR / frequencies) * coth(
        HBAR * frequencies / (2 * BOLTZMANN * temperature)
    )


def _debye_waller_factor(modes, temperatures, atm_index):
    """ Calculate a Debye-Waller factor for one atom. """
    # This calculation assumes that the Debye-Waller factor is isotropic
    # i.e. it only depends on the magnitude of q and polarizations ek
    # The anisotropic factor is much more computationally expensive
    n = ncells(modes["LA"].crystal)

    # The sum happens for all q's, but the polarization vectors of a single
    # Brillouin zone. This is annoying to keep track of. Since polarization
    # vectors are the same across different Brillouin zones, we sum over all
    # zones. This leads to "double" counting, hence a correction factor based
    # on the number of Brilluoin zones `nzones`.
    hkls, *_ = unique_by_rows(modes["LA"].hkls)
    nzones = hkls.shape[0]

    def accumulator(mode):
        # This is the sum over k
        # Really, it is a sum over all q, and it will be corrected
        # in the parent function.
        temp = temperatures[mode.name]
        return (
            np.sum(
                (1 / n)
                * phonon_amplitude(mode.frequencies, temp)
                * np.linalg.norm(mode.polarizations[:, atm_index, :], axis=1) ** 2
            )
            / nzones
        )

    # This is the sum over modes
    return ns.sum(accumulator(m) for m in modes.values())


def debye_waller_factors(modes, temperatures=None):
    """
    Compute the debye-waller factor based on all mode information.
    These modes are assumed to have been expanded, i.e. represent mode information
    over the entire detector range.

    For performance reasons, we consider the isotropic Debye-Waller effect.

    Parameters
    ----------
    modes : dict[str, Mode]
    temperatures : dict[str, float] or None, optional
        Mode temperatures [K]. Defaults to room temperature.

    Returns
    -------
    factors : iterable of narrays
        One factor for every atom in the unit cell.

    References
    ----------
    Xu and Chiang (2005) Eq. 19 (anisotropic) and Eq. 20 (isotropic)
    """
    if temperatures is None:
        temperatures = {m: 300 for m in modes.keys()}

    # We loop through atoms in order that they are visible in the PWSCF file
    # That's the tag properties on Atom objects
    atoms = sorted(modes["LA"].crystal, key=lambda a: a.tag)
    q2 = np.linalg.norm(modes["LA"].q_points, axis=1) ** 2
    prefactor = lambda atm: q2 / (12 * atm.mass * AMU_TO_KG)

    # TODO: should we filter modes out-of-planes?

    # Parallelizing this calculation is actually slower
    # The correction factor `nzones` accounts for the "double" counting
    # of multiple BZ in the sum
    return [
        prefactor(atom) * _debye_waller_factor(modes, temperatures, index)
        for index, atom in enumerate(atoms)
    ]


def one_phonon_structure_factor(mode, dw_factors):
    """ 
    Compute the one-phonon structure factor associated with a phonon mode. 
    
    Parameters
    ----------
    mode : PhononMode
        Mode defined at `N` q-points.
    dw_factors : iterable of ndarray, shapes (N,)
        Debye-Waller factors, at every q-point of `mode`, for each atom in the unit cell. 
        These factors must be computed separately because they require knowledge of many modes.
    
    Returns
    -------
    oneph: ndarray, shape (N,)
        One-phonon structure factor for `mode`.
    """
    qpoints, polarizations, hkls, crystal = (
        mode.q_points,
        mode.polarizations,
        mode.hkls,
        mode.crystal,
    )

    assert dw_factors[0].shape == (qpoints.shape[0],)
    assert qpoints.shape == hkls.shape
    assert polarizations[:, 0, :].shape == qpoints.shape

    q_norm = np.linalg.norm(qpoints, axis=1, keepdims=True)

    # We loop through atoms in order that they are visible in the PWSCF file
    atoms = sorted(crystal, key=lambda a: a.tag)
    accumulator = np.zeros(shape=(qpoints.shape[0], 1), dtype=np.complex)
    for atm_index, atm in enumerate(atoms):
        # Accumulator is built in pieces
        # because all of these arrays are pretty big
        arg = np.ones_like(
            accumulator, dtype=np.complex
        )  # because polarization are complex vectors

        arg *= np.exp(-1 * dw_factors[atm_index].reshape(-1, 1))
        arg *= affe(atm, q_norm) / np.sqrt(atm.mass)
        arg *= rowdot(qpoints, polarizations[:, atm_index, :])
        accumulator += arg

    return np.nan_to_num(accumulator)
