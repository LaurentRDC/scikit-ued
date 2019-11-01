# -*- coding: utf-8 -*-
"""
Calculation of one-phonon structure factors
"""

from .phonons import PhononMode, rowdot
import numpy as np
from ..simulation import affe


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

    # Assuming that the crystals was parsed from a PWSCF results file,
    # We loop through atoms in order that they are visible
    try:
        atoms = sorted(crystal, key=lambda a: a.tag)
    except TypeError: # if atoms don't have a tag value
        atoms = crystal
        
    accumulator = np.zeros(shape=(qpoints.shape[0], 1), dtype=np.complex)
    for atm_index, atm in enumerate(atoms):
        # Accumulator is built in pieces
        # because all of these arrays can be pretty big
        arg = np.ones_like(
            accumulator, dtype=np.complex
        )  # because polarization are complex vectors

        arg *= np.exp(-1 * dw_factors[atm_index].reshape(-1, 1))
        arg *= affe(atm, q_norm) / np.sqrt(atm.mass)
        arg *= rowdot(qpoints, polarizations[:, atm_index, :])
        accumulator += arg

    return np.nan_to_num(accumulator)
