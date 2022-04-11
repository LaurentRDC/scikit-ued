# -*- coding: utf-8 -*-
"""
Voigt and pseudo-voigt curves, as well as related Gaussian and Lorentzian functions
"""
from functools import lru_cache

import numpy as np
from numpy import pi


def gaussian(coordinates, center, fwhm=None, std=None):
    """
    Unit integral Gaussian function.

    Parameters
    ----------
    coordinates : ndarray or list of ndarrays
            Can be either a list of ndarrays, as a meshgrid coordinates list, or a
            single ndarray for 1D computation
    center : array_like
            Center of the gaussian. Should be the same shape as `coordinates.ndim`.
    fwhm : float or None, optional
            Full-width at half-max of the function. Either `std` or `fwhm` must be provided.
    std : float or None, optional
            Standard deviation of the function. Either `std` or `fwhm` must be provided.

    Returns
    -------
    out : ndarray
            Gaussian function of unit integral.

    Raises
    ------
    ValueError : If fwhm and std are not provided.

    Notes
    -----
    In the case where both `std` and `fwhm` are given, `fwhm` takes precedence.

    Examples
    --------
    >>> import numpy as np
    >>> from skued import gaussian
    >>>
    >>> span = np.arange(-10, 10, 0.1)
    >>> xx, yy = np.meshgrid(span, span)
    >>> g = gaussian( coordinates = [xx,yy], center = [0,0], std = 1)
    >>> g.shape == xx.shape
    True
    >>> np.sum(g)*0.1**2   # Integral should be close to unity (spacing = 0.1)
    1.0000000000000075
    """
    if all([fwhm is None, std is None]):
        raise ValueError("Either fwhm or std has to be provided")

    if fwhm is not None:
        std = fwhm / (2 * np.sqrt(2 * np.log(2)))

        # 1D is a special case, as coordinates are not given as a list of arrays
    if not isinstance(coordinates, (list, tuple)):  # iterable but not ndarray
        return (
            1
            / (std * np.sqrt(2 * pi))
            * np.exp(-((coordinates - center) ** 2) / (2 * std * std))
        )

        # Computation
    dim = len(coordinates)
    exponent = sum([(x - c) ** 2 for x, c in zip(coordinates, center)]) / (
        2 * std * std
    )
    factor = 1 / (std * np.sqrt(2 * pi)) ** dim
    return factor * np.exp(-exponent)


def lorentzian(coordinates, center, fwhm):
    """
    Unit integral Lorenzian function.

    Parameters
    ----------
    coordinates : array-like
            Can be either a list of ndarrays, as a meshgrid coordinates list, or a
            single ndarray for 1D computation
    center : array-like
            Center of the lorentzian. Should be the same shape as `coordinates.ndim`.
    fwhm : float
            Full-width at half-max of the function.

    Returns
    -------
    out : ndarray
            Lorentzian function of unit integral.

    Notes
    -----
    The functional form of the Lorentzian is given by:

    .. math::

            L(x) = \\frac{1}{\pi} \\frac{(\gamma/2)}{(x-c)^2 + (\gamma/2)^2}

    where :math:`\gamma` is the full-width at half-maximum, and :math:`c` is the
    center.

    For n dimensions, the functional form of the Lorentzian is given by:

    .. math::

            L(x_1, ..., x_n) = \\frac{1}{n \pi} \\frac{(\gamma/2)}{(\sum_i{(x_i - c_i)^2} + (\gamma/2)^2)^{\\frac{1+n}{2}}}

    Example
    -------
    >>> import numpy as np
    >>> from skued import lorentzian
    >>>
    >>> span = np.arange(-15, 15, 0.1)
    >>> xx, yy = np.meshgrid(span, span)
    >>> l = lorentzian( coordinates = [xx,yy], center = [0,0], fwhm = 1)
    >>> l.shape == xx.shape
    True
    >>> np.sum(l)*0.1**2          # Integral should be unity (spacing = 0.1)
    0.9700030627398781
    """
    width = 0.5 * fwhm

    # 1D is a special case, as coordinates are not given as a list of arrays
    if not isinstance(coordinates, (list, tuple)):  # iterable but not ndarray
        return (width / pi) / ((coordinates - center) ** 2 + width**2)

    dim = len(coordinates)
    core = width / (
        (sum([(x - c) ** 2 for x, c in zip(coordinates, center)]) + width**2)
    ) ** ((dim + 1) / 2)
    factor = 1 / (dim * pi)
    return factor * core


@lru_cache(maxsize=16)
def _pseudo_voigt_mixing_factor(width_l, width_g):
    """
    Returns the proportion of Lorentzian for the computation of a pseudo-Voigt profile.
    pseudoVoigt = (1 - eta) Gaussian + eta * Lorentzian

    Parameters
    ----------
    width_l, width_g : float
        FWHM for the Gaussian and Lorentzian parts, respectively.

    Returns
    -------
    eta : numerical
        Proportion of Lorentzian. Between 0 and 1
    """
    # Fast formula (see paper in pseudo_voigt docstrings)
    # This assumes width_g and width_l are the Gaussian FWHM and Lorentzian FWHM
    gamma = (
        width_g**5
        + 2.69269 * width_l * (width_g**4)
        + 2.42843 * (width_l**2) * (width_g**3)
        + 4.47163 * (width_l**3) * (width_g**2)
        + 0.07842 * (width_l**4) * width_g
        + width_l**5
    ) ** (1 / 5)

    # Proportion of the Voigt that should be Lorentzian
    return (
        1.36603 * (width_l / gamma)
        - 0.47719 * (width_l / gamma) ** 2
        + 0.11116 * (width_l / gamma) ** 3
    )


def pseudo_voigt(coordinates, center, fwhm_g, fwhm_l):
    """
    Unit integral pseudo-Voigt profile. Deviation from real Voigt
    by less than 1% [1]_.

    Parameters
    ----------
    coordinates : array_like
            Can be either a list of ndarrays, as a meshgrid coordinates list, or a
            single ndarray for 1D computation
    center : array_like
            Center of the pseudo-voigt. Should be the same shape as `coordinates.ndim`.
    fwhm_g, fwhm_l : float
                    Full-width at half-max of the Gaussian and Lorentzian parts respectively.

    Returns
    -------
    out : ndarray
            Pseudo-Voigt profile of unit integral.

    Example
    --------
    >>> import numpy as np
    >>> from skued import pseudo_voigt
    >>>
    >>> span = np.arange(-5, 5, 0.01)
    >>> xx, yy = np.meshgrid(span, span)
    >>> pV = pseudo_voigt( coordinates = [xx,yy], center = [0,0], fwhm_g = 0.1, fwhm_l = 0.1)
    >>> pV.shape == xx.shape
    True
    >>> np.trapz(np.trapz(pV, span[None, :], axis=1), span, axis=0) # Integral should be unity
    0.9938490810821982

    References
    ----------
    .. [1] T. Ida et al., Extended pseudo-Voigt function for approximating the Voigt profile.
            J. of Appl. Cryst. (2000) vol. 33, pp. 1311-1316
    """
    eta = _pseudo_voigt_mixing_factor(fwhm_g, fwhm_l)
    return (1 - eta) * gaussian(coordinates, center, fwhm_g) + eta * lorentzian(
        coordinates, center, fwhm_l
    )
