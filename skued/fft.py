# -*- coding: utf-8 -*-
"""
Common module for fast fourier transforms.
"""

import scipy.fft
from functools import wraps
from scipy.fft import _pocketfft
from os import cpu_count

CPU_COUNT = cpu_count()


class SkuedPocketFFTBackend:
    """
    FFT backend entirely based on Scipy's PocketFFT, with better defaults.

    The speed of two-dimensional transforms may be improved by up to 50% with these
    different defaults.

    See also
    --------
    with_skued_fft
    """

    __ua_domain__ = "numpy.scipy.fft"

    @staticmethod
    def __ua_function__(method, args, kwargs):
        fn = getattr(_pocketfft, method.__name__, None)

        if fn is None:
            return NotImplemented
        workers = kwargs.pop("workers", CPU_COUNT)
        return fn(*args, workers=workers, **kwargs)


def with_skued_fft(f):
    """Ensure the use of the SkuedPocketFFTBackend whenever the `scipy.fft` module is used."""

    @wraps(f)
    def newf(*args, **kwargs):
        with scipy.fft.set_backend(SkuedPocketFFTBackend):
            return f(*args, **kwargs)

    return newf
