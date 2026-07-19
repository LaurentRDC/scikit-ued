# -*- coding: utf-8 -*-
"""
Common module for fast fourier transforms.
"""

import scipy.fft
from functools import wraps


def with_scipy_fft(f):

    @wraps(f)
    def newf(*args, **kwargs):
        with scipy.fft.set_backend("scipy"):
            return f(*args, **kwargs)

    return newf
