# -*- coding: utf-8 -*-

from collections.abc import Sized
from functools import partial

import numpy as np
from scipy.signal import correlate

from npstreams import peek, array_stream


def time_shift(trace, reference, method = 'auto'):
    """ 
    Measure the time shift between a time trace and a reference trace
    by cross correlation.
    
    Parameters
    ----------
    trace : array-like
        Time trace. Must be the same lengths as ``reference``.
    reference : array-like
        Reference trace.
    method : str {'auto', 'fft', 'direct'}, optional
        A string indicating which method to use to calculate the correlation.
    
    Returns
    -------
    shift : int or float
        Index shift between the two traces, with respect to
        the ``reference`` trace. Note that the shift can be fractional.
    
    Raises
    ------
    ValueError : if ``trace`` and ``reference`` do not have the same shape.

    See Also
    --------
    time_shifts : measure time-shift between multiple traces and a reference.
    scipy.signal.choose_conv_method : contains more documentation on ``method``
    """
    trace, reference = np.atleast_1d(trace, reference)
    if trace.shape != reference.shape:
        raise ValueError('Time trace and reference trace are expected to have the same shape, be received \
                         a time-trace of shape {} and a reference trace of shape {}'.format(trace.shape, reference.shape))

    xcorr = correlate(trace, reference, mode = 'same', method = method)

    # Generalize to the average of multiple maxima
    maxima = np.transpose(np.nonzero(xcorr == xcorr.max())) 
    return np.mean(maxima) - int(xcorr.shape[0]/2)

@array_stream
def time_shifts(traces, reference = None, method = 'auto'):
    """
    Measure the time shifts between time traces and a reference by cross-correlation.

    Parameters
    ----------
    traces : iterable of ndarrays
        Time traces. These time-traces should be physically equivalent. Generators of time traces
        are also supported. All traces and ``reference`` must have the same shape.
    reference : `~numpy.ndarray` or None, optional
        If provided, the time-zero shift between the traces in ``traces`` will be measured
        with respect to ``reference``. Otherwise, the first trace in ``traces`` will be used 
        as a reference.
    method : str {'auto', 'fft', 'direct'}, optional
        A string indicating which method to use to calculate the correlation.
    
    Returns
    -------
    shifts : `~numpy.ndarray`, ndim 1, dtype float
        Time shifts as indices (possibly fractional). The length of ``shifts`` is always 
        equal to the number of time-traces; in the case where ``reference = None``, the first
        shifts will always be identically zero.
    
    Raises
    ------
    ValueError : if not all traces have the same shape.

    See Also
    --------
    time_shift : measure time-shift between a single trace and a reference.
    """
    # fromiter can preallocate the full array if the number of traces
    # is known in advance
    try:
        count = len(traces)
    except TypeError:
        count = -1

    traces = iter(traces)

    if reference is None:
        reference, traces = peek(traces)
    reference = np.atleast_1d(reference)

    kwargs = {'reference': reference, 'method': method}

    shifts = map(partial(time_shift, **kwargs), traces)
    return np.fromiter(shifts, dtype = np.float, count = count)