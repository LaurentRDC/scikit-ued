# -*- coding: utf-8 -*-

from functools import partial

import numpy as np
from scipy.signal import correlate

from npstreams import array_stream, peek


def time_shift(trace, reference):
    """ 
    Measure the time shift between a time trace and a reference trace
    by cross correlation.
    
    Parameters
    ----------
    trace : array-like
        Time trace. Must be the same lengths as ``reference``.
    reference : array-like
        Reference trace.
    
    Returns
    -------
    shift : int
        Index shift between the two traces, with respect to
        the ``reference`` trace. 
    
    Raises
    ------
    ValueError : if ``trace`` and ``reference`` do not have the same shape.

    See Also
    --------
    time_shifts : measure time-shift between multiple traces and a reference.
    """
    trace, reference = np.atleast_1d(trace, reference)
    if trace.shape != reference.shape:
        raise ValueError('Time trace and reference trace are expected to have the same shape, be received \
                         a time-trace of shape {} and a reference trace of shape {}'.format(trace.shape, reference.shape))

    xcorr = correlate(trace, reference, mode = 'same')

    # Generalize to the average of multiple maxima
    maxima = np.transpose(np.nonzero(xcorr == xcorr.max())) 
    return round(np.mean(maxima) - int(xcorr.shape[0]/2))

@array_stream
def time_shifts(traces, reference = None):
    """
    Measure the time shifts between time traces and a reference.

    Parameters
    ----------
    traces : iterable of ndarrays
        Time traces. These time-traces should be physically equivalent. Generators of time traces
        are also supported. All traces and ``reference`` must have the same shape.
    reference : `~numpy.ndarray` or None, optional
        If provided, the time-zero shift between the traces in ``traces`` will be measured
        with respect to ``reference``. Otherwise, the first trace in ``traces`` will be used 
        as a reference.
    
    Returns
    -------
    shifts : `~numpy.ndarray`, ndim 1, dtype int
        Time shifts as indices.
    
    Raises
    ------
    ValueError : if not all traces have the same shape.

    See Also
    --------
    time_shift : measure time-shift between a single trace and a reference.
    """
    traces = iter(traces)

    if reference is None:
        reference = next(traces)
    reference = np.atleast_1d(reference)

    shifts = map(partial(time_shift, reference = reference), traces)
    return np.fromiter(shifts, dtype = np.int)
