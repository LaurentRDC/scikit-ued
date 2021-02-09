# -*- coding: utf-8 -*-

from functools import lru_cache, partial

import numpy as np
from npstreams import array_stream, peek
from scipy.signal import correlate


# Save the normalization of correlations so that identical
# autocorrelations are saved.
@lru_cache(maxsize=128)
def __xcorr_normalization(size, dtype):
    arr = np.ones(shape=(size,), dtype=dtype)
    return correlate(arr, arr, mode="full")


def register_time_shift(trace, reference, method="auto"):
    """
    Measure the time shift between a time trace and a reference trace
    by normalized cross correlation.

    .. versionadded:: 1.0.1.1

    Parameters
    ----------
    trace : array-like, shape (N,)
        Time trace. Must be the same lengths as ``reference``. Only 1D traces are supported.
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
    ValueError : if ``trace`` is not a 1D array

    See Also
    --------
    register_time_shifts : measure time-shift between multiple traces and a reference.
    scipy.signal.choose_conv_method : contains more documentation on ``method``
    """
    trace, reference = np.atleast_1d(trace, reference)
    if trace.shape != reference.shape:
        raise ValueError(
            f"Time trace and reference trace are expected to have the same shape, be received \
                         a time-trace of shape {trace.shape} and a reference trace of shape {reference.shape}"
        )

    if trace.ndim > 1:
        raise ValueError(
            f"Expected 1D time traces, but received traces of shape {trace.shape}"
        )

    trace = trace - trace.mean()
    reference = reference - reference.mean()

    # Normalized cross-correlation
    # Note : we use an external function to calculate normalization
    #        so that it can be efficiently cached
    xcorr = correlate(
        trace, reference, mode="full", method="auto"
    ) / __xcorr_normalization(trace.size, trace.dtype)

    # Generalize to the average of multiple maxima
    maxima = np.transpose(np.nonzero(xcorr == xcorr.max()))
    return np.mean(maxima) - int(xcorr.shape[0] / 2)


@array_stream
def register_time_shifts(traces, reference=None, method="auto"):
    """
    Measure the time shifts between time traces and a reference by cross-correlation.

    .. versionadded:: 1.0.1.1

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
    ValueError : if traces are not 1D arrays

    See Also
    --------
    register_time_shift : measure time-shift between a single trace and a reference.
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

    kwargs = {"reference": reference, "method": method}

    shifts = map(partial(register_time_shift, **kwargs), traces)
    return np.fromiter(shifts, dtype=float, count=count)
