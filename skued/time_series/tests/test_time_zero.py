import numpy as np
from scipy.ndimage import shift as scipy_shift

from skued import register_time_shift, register_time_shifts
import pytest

np.random.seed(23)


def test_time_shift_trivial():
    """Test that the time-shift between two identical traces is zero."""
    trace1 = np.sin(2 * np.pi * np.linspace(0, 10, 64))
    trace2 = np.array(trace1, copy=True)
    shift = register_time_shift(trace1, trace2)
    assert shift == 0

    trace1 = np.sin(2 * np.pi * np.linspace(0, 10, 65))
    trace2 = np.array(trace1, copy=True)
    shift = register_time_shift(trace1, trace2)
    assert shift == 0


def test_time_shift_no_noise():
    """Test measuring the time-shift between traces shifted from one another, without added noise"""
    trace1 = np.sin(2 * np.pi * np.linspace(0, 10, 64))
    trace2 = np.roll(trace1, 5)
    shift = register_time_shift(trace1, trace2)
    assert shift == -5


def test_time_shift_with_noise():
    """Test measuring the time-shift between traces shifted from one another, with added 6% gaussian noise"""
    trace1 = np.sin(2 * np.pi * np.linspace(0, 10, 64))
    trace2 = scipy_shift(trace1, 5)

    trace1 = trace1[6:-6]
    trace2 = trace2[6:-6]

    trace1 += 0.03 * np.random.random(size=trace1.shape)
    trace2 += 0.03 * np.random.random(size=trace2.shape)
    shift = register_time_shift(trace1, trace2)
    assert shift == -5


def test_time_shift_different_lengths():
    """Test that register_time_shift() raises an exception if the reference and trace do not have the same shape"""
    with pytest.raises(ValueError):
        trace1 = np.empty((16,))
        trace2 = np.empty((8,))
        register_time_shift(trace1, trace2)


def test_time_shift_not1d():
    """Test that register_time_shift() raises an exception if the reference or trace are not 1D"""
    with pytest.raises(ValueError):
        trace1 = np.empty((16, 45))
        trace2 = np.empty((8,))
        register_time_shift(trace1, trace2)

    with pytest.raises(ValueError):
        trace1 = np.empty((16,))
        trace2 = np.empty((8, 2))
        register_time_shift(trace1, trace2)


def test_time_shifts_trivial():
    """Test that the time-shifts between identical time traces"""
    traces = [np.sin(2 * np.pi * np.linspace(0, 10, 64)) for _ in range(10)]
    shifts = register_time_shifts(traces)
    assert np.allclose(shifts, np.zeros_like(shifts))

    traces = [np.sin(2 * np.pi * np.linspace(0, 10, 31)) for _ in range(10)]
    shifts = register_time_shifts(traces)
    assert np.allclose(shifts, np.zeros_like(shifts))


def test_time_shifts_output_shape():
    """Test the output shape"""
    traces = [np.sin(2 * np.pi * np.linspace(0, 10, 64) + i) for i in range(10)]
    shifts = register_time_shifts(traces)
    assert shifts.shape == (len(traces),)
    # The first shift should then be zero
    # because it is the shift between the reference and it
    assert shifts[0] == 0

    traces = [np.sin(2 * np.pi * np.linspace(0, 10, 64) + i) for i in range(10)]
    shifts = register_time_shifts(traces, reference=np.array(traces[0], copy=True))
    assert shifts.shape == (len(traces),)
