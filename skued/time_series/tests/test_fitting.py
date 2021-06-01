from random import random, seed

import numpy as np

from skued import biexponential, exponential, with_irf

seed(23)


def test_exponential_tzero_limits():
    """Test that the output of ``exponential`` has the correct time-zero"""
    tzero = 10 * (random() - 0.5)  # between -5 and 5
    amp = 5 * random() + 5  # between 5 and 10
    tconst = random() + 0.3  # between 0.3 and 1.3

    t = np.arange(-10, 50, step=0.3)
    I = exponential(t, tzero=tzero, amp=amp, tconst=tconst)

    # Check that all values before time-zero are the amplitude
    assert np.all(np.equal(I[t < tzero], amp))
    assert np.all(np.less(I[t > tzero], amp))


def test_exponential_positivity():
    """Test that the output of ``exponential`` is always positive."""
    tzero = 10 * (random() - 0.5)  # between -5 and 5
    amp = 5 * random() + 5  # between 5 and 10
    tconst = random() + 0.3  # between 0.3 and 1.3

    t = np.arange(-10, 50, step=0.3)
    I = exponential(t, tzero=tzero, amp=amp, tconst=tconst)

    assert np.all(I > 0)


def test_exponential_amplitude():
    """Test that the output of ``exponential`` is at most ``amp``."""
    tzero = 10 * (random() - 0.5)  # between -5 and 5
    amp = 5 * random() + 5  # between 5 and 10
    tconst = random() + 0.3  # between 0.3 and 1.3

    t = np.arange(-10, 50, step=0.3)
    I = exponential(t, tzero=tzero, amp=amp, tconst=tconst)

    assert np.all(np.less_equal(I, amp))


def test_exponential_offset():
    """Test that the output of ``exponential`` is at lest ``offset``."""
    tzero = 10 * (random() - 0.5)  # between -5 and 5
    amp = 5 * random() + 5  # between 5 and 10
    tconst = random() + 0.3  # between 0.3 and 1.3

    offset = 15
    t = np.arange(-10, 50, step=0.3)
    I = exponential(t, tzero=tzero, amp=amp, tconst=tconst, offset=offset)

    assert np.all(np.greater_equal(I, offset))


def test_biexponential_tzero_limits():
    """Test that the output of ``biexponential`` has the correct time-zero"""
    tzero = 10 * (random() - 0.5)  # between -5 and 5
    amp1 = 5 * random() + 5  # between 5 and 10
    tconst1 = random() + 0.3  # between 0.3 and 1.3
    amp2 = 5 * random() + 5  # between 5 and 10
    tconst2 = random() + 0.3  # between 0.3 and 1.3

    t = np.arange(-10, 50, step=0.3)
    I = biexponential(
        t,
        tzero=tzero,
        amp1=amp1,
        amp2=amp2,
        tconst1=tconst1,
        tconst2=tconst2,
    )

    # Check that all values before time-zero are the amplitude
    assert np.all(np.equal(I[t < tzero], amp1 + amp2))
    assert np.all(np.less(I[t > tzero], amp1 + amp2))


def test_biexponential_positivity():
    """Test that the output of ``biexponential`` is always positive."""
    tzero = 10 * (random() - 0.5)  # between -5 and 5
    amp1 = 5 * random() + 5  # between 5 and 10
    tconst1 = random() + 0.3  # between 0.3 and 1.3
    amp2 = 5 * random() + 5  # between 5 and 10
    tconst2 = random() + 0.3  # between 0.3 and 1.3

    t = np.arange(-10, 50, step=0.3)
    I = biexponential(
        t,
        tzero=tzero,
        amp1=amp1,
        amp2=amp2,
        tconst1=tconst1,
        tconst2=tconst2,
    )

    assert np.all(I > 0)


def test_biexponential_amplitude():
    """Test that the output of ``biexponential`` is at most ``amp1 + amp2``."""
    tzero = 10 * (random() - 0.5)  # between -5 and 5
    amp1 = 5 * random() + 5  # between 5 and 10
    tconst1 = random() + 0.3  # between 0.3 and 1.3
    amp2 = 5 * random() + 5  # between 5 and 10
    tconst2 = random() + 0.3  # between 0.3 and 1.3

    t = np.arange(-10, 50, step=0.3)
    I = biexponential(
        t,
        tzero=tzero,
        amp1=amp1,
        amp2=amp2,
        tconst1=tconst1,
        tconst2=tconst2,
    )

    assert np.all(np.less_equal(I, amp1 + amp2))


def test_biexponential_offset():
    """Test that the output of ``biexponential`` is at least ``offset``."""
    tzero = 10 * (random() - 0.5)  # between -5 and 5
    amp1 = 5 * random() + 5  # between 5 and 10
    tconst1 = random() + 0.3  # between 0.3 and 1.3
    amp2 = 5 * random() + 5  # between 5 and 10
    tconst2 = random() + 0.3  # between 0.3 and 1.3

    offset = 15
    t = np.arange(-10, 50, step=0.3)
    I = biexponential(
        t,
        tzero=tzero,
        amp1=amp1,
        amp2=amp2,
        tconst1=tconst1,
        tconst2=tconst2,
        offset=offset,
    )

    assert np.all(np.greater_equal(I, offset))


def test_biexponential_against_exponential():
    """Test that ``biexponential`` reduces to ``exponential`` for appropriate parameters"""
    tzero = 10 * (random() - 0.5)  # between -5 and 5
    amp1 = 5 * random() + 5  # between 5 and 10
    tconst1 = random() + 0.3  # between 0.3 and 1.3
    amp2 = 5 * random() + 5  # between 5 and 10
    tconst2 = random() + 0.3  # between 0.3 and 1.3

    t = np.arange(-10, 50, step=0.3)
    offset = 2
    exp = exponential(t, tzero, amp1, tconst1, offset=offset)
    biexp = biexponential(t, tzero, amp1, 0, tconst1, 1, offset=offset)

    assert np.allclose(exp, biexp)


def test_with_irf_trivial_constant_spacing():
    """Test with_irf with a trivial IRF, with constant spacing"""
    params = (0, 1, 3)
    times = np.linspace(-5, 15, 256)
    data = exponential(times, *params)

    @with_irf(0.00001)  # vanishingly small irf
    def exponential_with_irf(time, *args, **kwargs):
        return exponential(time, *args, **kwargs)

    conv = exponential_with_irf(times, *params)

    assert np.allclose(data, conv)


def test_with_irf_trivial_nonconstant_spacing():
    """Test with_irf with a trivial IRF, with non-constant spacing"""
    # Note that the spacing of the steps is important for this test
    # If the array `times` is of even length, then the convolution will result
    # in one time-step shift
    params = (0, 1, 3)
    times = np.concatenate(
        (
            np.arange(-10, -2, step=1),
            np.arange(-2, 2, step=0.04),
            np.arange(2, 10, step=1),
        )
    )
    data = exponential(times, *params)

    @with_irf(0.00001)  # vanishingly small irf
    def exponential_with_irf(time, *args, **kwargs):
        return exponential(time, *args, **kwargs)

    conv = exponential_with_irf(times, *params)

    assert np.allclose(data, conv)
