import unittest
from random import random, seed

import numpy as np

from skued import biexponential, exponential, with_irf

seed(23)


class TestExponentialDecay(unittest.TestCase):
    def setUp(self):
        self.tzero = 10 * (random() - 0.5)  # between -5 and 5
        self.amp = 5 * random() + 5  # between 5 and 10
        self.tconst = random() + 0.3  # between 0.3 and 1.3

    def test_tzero_limits(self):
        """ Test that the output of ``exponential`` has the correct time-zero """
        t = np.arange(-10, 50, step=0.3)
        I = exponential(t, tzero=self.tzero, amp=self.amp, tconst=self.tconst)

        # Check that all values before time-zero are the amplitude
        self.assertTrue(np.all(np.equal(I[t < self.tzero], self.amp)))
        self.assertTrue(np.all(np.less(I[t > self.tzero], self.amp)))

    def test_positivity(self):
        """ Test that the output of ``exponential`` is always positive. """
        t = np.arange(-10, 50, step=0.3)
        I = exponential(t, tzero=self.tzero, amp=self.amp, tconst=self.tconst)

        self.assertTrue(np.all(I > 0))

    def test_amplitude(self):
        """ Test that the output of ``exponential`` is at most ``amp``. """
        t = np.arange(-10, 50, step=0.3)
        I = exponential(t, tzero=self.tzero, amp=self.amp, tconst=self.tconst)

        self.assertTrue(np.all(np.less_equal(I, self.amp)))

    def test_offset(self):
        """ Test that the output of ``exponential`` is at lest ``offset``. """
        offset = 15
        t = np.arange(-10, 50, step=0.3)
        I = exponential(
            t, tzero=self.tzero, amp=self.amp, tconst=self.tconst, offset=offset
        )

        self.assertTrue(np.all(np.greater_equal(I, offset)))


class TestBiExponentialDecay(unittest.TestCase):
    def setUp(self):
        self.tzero = 10 * (random() - 0.5)  # between -5 and 5
        self.amp1 = 5 * random() + 5  # between 5 and 10
        self.tconst1 = random() + 0.3  # between 0.3 and 1.3
        self.amp2 = 5 * random() + 5  # between 5 and 10
        self.tconst2 = random() + 0.3  # between 0.3 and 1.3

    def test_tzero_limits(self):
        """ Test that the output of ``biexponential`` has the correct time-zero """
        t = np.arange(-10, 50, step=0.3)
        I = biexponential(
            t,
            tzero=self.tzero,
            amp1=self.amp1,
            amp2=self.amp2,
            tconst1=self.tconst1,
            tconst2=self.tconst2,
        )

        # Check that all values before time-zero are the amplitude
        self.assertTrue(np.all(np.equal(I[t < self.tzero], self.amp1 + self.amp2)))
        self.assertTrue(np.all(np.less(I[t > self.tzero], self.amp1 + self.amp2)))

    def test_positivity(self):
        """ Test that the output of ``biexponential`` is always positive. """
        t = np.arange(-10, 50, step=0.3)
        I = biexponential(
            t,
            tzero=self.tzero,
            amp1=self.amp1,
            amp2=self.amp2,
            tconst1=self.tconst1,
            tconst2=self.tconst2,
        )

        self.assertTrue(np.all(I > 0))

    def test_amplitude(self):
        """ Test that the output of ``biexponential`` is at most ``amp1 + amp2``. """
        t = np.arange(-10, 50, step=0.3)
        I = biexponential(
            t,
            tzero=self.tzero,
            amp1=self.amp1,
            amp2=self.amp2,
            tconst1=self.tconst1,
            tconst2=self.tconst2,
        )

        self.assertTrue(np.all(np.less_equal(I, self.amp1 + self.amp2)))

    def test_offset(self):
        """ Test that the output of ``biexponential`` is at least ``offset``. """
        offset = 15
        t = np.arange(-10, 50, step=0.3)
        I = biexponential(
            t,
            tzero=self.tzero,
            amp1=self.amp1,
            amp2=self.amp2,
            tconst1=self.tconst1,
            tconst2=self.tconst2,
            offset=offset,
        )

        self.assertTrue(np.all(np.greater_equal(I, offset)))

    def test_against_exponential(self):
        """ Test that ``biexponential`` reduces to ``exponential`` for appropriate parameters """
        t = np.arange(-10, 50, step=0.3)
        offset = 2
        exp = exponential(t, self.tzero, self.amp1, self.tconst1, offset=offset)
        biexp = biexponential(
            t, self.tzero, self.amp1, 0, self.tconst1, 1, offset=offset
        )

        self.assertTrue(np.allclose(exp, biexp))


class TestWithIrf(unittest.TestCase):
    def test_trivial_constant_spacing(self):
        """ Test with_irf with a trivial IRF, with constant spacing """
        params = (0, 1, 3)
        times = np.linspace(-5, 15, 256)
        data = exponential(times, *params)

        @with_irf(0.00001)  # vanishingly small irf
        def exponential_with_irf(time, *args, **kwargs):
            return exponential(time, *args, **kwargs)

        conv = exponential_with_irf(times, *params)

        self.assertTrue(np.allclose(data, conv))

    def test_trivial_nonconstant_spacing(self):
        """ Test with_irf with a trivial IRF, with non-constant spacing """
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

        self.assertTrue(np.allclose(data, conv))


if __name__ == "__main__":
    unittest.main()
