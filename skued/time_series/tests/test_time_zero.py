
from .. import register_time_shift, register_time_shifts
import unittest
import numpy as np

from scipy.ndimage.interpolation import shift as scipy_shift

np.random.seed(23)

class TestTimeShift(unittest.TestCase):

    def test_trivial(self):
        """ Test that the time-shift between two identical traces is zero. """
        with self.subTest('Even length'):
            trace1 = np.sin(2*np.pi*np.linspace(0, 10, 64))
            trace2 = np.array(trace1, copy = True)
            shift = register_time_shift(trace1, trace2)
            self.assertEqual(shift, 0)
        
        with self.subTest('Odd length'):
            trace1 = np.sin(2*np.pi*np.linspace(0, 10, 65))
            trace2 = np.array(trace1, copy = True)
            shift = register_time_shift(trace1, trace2)
            self.assertEqual(shift, 0)
    
    def test_shift_no_noise(self):
        """ Test measuring the time-shift between traces shifted from one another, without added noise """
        trace1 = np.sin(2*np.pi*np.linspace(0, 10, 64))
        trace2 = np.roll(trace1, 5)
        shift = register_time_shift(trace1, trace2)
        self.assertEqual(shift, -5)

    def test_shift_with_noise(self):
        """ Test measuring the time-shift between traces shifted from one another, with added 10% gaussian noise """
        trace1 = np.sin(2*np.pi*np.linspace(0, 10, 64))
        trace2 = scipy_shift(trace1, 5)

        trace1 = trace1[6:-6]
        trace2 = trace2[6:-6]

        trace1 += 0.05*np.random.random(size = trace1.shape)
        trace2 += 0.05*np.random.random(size = trace2.shape)
        shift = register_time_shift(trace1, trace2)
        self.assertEqual(shift, -5)
    
    def test_shift_different_lengths(self):
        """ Test that register_time_shift() raises an exception if the reference and trace do not have the same shape """
        with self.assertRaises(ValueError):
            trace1 = np.empty((16,))
            trace2 = np.empty((8,))
            register_time_shift(trace1, trace2)

    def test_shift_not1d(self):
        """ Test that register_time_shift() raises an exception if the reference or trace are not 1D """
        with self.subTest('Trace not 1D'):
            with self.assertRaises(ValueError):
                trace1 = np.empty((16,45))
                trace2 = np.empty((8,))
                register_time_shift(trace1, trace2)

        with self.subTest('Reference not 1D'):
            with self.assertRaises(ValueError):
                trace1 = np.empty((16,))
                trace2 = np.empty((8,2))
                register_time_shift(trace1, trace2)

class TestTimeShifts(unittest.TestCase):
    
    def test_trivial(self):
        """ Test that the time-shifts between identical time traces """
        with self.subTest('Even lengths'):
            traces = [np.sin(2*np.pi*np.linspace(0, 10, 64)) for _ in range(10)]
            shifts = register_time_shifts(traces)
            self.assertTrue(np.allclose(shifts, np.zeros_like(shifts)))

        with self.subTest('Odd lengths'):
            traces = [np.sin(2*np.pi*np.linspace(0, 10, 31)) for _ in range(10)]
            shifts = register_time_shifts(traces)
            self.assertTrue(np.allclose(shifts, np.zeros_like(shifts)))
    
    def test_output_shape(self):
        """ Test the output shape """
        with self.subTest('reference = None'):
            traces = [np.sin(2*np.pi*np.linspace(0, 10, 64) + i) for i in range(10)]
            shifts = register_time_shifts(traces)
            self.assertTupleEqual(shifts.shape, (len(traces), ))
            # The first shift should then be zero
            # because it is the shift between the reference and itself
            self.assertEqual(shifts[0], 0)

        with self.subTest('reference is not None'):
            traces = [np.sin(2*np.pi*np.linspace(0, 10, 64) + i) for i in range(10)]
            shifts = register_time_shifts(traces, reference = np.array(traces[0], copy = True))
            self.assertTupleEqual(shifts.shape, (len(traces), ))



if __name__ == '__main__':
    unittest.main()