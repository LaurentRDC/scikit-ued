from .. import nfft, nfftfreq
import unittest
import numpy as np

np.random.seed(23)

class Testnfftfreq(unittest.TestCase):
    
    def test_shape_even(self):
        """ Test that the nfftfreq function returns expected shape """
        freqs = nfftfreq(16)
        self.assertTupleEqual(freqs.shape, (16,))

    def test_shape_odd(self):
        """ Test that the nfftfreq function returns expected shape """
        freqs = nfftfreq(13)
        self.assertTupleEqual(freqs.shape, (13,))

class Testnfft(unittest.TestCase):

    def test_against_fft(self):
        """ Test against goold ol' FFT on evenly spaced data """

        x = np.linspace(0, 10, num = 128)
        y = np.sin(2*np.pi*x)

        k = np.fft.fftfreq(y.size)
        knfft = nfftfreq(len(x), df = 1)

        dft = np.fft.fft(y)
        from_nfft = nfft(x,y, M = len(x))

if __name__ == '__main__':
    unittest.main()