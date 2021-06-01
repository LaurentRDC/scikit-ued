import numpy as np

from skued import nfft, nfftfreq

np.random.seed(23)


def test_nfft_shape_even():
    """Test that the nfftfreq function returns expected shape"""
    freqs = nfftfreq(16)
    assert freqs.shape == (16,)


def test_nfft_shape_odd():
    """Test that the nfftfreq function returns expected shape"""
    freqs = nfftfreq(13)
    assert freqs.shape == (13,)


def test_nfftfreq_against_fft():
    """Test against goold ol' FFT on evenly spaced data"""

    x = np.linspace(0, 10, num=128)
    y = np.sin(2 * np.pi * x)

    k = np.fft.fftfreq(y.size)
    knfft = nfftfreq(len(x), df=1)

    dft = np.fft.fft(y)
    from_nfft = nfft(x, y, M=len(x))
