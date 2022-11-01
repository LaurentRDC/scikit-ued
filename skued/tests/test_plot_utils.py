# -*- coding: utf-8 -*-


from skued import rgb_sweep, spectrum_colors, spectrum_cmap
import matplotlib


def test_spectrum_colors_on_ints():
    """Test spectrum_colors on an int"""
    colors = spectrum_colors(10)
    assert len(list(colors)) == 10


def test_spectrum_colorson_sized_iterable():
    """Test on iterable that has a __len__ attribute: list, tuple, etc."""
    colors = spectrum_colors([1, 2, 3, 4, 5])
    assert len(list(colors)) == 5


def test_spectrum_colorson_unsized_iterable():
    """Test spectrum_colors on unsized_iterable (e.g. generator)"""
    colors = spectrum_colors(range(0, 10))
    assert len(list(colors)) == 10


def test_spectrum_colorson_length_1_iterable():
    """Test that spectrum_colors is working on single-length iterables"""
    assert list(spectrum_colors(1)) == list(spectrum_colors([0]))


def test_spectrum_cmap_in_matplotlib_get_cmap():
    """Test that the "spectrum" colormap is added to Matplotlib's colormaps"""
    assert matplotlib.colormaps["spectrum"] == spectrum_cmap


def test_rgb_sweep_on_ints():
    """Test the number of elements yielded from rgb_sweep"""
    colors = rgb_sweep(10, source=(1, 0, 0), dest=(0, 1, 0))
    assert len(list(colors)) == 10


def test_rgb_sweepsource_equal_to_dest():
    """Test that rgb_sweep still works if source is equal to destination."""
    colors = rgb_sweep(10, source=(1, 0, 0), dest=(1, 0, 0))
    assert len(list(colors)) == 10


def test_rgb_sweephex_colors():
    """Test that rgb_sweep works on hexadecimal strings"""
    colors = list(rgb_sweep(10, source="#ffffff", dest="#000000"))

    assert colors[0] == (1.0, 1.0, 1.0)
    assert colors[-1] == (0.0, 0.0, 0.0)
