# -*- coding: utf-8 -*-
import os

from datetime import datetime
from pathlib import Path

import numpy as np
from skimage.io import imsave
import tempfile

from skued import diffread, dmread, imibread, mibheader, mibread
from skued.utils import suppress_warnings

TEST_MIB = Path(__file__).parent / "data" / "test.mib"
TEST_MIB_MULTI = Path(__file__).parent / "data" / "test_multi.mib"

TEST_DM3 = Path(__file__).parent / "data" / "bismuth.dm3"
TEST_DM4 = Path(__file__).parent / "data" / "bismuth.dm4"

TEST_PNG = Path(__file__).parent / "data" / "png_test.png"


def test_diffread_on_merlin_image_binary():
    """Test diffread() on Merlin Image Binary (.mib)"""
    im = diffread(TEST_MIB)
    assert im.shape == (256, 256)
    # Sometimes the datatype of numpy.dstack changes slightly
    # See https://github.com/numpy/numpy/issues/21914
    assert im.dtype in {np.dtype("uint16"), np.dtype(">u2")}


def test_diffread_on_dm3_vs_dm4_image():
    """Test that diffread() works on DM3 images"""
    im3 = diffread(TEST_DM3)
    im4 = diffread(TEST_DM4)

    assert im3.shape == (256, 256)
    assert im3.dtype == np.dtype("int8")

    assert im4.shape == (256, 256)
    assert im4.dtype == np.dtype("int8")

    assert np.allclose(im3, im4)


def test_diffread_on_tiff():
    """Test diffread() on tiff files"""
    im = np.random.randint(0, 127, size=(512, 512))
    path = Path(".\\test_tif.tif")

    # Annoying low contrast warning
    with suppress_warnings():
        imsave(str(path), im)

    from_skued = diffread(path)
    assert np.allclose(im, from_skued)
    os.remove(path)


def test_diffread_on_npy():
    """Test diffread() on *.npy files"""
    arr = np.random.random(size=(128, 128))
    with tempfile.TemporaryDirectory() as tmpdir:
        np.save(Path(tmpdir) / "skued_test.npy", arr)
        arr2 = diffread(Path(tmpdir) / "skued_test.npy")
    assert np.allclose(arr, arr2)


def test_diffread_on_skimage_png():
    """Test the last resort of using skimage.io for pngs"""
    from_skimage = diffread(TEST_PNG)

    assert from_skimage.shape == (256, 256)
    assert np.allclose(from_skimage, np.ones_like(from_skimage))


def test_mib_header():
    """Test that header parsing of MIB files is working as intended"""
    header = mibheader(TEST_MIB)

    true_value = {
        "ID": "MQ1",
        "seq_num": 1,
        "offset": 384,
        "nchips": 1,
        "shape": (256, 256),
        "dtype": np.dtype(">u2"),
        "timestamp": datetime(2018, 1, 19, 20, 55, 10, 966026).timestamp(),
    }

    assert header == true_value


def test_imibread():
    """Test the generator version of mibread()"""
    gen = imibread(TEST_MIB)
    arr = next(gen)
    assert arr.shape == (256, 256)
    assert arr.dtype == np.dtype(">u2")


def test_mibread():
    """Test that the array extracted from a test MIB files has the
    expected attributes"""
    arr = mibread(TEST_MIB)
    assert arr.shape == (256, 256)
    # Sometimes the datatype of numpy.dstack changes slightly
    # See https://github.com/numpy/numpy/issues/21914
    assert arr.dtype in {np.dtype("uint16"), np.dtype(">u2")}


def test_mibread_multi():
    """Test that the array extracted from a test MIB file containing
    multiple images has the expected attributes"""
    arr = mibread(TEST_MIB_MULTI)
    assert arr.shape == (256, 256, 9)
    assert arr.dtype == np.dtype(">u1")
