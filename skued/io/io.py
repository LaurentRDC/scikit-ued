# -*- coding: utf-8 -*-

from pathlib import Path

import numpy as np

from .dm import dmread
from .merlin import mibread


def diffread(fname):
    """
    Load an image from a file. Supported file formats are:

        * Merlin Image Binary (`*.mib`)

        * Digital Micrograph (`*.dm3`, `*.dm4`)

            .. versionadded:: 1.0.1.0

        * TIFF (`*.tif`, `*.tiff`)

        * NumPy (`*.npy`)

        * All file formats supported by Scikit-image

    Parameters
    ----------
    fname : path-like
        Image file name.

    Returns
    -------
    img_array : `~numpy.ndarray`, ndim 2
        Diffraction image. Color images are flattened to
        grayscale.

    See Also
    --------
    skimage.io.imread : load images from files
    skued.mibread : Read single- and multi-image Merlin Image Binary files
    """
    fname = Path(fname)

    if fname.suffix == ".mib":
        return mibread(fname)
    elif fname.suffix in {".dm3", ".dm4"}:
        return dmread(fname)
    elif fname.suffix == ".npy":
        return np.load(fname, fix_imports=False)

    # Worst-case scenario, we need skimage.io
    # since this is a slow import, we import it here
    import skimage.io

    return skimage.io.imread(fname, as_gray=True)
