# -*- coding: utf-8 -*-

from pathlib import Path

import tifffile

from .dm import dmread
from .merlin import mibread

try:
    import pyqtgraph as pg
except ImportError:
    WITH_PYQTGRAPH = False
else:
    WITH_PYQTGRAPH = True

def diffread(fname):
    """
    Load an image from a file. Supported file formats are:

        * Merlin Image Binary (`*.mib`)

        * Digital Micrograph (`*.dm3`, `*.dm4`) 
        
            .. versionadded:: 1.0.1.0

        * TIFF (`*.tif`, `*.tiff`)

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
    fname = str(fname)  # In case of pathlib.Path

    if fname.endswith(('tiff', 'tif')):
        return tifffile.imread(fname)
    elif fname.endswith('.mib'):
        return mibread(fname)
    elif fname.endswith(('.dm3', '.dm4')):
        return dmread(fname)
    
    # Worst-case scenario, we need skimage.io
    # since this is a slow import, we import it here
    import skimage.io

    return skimage.io.imread(fname, as_gray = True)

def diffshow(image):
    """ 
    Display an image (from an array or from a file) in an interactive window.

    This function requires `PyQtGraph` to be importable. These
    are optional dependencies.

    Parameters
    ----------
    image : path-like or array-like
        Image file name or array-like. 
    
    Raises
    ------
    ImportError : if `PyQtGraph` is not available.

    Notes
    -----
    All file formats supported by ``skued.diffread`` are
    also supported by this function. 
    """
    if not WITH_PYQTGRAPH:
        raise ImportError('PyQtGraph is not installed.')

    if isinstance(image, (str, Path)):
        image = diffread(image)

    app = pg.QtGui.QApplication([])
    viewer = pg.ImageView()
    viewer.setImage(image)
    viewer.setWindowTitle('Scikit-UED Diffraction Viewer')
    viewer.show()
    app.exec_()
