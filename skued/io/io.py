# -*- coding: utf-8 -*-

from skimage.io import imread as skread
from tifffile import imread as tiffread

from .merlin import mibread

try:
    import pyqtgraph as pg
except ImportError:
    WITH_PYQTGRAPH = False
else:
    WITH_PYQTGRAPH = True
    import os
    os.environ['PYQTGRAPH_QT_LIB'] = 'PyQt5'
    from PyQt5.QtWidgets import QApplication

def diffread(fname):
    """
    Load an image from a file. This function is a generalization
    of scikit-image's ``imread`` function. ``Path`` objects are also supported.

    Parameters
    ----------
    fname : str or Path
        Image file name, e.g. ``test.mib`` or ``test.jpg``. 
    
    Returns
    -------
    img_array : `~numpy.ndarray`, ndim 2
        Diffraction image. Color images are flattened to
        grayscale.
    
    See Also
    --------
    skimage.io.imread : load images from files
    
    Notes
    -----
    Supported file formats are the same as Scikit-image, with
    the inclusion of Merlin Image Binary (.mib).
    """
    fname = str(fname)  # In case of pathlib.Path

    if fname.endswith(('tiff', 'tif')):
        return tiffread(fname)
    elif fname.endswith('.mib'):
        return mibread(fname)
    else:
        return skread(fname, as_grey = True)

def diffshow(image):
    """ 
    Display an image (from an array or from a file) in an interactive window.

    This function requires PyQtGraph and PyQt5 to be importable. These
    are optional dependencies.

    Parameters
    ----------
    image : str or array-like
        Image file name or array-like. 
    
    Raises
    ------
    ImportError : if PyQtGraph \ PyQt5 are not available

    Notes
    -----
    All file formats supported by ``skued.diffread`` are
    also supported by this function. 
    """
    if not WITH_PYQTGRAPH:
        raise ImportError('PyQtGraph is not installed')

    if isinstance(image, str):
        image = diffread(image)

    app = QApplication([])
    viewer = pg.ImageView()
    viewer.setImage(image)
    viewer.setWindowTitle('Scikit-UED Diffraction Viewer')
    viewer.show()
    app.exec_()