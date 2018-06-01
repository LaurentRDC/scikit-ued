# -*- coding: utf-8 -*-

from pathlib import Path
import skimage
import tifffile

from .merlin import mibread

try:
    import pyqtgraph as pg
except ImportError:
    WITH_PYQTGRAPH = False
else:
    WITH_PYQTGRAPH = True

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
        return tiffile.imread(fname)
    elif fname.endswith('.mib'):
        return mibread(fname)
    else:
        return skimage.io.imread(fname, as_grey = True)

def diffshow(image):
    """ 
    Display an image (from an array or from a file) in an interactive window.

    This function requires PyQtGraph to be importable. These
    are optional dependencies.

    Parameters
    ----------
    image : str or array-like
        Image file name or array-like. 
    
    Raises
    ------
    ImportError : if PyQtGraph is not available.

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
