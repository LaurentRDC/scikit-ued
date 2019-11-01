# -*- coding: utf-8 -*-
from contextlib import contextmanager
from pathlib import Path

from .io import diffread

try:
    import pyqtgraph as pg
except ImportError:
    WITH_PYQTGRAPH = False
else:
    WITH_PYQTGRAPH = True


@contextmanager
def rowmajor_axisorder():
    """
    Context manager that sets the PyQtGraph image axis order to row-major. 
    The environment is reset to the initial value after context close.
    """
    old_image_axis_order = pg.getConfigOption("imageAxisOrder")
    pg.setConfigOptions(imageAxisOrder="row-major")
    yield
    pg.setConfigOptions(imageAxisOrder=old_image_axis_order)


def diffshow(image):
    """ 
    Display an image (from an array or from a file) in an interactive window.

    This function requires `PyQtGraph` to be importable. These
    are optional dependencies.

    Parameters
    ----------
    image : path-like or array-like
        Image file name or array-like. All file formats supported 
        by ``skued.diffread`` are also supported by this function. 
    
    Raises
    ------
    ImportError : if `PyQtGraph` is not available.
    """
    if not WITH_PYQTGRAPH:
        raise ImportError("PyQtGraph is not installed.")

    if isinstance(image, (str, Path)):
        image = diffread(image)

    with rowmajor_axisorder():
        app = pg.QtGui.QApplication([])
        viewer = pg.ImageView()
        viewer.setImage(image)
        viewer.setWindowTitle("scikit-ued image viewer")
        viewer.show()
        app.exec_()
