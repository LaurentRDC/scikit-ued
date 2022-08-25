# -*- coding: utf-8 -*-
import warnings
from contextlib import contextmanager
from pathlib import Path

from .io import diffread

try:
    import pyqtgraph as pg
except ImportError:
    WITH_PYQTGRAPH = False
else:
    WITH_PYQTGRAPH = True


if WITH_PYQTGRAPH:

    # This is weird, but PyQtGraph is an optional dependency
    # Therefore, we cannot define this class unless PyQtGraph is importable
    if hasattr(pg, "QtWidgets"):
        BASE_CLASS = pg.QtWidgets
    elif hasattr(pg, "QtGui"):
        BASE_CLASS = pg.QtGui

    class Diffshow(BASE_CLASS.QWidget):
        """
        Widget containing a main viewer, plus some cursor information.

        Parameters
        ----------
        image : ndarray
        """

        def __init__(self, image, **kwargs):
            super().__init__(**kwargs)
            self.viewer = pg.ImageView()

            with warnings.catch_warnings():
                # Pesky FutureWarning from PyQtGraph
                warnings.simplefilter("ignore")
                self.viewer.setImage(image)

            self.cursor_info = pg.QtGui.QLabel("")
            self.cursor_info.setAlignment(pg.QtCore.Qt.AlignCenter)

            self.__cursor_proxy = pg.SignalProxy(
                self.viewer.scene.sigMouseMoved,
                rateLimit=60,
                slot=self.update_cursor_info,
            )

            self.setWindowTitle("scikit-ued image viewer")

            layout = pg.QtGui.QVBoxLayout()
            layout.addWidget(self.viewer)
            layout.addWidget(self.cursor_info)
            self.setLayout(layout)

        def update_cursor_info(self, event):
            """Determine cursor information from mouse event."""
            mouse_point = self.viewer.getView().mapSceneToView(event[0])
            i, j = int(mouse_point.y()), int(mouse_point.x())
            try:
                val = self.viewer.getImageItem().image[i, j]
            except IndexError:
                val = 0
            self.cursor_info.setText(
                f"Position: ({i},{j}) | Pixel value: {val:.2f} cnts"
            )


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
        viewer = Diffshow(image)
        viewer.show()
        app.exec_()
