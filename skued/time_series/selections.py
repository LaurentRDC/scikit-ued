# -*- coding: utf-8 -*-
"""
Selection masks for assembling time-series.
"""

from abc import ABCMeta, abstractmethod

import matplotlib.patches as mpatches
import numpy as np


class Selection(metaclass=ABCMeta):
    """
    Abstract base class for time-series selection masks.

    In the context of ultrafast electron/x-ray scattering, time-series are
    assembled by integrating over a portion of scattering patterns for each
    time-delay. This class is the generalization of selecting a rectangular
    area of scattering patterns to arbitrary patterns, e.g. disks, torii, etc.

    .. versionadded:: 2.0.2

    Parameters
    ----------
    shape : 2-tuple
        Shape of scattering patterns on which the selection will be applied.

    See also
    --------
    RectSelection : rectangular selection
    DiskSelection : circular disk selection
    RingSelection : ring selection, i.e. 2-torus
    ArbitrarySelection : arbitrary selection
    """

    def __init__(self, shape, *args, **kwargs):
        self.shape = shape

    @abstractmethod
    def __array__(self, *args, **kwargs):
        """Cast as a NumPy array."""
        pass

    def mpatch(self, *args, **kwargs):
        """
        Matplotlib patch associated with this selection.
        keyword arguments are passed to the appropriate `Matplotlib.patches.Patch`
        subclass.

        By default, a patch drawing a rectangle around the bounding box is used.

        .. versionadded:: 2.0.3
        """
        top, bottom, left, right = self.bounding_box

        return mpatches.Rectangle(
            xy=(left, top), width=right - left, height=bottom - top, angle=0, **kwargs
        )

    # The method below should be specialized for subclasses.
    @property
    def bounding_box(self):
        """
        Returns the array bounding box.

        Returns
        -------
        r1, r2 : int
            Row-wise bounds
        c1, c2 : int
            Column-wise bounds
        """
        selection = self.__array__()

        # The idea is to add values along rows and columns.
        # Since True ~ 1, and False ~ 0, we can determine
        # the bounding box by finding minimum and maximum indices
        # where the projection is nonzero
        # Warning : nonzero() returns a tuple!
        (row_nonzero,) = np.sum(selection, axis=1).nonzero()
        (col_nonzero,) = np.sum(selection, axis=0).nonzero()

        # In case of empty selection (i.e. all False), min and max will fail
        # Therefore, we replace with trivial value
        if row_nonzero.size == 0:
            row_nonzero = np.array([0])
        if col_nonzero.size == 0:
            col_nonzero = np.array([0])

        r1, r2 = np.min(row_nonzero), np.max(row_nonzero) + 1
        c1, c2 = np.min(col_nonzero), np.max(col_nonzero) + 1
        return (r1, r2, c1, c2)


class ArbitrarySelection(Selection):
    """
    Arbirary selection mask, represented by a boolean array.

    .. versionadded:: 2.0.2

    Parameters
    ----------
    array : ndarray, ndim 2, dtype bool
        Boolean array that evaluates to `True` on valid selections.
    """

    def __init__(self, array):
        self._array = np.asarray(array, dtype=bool)
        super().__init__(shape=self._array.shape)

    def __array__(self, *args, **kwargs):
        return self._array


class RectSelection(Selection):
    """
    Rectangular selection mask. Note that rectangular bounds are *inclusive*,
    contrary to normal numpy index selections.

    .. versionadded:: 2.0.2

    Parameters
    ----------
    shape : 2-tuple
        Shape of the scattering patterns from which data will be selected.
    r1, r2 : int
        Row indices that bound the selection, inclusively.
    c1, c2 : int
        Column indices that bound the selection, inclusively.
    """

    def __init__(self, shape, r1, r2, c1, c2):
        super().__init__(shape=shape)
        self._bbox = (r1, r2 + 1, c1, c2 + 1)

    @property
    def bounding_box(self):
        """
        Returns the array bounding box.

        Returns
        -------
        r1, r2 : int
            Row-wise bounds
        c1, c2 : int
            Column-wise bounds
        """
        return self._bbox

    def __array__(self, *args, **kwargs):
        arr = np.zeros(shape=self.shape, dtype=bool)
        r1, r2, c1, c2 = self._bbox
        arr[r1:r2, c1:c2] = True
        return arr


class DiskSelection(Selection):
    """
    Disk selection mask.

    .. versionadded:: 2.0.2

    Parameters
    ----------
    shape : 2-tuple
        Shape of the scattering patterns from which data will be selected.
    center : 2-tuple of ints
        Center (row, col) of the selection.
    radius : float
        Radius of the selection.
    """

    def __init__(self, shape, center, radius):
        super().__init__(shape=shape)
        self._center = center
        self._radius = radius

    @property
    def bounding_box(self):
        """
        Returns the array bounding box.

        Returns
        -------
        r1, r2 : int
            Row-wise bounds
        c1, c2 : int
            Column-wise bounds
        """
        rc, cc = self._center
        return (
            rc - self._radius,
            rc + self._radius + 1,
            cc - self._radius,
            cc + self._radius + 1,
        )

    def __array__(self, *args, **kwargs):
        center_row, center_col = self._center
        selection = np.zeros(shape=self.shape, dtype=bool)
        cc, rr = np.meshgrid(
            np.arange(0, self.shape[0], dtype=int) - center_col,
            np.arange(0, self.shape[1], dtype=int) - center_row,
        )
        distance = np.sqrt(rr**2 + cc**2)
        selection[distance <= self._radius] = True
        return selection

    def mpatch(self, **kwargs):
        """
        Circular patch. Keyword arguments are passed
        to `matplotlib.patches.Circle`.

        .. versionadded:: 2.0.3

        Returns
        -------
        patch : matplotlib.patches.Circle
        """
        y, x = self._center

        return mpatches.Circle(xy=(x, y), radius=self._radius, angle=0, **kwargs)


class RingSelection(Selection):
    """
    Ring selection mask, i.e. 2-torus.

    .. versionadded:: 2.0.2

    Parameters
    ----------
    shape : 2-tuple
        Shape of the scattering patterns from which data will be selected.
    center : 2-tuple of ints
        Center (row, col) of the selection.
    inner_radius : float
        Inner radius of the selection.
    outer_radius : float
        Outer radius of the selection.
    """

    def __init__(self, shape, center, inner_radius, outer_radius):
        if inner_radius > outer_radius:
            raise ValueError("Inner radius cannot be larger than outer radius.")

        super().__init__(shape=shape)
        self._center = center
        self._inner_radius = inner_radius
        self._outer_radius = outer_radius

    @property
    def bounding_box(self):
        """
        Returns the array bounding box.

        Returns
        -------
        r1, r2 : int
            Row-wise bounds
        c1, c2 : int
            Column-wise bounds
        """
        rc, cc = self._center
        return (
            rc - self._outer_radius,
            rc + self._outer_radius + 1,
            cc - self._outer_radius,
            cc + self._outer_radius + 1,
        )

    def __array__(self, *args, **kwargs):
        center_row, center_col = self._center
        selection = np.zeros(shape=self.shape, dtype=bool)
        cc, rr = np.meshgrid(
            np.arange(0, self.shape[0], dtype=int) - center_col,
            np.arange(0, self.shape[1], dtype=int) - center_row,
        )
        distance = np.sqrt(rr**2 + cc**2)
        selection[
            np.logical_and(
                distance >= self._inner_radius, distance <= self._outer_radius
            )
        ] = True
        return selection

    # TODO: make new patch class
    def mpatch(self, **kwargs):
        """
        Toroidal patch. Keyword arguments are passed
        to `matplotlib.patches.Circle`.

        .. versionadded:: 2.0.3

        Returns
        -------
        inner : matplotlib.patches.Circle
        outer : matplotlib.patches.Circle
        """
        y, x = self._center

        inner_circ = mpatches.Circle(xy=(x, y), radius=self._inner_radius, **kwargs)
        outer_circ = mpatches.Circle(xy=(x, y), radius=self._outer_radius, **kwargs)

        return inner_circ, outer_circ


class RingArcSelection(Selection):
    """
    Selection patch for a partial 2-torus.

    .. versionadded:: 2.0.5

    Parameters
    ----------
    shape : 2-tuple
        Shape of the scattering patterns from which data will be selected.
    center : 2-tuple of ints
        Center (row, col) of the selection.
    inner_radius : float
        Inner radius of the selection.
    outer_radius : float
        Outer radius of the selection.
    angle : float
        Rotation of the ring in degrees.
    theta1, theta2 : float
        Starting and ending angles of the 2-torus in degrees, relative to ``angle``.
    """

    def __init__(
        self, shape, center, inner_radius, outer_radius, angle=0, theta1=0, theta2=360
    ):
        if inner_radius > outer_radius:
            raise ValueError("Inner radius cannot be larger than outer radius.")

        super().__init__(shape=shape)
        self._center = center
        self._inner_radius = inner_radius
        self._outer_radius = outer_radius
        self._angle = angle
        self._theta1 = theta1
        self._theta2 = theta2

    @property
    def bounding_box(self):
        """
        Returns the array bounding box.

        Returns
        -------
        r1, r2 : int
            Row-wise bounds
        c1, c2 : int
            Column-wise bounds
        """
        rc, cc = self._center
        return (
            rc - self._outer_radius,
            rc + self._outer_radius + 1,
            cc - self._outer_radius,
            cc + self._outer_radius + 1,
        )

    def __array__(self, *args, **kwargs):
        center_row, center_col = self._center
        selection = np.zeros(shape=self.shape, dtype=bool)
        cc, rr = np.meshgrid(
            np.arange(0, self.shape[0], dtype=int) - center_col,
            np.arange(0, self.shape[1], dtype=int) - center_row,
        )
        distance = np.sqrt(rr**2 + cc**2)
        angle = np.rad2deg(np.arctan2(rr, cc)) + self._angle
        angle[:] = np.mod(angle, 360)

        distance_criteria = np.logical_and(
            distance >= self._inner_radius, distance <= self._outer_radius
        )
        angle_criteria = np.logical_and(angle >= self._theta1, angle <= self._theta2)
        selection[np.logical_and(angle_criteria, distance_criteria)] = True
        return selection

    def mpatch(self, **kwargs):
        """
        Partial toroidal patch. Keyword arguments are passed
        to `matplotlib.patches.Arc`.

        Returns
        -------
        inner : matplotlib.patches.Circle
        outer : matplotlib.patches.Circle
        """
        y, x = self._center

        arc = lambda radius: mpatches.Arc(
            xy=(x, y),
            width=2 * radius,
            height=2 * radius,
            angle=self._angle,
            theta1=self._theta1,
            theta2=self._theta2,
            **kwargs
        )

        inner_arc = arc(self._inner_radius)
        outer_arc = arc(self._outer_radius)

        return inner_circ, outer_circ
