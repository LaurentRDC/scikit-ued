# -*- coding: utf-8 -*-


import numpy as np

from skued import (
    RectSelection,
    DiskSelection,
    RingSelection,
    RingArcSelection,
    ArbitrarySelection,
)


def test_selection_rect_array():
    selection = RectSelection((8, 8), 0, 1, 2, 3)
    toarr = np.array(selection)
    assert toarr.shape == (8, 8)
    assert toarr.dtype == bool


def test_selection_disk_array():
    selection = DiskSelection((8, 8), center=(4, 4), radius=2)
    toarr = np.array(selection)
    assert toarr.shape == (8, 8)
    assert toarr.dtype == bool


def test_selection_ring_array():
    selection = RingSelection((8, 8), center=(4, 4), inner_radius=0, outer_radius=2)
    toarr = np.array(selection)
    assert toarr.shape == (8, 8)
    assert toarr.dtype == bool


def test_selection_ring_arc_array():
    fullarc = RingArcSelection(
        (8, 8),
        center=(4, 3),
        inner_radius=1,
        outer_radius=3,
        angle=0,
        theta1=0,
        theta2=360,
    )
    toarr = np.array(fullarc)
    assert toarr.shape == (8, 8)
    assert toarr.dtype == bool

    ring = RingSelection((8, 8), center=(4, 3), inner_radius=1, outer_radius=3)

    assert np.allclose(fullarc, ring)


def test_selection_arbitrary_array():
    selection = ArbitrarySelection(np.random.choice([True, False], size=(8, 8)))
    toarr = np.array(selection)
    assert toarr.shape == (8, 8)
    assert toarr.dtype == bool


def test_selection_bbox_rect_selection():

    selection = RectSelection((8, 8), 0, 0, 0, 0)
    assert selection.bounding_box == (0, 1, 0, 1)

    selection = RectSelection((8, 8), 3, 7, 1, 3)
    assert selection.bounding_box == (3, 8, 1, 4)


def test_selection_bbox_disk_selection():
    selection = DiskSelection((8, 8), center=(4, 4), radius=0)
    assert selection.bounding_box == (4, 5, 4, 5)

    selection = DiskSelection((8, 8), center=(1, 2), radius=1)
    assert selection.bounding_box == (0, 3, 1, 4)


def test_selection_bbox_ring_selection():

    # A ring selection should have the same bounding box as a disk selection
    # with equivalent outer radius
    disk = DiskSelection((16, 16), center=(8, 7), radius=3)
    ring = RingSelection((16, 16), center=(8, 7), inner_radius=1, outer_radius=3)

    assert disk.bounding_box == ring.bounding_box


def test_selection_bbox_ring_arc_selection():
    ring = RingSelection((16, 16), center=(8, 7), inner_radius=1, outer_radius=3)
    fullarc = RingArcSelection(
        (16, 16),
        center=(8, 7),
        inner_radius=1,
        outer_radius=3,
        angle=0,
        theta1=0,
        theta2=360,
    )
    assert ring.bounding_box == fullarc.bounding_box


def test_selection_bbox_arbitrary_selection():

    arr = np.ones(shape=(8, 8), dtype=bool)  # selection everywhere
    selection = ArbitrarySelection(arr)

    assert selection.bounding_box == (0, 8, 0, 8)
