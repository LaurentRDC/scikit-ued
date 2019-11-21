# -*- coding: utf-8 -*-

import unittest

import numpy as np

from skued import (
    Selection,
    RectSelection,
    DiskSelection,
    RingSelection,
    ArbitrarySelection,
)


class TestSelectionArray(unittest.TestCase):
    def test_rect_array(self):
        selection = RectSelection((8, 8), 0, 1, 2, 3)
        toarr = np.array(selection)
        self.assertEqual(toarr.shape, (8, 8))
        self.assertEqual(toarr.dtype, np.bool)

    def test_disk_array(self):
        selection = DiskSelection((8, 8), center=(4, 4), radius=2)
        toarr = np.array(selection)
        self.assertEqual(toarr.shape, (8, 8))
        self.assertEqual(toarr.dtype, np.bool)

    def test_ring_array(self):
        selection = RingSelection((8, 8), center=(4, 4), inner_radius=0, outer_radius=2)
        toarr = np.array(selection)
        self.assertEqual(toarr.shape, (8, 8))
        self.assertEqual(toarr.dtype, np.bool)

    def test_arbitrary_array(self):
        selection = ArbitrarySelection(np.random.choice([True, False], size=(8, 8)))
        toarr = np.array(selection)
        self.assertEqual(toarr.shape, (8, 8))
        self.assertEqual(toarr.dtype, np.bool)


class TestSelectionBBox(unittest.TestCase):
    """ Test the Selection.bounding_box property for all Selection subclasses """

    def test_rect_selection(self):

        with self.subTest("Trivial"):
            selection = RectSelection((8, 8), 0, 0, 0, 0)
            self.assertTupleEqual(selection.bounding_box, (0, 1, 0, 1))

        with self.subTest("Non-trivial"):
            selection = RectSelection((8, 8), 3, 7, 1, 3)
            self.assertTupleEqual(selection.bounding_box, (3, 8, 1, 4))

    def test_disk_selection(self):
        with self.subTest("Trivial"):
            selection = DiskSelection((8, 8), center=(4, 4), radius=0)
            self.assertTupleEqual(selection.bounding_box, (4, 5, 4, 5))

        with self.subTest("Non-trivial"):
            selection = DiskSelection((8, 8), center=(1, 2), radius=1)
            self.assertTupleEqual(selection.bounding_box, (0, 3, 1, 4))

    def test_ring_selection(self):

        # A ring selection should have the same bounding box as a disk selection
        # with equivalent outer radius
        disk = DiskSelection((16, 16), center=(8, 7), radius=3)
        ring = RingSelection((16, 16), center=(8, 7), inner_radius=1, outer_radius=3)

        self.assertEqual(disk.bounding_box, ring.bounding_box)

    def test_arbitrary_selection(self):

        arr = np.ones(shape=(8, 8), dtype=np.bool)  # selection everywhere
        selection = ArbitrarySelection(arr)

        self.assertEqual(selection.bounding_box, (0, 8, 0, 8))


if __name__ == "__main__":
    unittest.main()
