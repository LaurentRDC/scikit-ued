# -*- coding: utf-8 -*-
"""
Module concerned with determination of BZs
=====================================================
"""
import numpy as np
from scipy.spatial import Voronoi
import sys
from .center import autocenter
from ..time_series import DiskSelection


class VoronoiRegion:
    """
    This class creates a region of points defined by the vertices of the Voronoi regions of a given pattern.
    It is made to make calculations later easier to manage

    Parameters
    ----------

    region_id: type(int)
        Defines the ID number of this particular Voronoi region
    """

    def __init__(self, region_id):
        self.id = region_id
        self.vertices = []
        self.is_inf = False
        self.point_inside = None
        self.visible_vertices = []
        self.is_visible = False

    def __str__(self):
        text = f"region id={self.id}"
        if self.point_inside:
            point_idx, point = self.point_inside
            text = f"{text}[point:{point}(point_id:{point_idx})]"
        text += ", vertices: "
        if self.is_inf:
            text += "(inf)"
        for v in self.vertices:
            text += f"{v}"
        return text

    def __repr__(self):
        return str(self)

    def add_vertex(self, vertex, vertices):
        if vertex == -1:
            self.is_inf = True
        else:
            point = vertices[vertex]
            self.vertices.append(point)

    def add_visible_vertex(self, vertex):
        self.visible_vertices.append(vertex)

    def set_center(self, center):
        self.center = center


class brillouin_zones:
    def __init__(self, image, mask, peaks, center=None, optimization_radius=None):
        """
        Generate Brillouin zone projections in the particular 2D geometry based on
        Bragg peak locations.

        Parameters
        ----------
        image : `~numpy.ndarray`, shape (M,N)
            Image to be aligned.
        mask : `~numpy.ndarray`, shape (M,N)
            Mask that evaluates to True on valid pixels of the array `image`.
        center : `~numpy.ndarray`, shape (2,), optional
            center of the image. Else, calls autocentering functionality
        optimization_radius : float
            Number of pixels around each Bragg peak to perform self-optimization.
            If not specified, auto-optimization skipped.
        """
        self._image = image
        self._mask = mask
        self.visible_BZs = None
        try:
            peaks = peaks.reshape(-1, 2)
        except ValueError:
            print("Peaks seem to have more than 2 dimensions. Truncating...")
            peaks = peaks[:, :2]
            try:
                peaks = peaks.reshape(-1, 2)
            except ValueError as e:
                print(f"Error loading peaks: {e}")
                sys.exit(2)
        if center is None:
            self.center = autocenter(image, mask)
        else:
            self.center = center
        peaks = np.vstack((self.center, peaks))  # self.center[::-1]
        new_peaks = np.array(peaks)
        if optimization_radius is not None:
            skipped_peaks = list()
            # assert type(optimization_radius) == float
            for idx, peak in enumerate(peaks):
                if (
                    idx > 0
                ):  # do not optimize center, which is masked anyway most likely
                    r, c = peaks[idx]
                    peak = np.asarray(peak).astype(int)
                    disk = DiskSelection(
                        shape=self._image.shape,
                        center=peak,
                        radius=optimization_radius,
                    )
                    r1, r2, c1, c2 = tuple([int(d) for d in disk.bounding_box])
                    region = image[r1:r2, c1:c2]
                    try:
                        true_peak_idx_local = np.where(region == region.max())
                        true_peak_idx_global = np.where(
                            self._image == region[true_peak_idx_local]
                        )
                        new_r, new_c = (
                            true_peak_idx_global[0][0],
                            true_peak_idx_global[1][0],
                        )
                    except:
                        # if idx != 0:
                        skipped_peaks.append(idx)
                        new_r, new_c = r, c
                    new_peaks[idx] = np.array((new_r, new_c)).astype(int)
            print(f"Could not auto-adjust peak #s: {skipped_peaks}.")
        self.bragg_peaks = new_peaks

        self.__vor = Voronoi(new_peaks)  # Voronoi(new_peaks[:, [1, 0]])
        self.__voronoi_regions = []

        for i, point_region in enumerate(self.__vor.point_region):
            region = self.__vor.regions[point_region]
            vr = VoronoiRegion(point_region)
            for r in region:
                vr.add_vertex(r, self.__vor.vertices)
            vr.point_inside = (i, self.__vor.points[i])
            vr.set_center(self.__vor.points[i])
            self.__voronoi_regions.append(vr)
        return

    def getVisibleBZs(self, symmetry=None, bbox=None):
        """
        Determines which BZs have the following criteria:
            (i) Does the BZ have the number of corners given by `symmetry`?
            (ii) If bbox, does the BZ lie within the `bbox` width?

        Parameters
        ----------

        symmetry : int, optional
            Specifies number of vertices each BZ must have

        bbox : float, optional
            Further require that the BZ be within `bbox` from the center of the image

        Returns
        -------
        visible_BZs : list of :class: VoronoiRegions
        """
        if bbox is None:
            bbox = int(0.49 * self._image.shape[0])
        else:
            assert type(bbox) == int
        self._symmetry = symmetry
        for r in self.__voronoi_regions:
            if not r.is_inf:
                verts = np.array(r.vertices).reshape(-1, 2)
                COND1 = np.all(
                    np.sqrt(np.sum((verts - self.center) ** 2, axis=1)) < bbox
                )
                COND2 = (
                    (verts.shape[0] == self._symmetry)
                    if self._symmetry is not None
                    else True
                )
                if COND1 and COND2:
                    r.add_visible_vertex(verts)
            r.visible_vertices = np.array(r.visible_vertices).reshape(-1, 2)
            if r.visible_vertices.size != 0:
                r.is_visible = True
        self.visible_BZs = [r for r in self.__voronoi_regions if r.is_visible]
        print(f"Determined BZ consistency: {1e2*self.determineConsistency():.2f}%")
        return self.visible_BZs

    def renderVisibleBZs(self, axis, BZs=None, **kwargs):
        """
        Renders BZs (specified manually as argument or all visible ones determined by a previous run
        of :func: `getVisibleBZs`) on a given axis. Passes kwargs to plotting routine.

        Parameters
        ----------

        axis : matplotlib `axis` object
            Specifies the axes onto which the BZs will be rendered

        BZs : list of :class: VoronoiRegions, optional
            List of BZ objects to plot. Defaults to all visible ones if not provided.

        """
        from matplotlib.patches import Polygon

        for r in self.visible_BZs if BZs is None else BZs:
            axis.add_patch(
                Polygon(np.array(r.vertices).reshape(-1, 2), fill=False, **kwargs)
            )

    def determineConsistency(self, BZs=None):
        """
        Determines the consistency of the BZ determination. Simply the
        ratio of the standard deviation of BZ areas over the mean of the areas.
        Areas are determined through a simple implementation of the Shoelace formula.

        Parameters
        ----------
        BZs : list of :class: VoronoiRegions
            Iterable of regions who will be considered for the consistency determination
            calculation

        Returns
        -------
        result : float
            Normalized consistency of the BZ generation

        References
        ----------
        https://en.wikipedia.org/wiki/Shoelace_formula
        """

        regions = self.visible_BZs if BZs is None else BZs
        self._areas = [
            0.5
            * np.abs(
                np.dot(r.visible_vertices[:, 0], np.roll(r.visible_vertices[:, 1], 1))
                - np.dot(r.visible_vertices[:, 1], np.roll(r.visible_vertices[:, 0], 1))
            )
            for r in regions
        ]
        return 1 - np.std(self._areas) / np.average(self._areas)

    def getEquivalentScatteringVectors(self, Qvector, use_visible=False):
        """
        For a given scattering vector `Q', determine all equivalent reduced scattering vectors
        based on the Bragg peaks determined to be visible by the methods of this class.

        Parameters
        ----------

        Qvector : `~numpy.ndarray`, shape (2,)
            The scattering vector to consider. The closest Bragg peak will be determined,
            then the relative angle to this peak, and then it will iterate over the other peaks
            and determine the same vectors who have radial and angular displacements identical
            to this vector relative to the other Bragg peaks.

        use_visible : bool
            If True, consider only Bragg peaks inside of visible BZs determined by
            a previous run of :meth:`getVisibleBZs()`. Else, use all Bragg peaks
            used as input for the initialization of this class.

        Returns
        -------

        QVectors : `~numpy.ndarray`, shape (N,2)
            The determined equivalent scattering vectors. The number of vectors
            is dependent on the set of peaks used determined by the `use_visible` parameter.
        """
        self.__Qvector = (
            Qvector - self.center
        )  # place scat vec in pixel space with shifted origin to (0,0)
        self.__bragg_peaks = (
            self.bragg_peaks - self.center
        )  # place undiffracted beam at origin in pixel space

        deltas = self.__bragg_peaks - self.__Qvector
        CLOSEST_PEAK_INDEX = np.argmin(np.einsum("ij,ij->i", deltas, deltas))
        REDUCED_WAVEVECTOR = (
            self.__Qvector - self.__bragg_peaks[CLOSEST_PEAK_INDEX]
        )  # this gives wavevector in pixel space with origin at (0,0)

        # Now determine relative angle of reduced wavevector to the closest bragg peak
        inner = np.inner(REDUCED_WAVEVECTOR, self.__bragg_peaks[CLOSEST_PEAK_INDEX])
        norms = np.linalg.norm(REDUCED_WAVEVECTOR) * np.linalg.norm(
            self.__bragg_peaks[CLOSEST_PEAK_INDEX]
        )
        ORIGINAL_ANG = np.arccos(np.clip(inner / norms, -1.0, 1.0))
        REDUCED_WAVEVECTOR_RADIUS = np.linalg.norm(REDUCED_WAVEVECTOR)
        ANGULAR_SEP = np.arctan2(REDUCED_WAVEVECTOR[1], REDUCED_WAVEVECTOR[0])
        QVectors = list()
        if use_visible:
            for r in self.visible_BZs:
                CURRENT_ANGLE = np.arctan2(
                    r.center[1] - self.center[1], r.center[0] - self.center[0]
                )
                COSANG = np.cos(ANGULAR_SEP + CURRENT_ANGLE - ORIGINAL_ANG)
                SINANG = np.sin(ANGULAR_SEP + CURRENT_ANGLE - ORIGINAL_ANG)
                QVectors.append(
                    r.center + REDUCED_WAVEVECTOR_RADIUS * np.array([COSANG, SINANG])
                )

        else:
            for peak in self.bragg_peaks:
                CURRENT_ANGLE = np.arctan2(
                    peak[1] - self.center[1], peak[0] - self.center[0]
                )
                COSANG = np.cos(ANGULAR_SEP + CURRENT_ANGLE - ORIGINAL_ANG)
                SINANG = np.sin(ANGULAR_SEP + CURRENT_ANGLE - ORIGINAL_ANG)
                QVectors.append(
                    peak + REDUCED_WAVEVECTOR_RADIUS * np.array([COSANG, SINANG])
                )
        QVectors = (
            np.array(sorted(QVectors, key=lambda p: np.linalg.norm(p - self.center)))
            .reshape(-1, 2)
            .astype(int)
        )
        return QVectors[1:, :]  # ignore vector around center
