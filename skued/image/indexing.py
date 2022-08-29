# -*- coding: utf-8 -*-
"""
Indexing single crystals 
------------------------
"""
from warnings import catch_warnings, simplefilter
from itertools import combinations
import numpy as np
from scipy.ndimage import gaussian_filter, laplace, shift
from skimage.registration._masked_phase_cross_correlation import cross_correlate_masked
import skimage.filters as filters
from skimage.measure import label, regionprops
from skimage.morphology import binary_erosion, disk

from ..fft import with_skued_fft
from .center import autocenter


@with_skued_fft
def bragg_peaks(im, mask=None, center=None, min_dist=None):
    """
    Extract the position of Bragg peaks in a single-crystal diffraction pattern.

    .. versionadded:: 2.1.3

    Parameters
    ----------
    im : ndarray, shape (N,M)
        Single-crystal diffraction pattern.
    mask : ndarray, shape (N,M), dtype bool, optional
        Mask that evaluates to `True` on valid pixels of ``im``.
    center : 2-tuple, optional
        Center of the diffraction pattern, in ``(row, col)`` format.
        If ``None``, the center will be determined via :func:`autocenter`.
    min_dist : float or None, optional
        Minimum distance between Bragg peaks (in pixel coordinates). Peaks that are closer
        than this distance will be considered the same peak, and only one of them will
        be returned. If `None` (default), the minimum distance is guessed based on the
        image size.

    Returns
    -------
    peaks : list of 2-tuples
        List of coordinates ``[row, col]`` for every detected peak, sorted
        in order of how close they are to the center of the image.

    References
    ----------
    Liu, Lai Chung. Chemistry in Action: Making Molecular Movies with Ultrafast
    Electron Diffraction and Data Science, Chapter 2. Springer Nature, 2020.
    """
    if mask is None:
        mask = np.ones(im.shape)
    if center is None:
        center = autocenter(im=im, mask=mask)

    im = np.array(im, copy=True, dtype=float)
    im -= im.min()

    with catch_warnings():
        simplefilter("ignore", category=RuntimeWarning)
        im /= gaussian_filter(input=im, sigma=min(im.shape) / 20, truncate=2)
    im = np.nan_to_num(im, copy=False)

    autocorr = np.abs(
        cross_correlate_masked(arr1=im, arr2=im, m1=mask, m2=mask, mode="same")
    )

    # The regions of interest are defined on the labels made
    # from the autocorrelation of the image. The center of the autocorr
    # is the center of the array; we need to correct the offset.
    # This also allows to use the mask on properties derived
    # from the autocorr
    autocorr = shift(
        autocorr,
        shift=np.asarray(center) - np.array(im.shape) / 2,
        order=1,
        mode="nearest",
    )

    laplacian = -1 * laplace(autocorr)
    threshold = filters.threshold_triangle(laplacian)
    regions = (laplacian > threshold) * mask

    # To prevent noise from looking like actual peaks,
    # we erode labels using a small selection area
    regions = binary_erosion(regions, footprint=disk(2))

    labels = label(regions, return_num=False)
    props = regionprops(label_image=labels, intensity_image=im)
    candidates = [
        prop for prop in props if not np.any(np.isnan(prop.weighted_centroid))
    ]

    # Some regions are very close to each other; we prune them!
    if min_dist is None:
        min_dist = min(im.shape) / 100

    peaks = list()
    for prop in candidates:
        pos = np.asarray(prop.weighted_centroid)
        if any((np.linalg.norm(peak - pos) < min_dist) for peak in peaks):
            continue
        else:
            peaks.append(pos)

    return sorted(peaks, key=lambda p: np.linalg.norm(p - center))


def bragg_peaks_persistence(
    im, mask=None, center=None, min_dist=None, bd_threshold=0.0, prominence=0.0
):
    """
    Extract the position of Bragg peaks in a single-crystal diffraction pattern
    using 2D persistence of the image landscape. If detected peaks fall above
    certain thresholds, they are deamed acceptable peaks and returned.

    .. versionadded:: 2.5.5

    Parameters
    ----------
    im : ndarray, shape (N,M)
        Single-crystal diffraction pattern.
    mask : ndarray, shape (N,M), dtype bool, optional
        Mask that evaluates to `True` on valid pixels of ``im``.
    center : 2-tuple, optional
        Center of the diffraction pattern, in ``(row, col)`` format.
        If ``None``, the center will be determined via :func:`autocenter`.
    min_dist : float or list (size 2) or None, optional
        Minimum distance between Bragg peaks (in pixel coordinates). Peaks that are closer
        than this distance will be considered the same peak, and only one of them will
        be returned. If passing only single number, minimum distance is determined by
        the degree of radial separation. Else, minimum distancing is determined by thresholding
        x and y coordinates separately with min_dist[0] and min_dist[1].
        If `None` (default), the minimum distance criterion is not enforced, and all peaks
        passing the prominence threshold are returned.
    bd_threshold : float, optional
        Specifies the persistence threshold for the determination of
        birth/death parameters in the persistence calculations
    prominence : float, optional
        Specifies the persistence threshold for the distinction between
        acceptable peaks and 'ghost' peaks. Lower values will yield more peaks.


    Returns
    -------
    peaks : ndarray, shape (N,2)
        Array of coordinates ``[row, col]`` for every detected peak, sorted
        in order of how close they are to the center of the image.
    birth_death : ndarray, shape (N,2)
        Array of birth, death attributes for every detected peak.
    birth_index_indices : list
        Indices of candidate peaks who pass the bd_threshold
        prominence criterion
    persistencies : list of floats
        List of persistencies of the detected peaks. Useful for determining
        which threshold to use to accurate peak determination.

    References
    ----------
    Huber, S. (2021). Persistent Homology in Data Science. In: Haber, P.,
    Lampoltshammer, T., Mayr, M., Plankensteiner, K. (eds) Data Science - Analytics
    and Applications. Springer Vieweg, Wiesbaden.
    https://doi.org/10.1007/978-3-658-32182-6_13

    Edelsbrunner, H. and John L Harer (2010). Computational Topology. In: American
    Mathematical Society.

    Example
    -------
    To generate the peak determination from a 2D array and visualize the results,
    as well as the birth-death persistence diagram:
    ```
    from skued import bragg_peak_persistence
    peaks, birth_death, indices = bragg_peak_persistence(image)
    axis.imshow(image, extent='lower')
    for i, peak in enumerate(peaks):
        x, y = peak
        ax.plot([x], [y], '.', c='b')
        ax.text(x, y+0.25, str(i+1), color='b')

    ax2.set_title("Peristence diagram")
    for i, (index, bd) in enumerate(zip(indices, birth_death)):
        x, y = bd
        ax2.plot([x], [y], '.', c='b')
        ax2.text(x, y+2, str(index+1), color='b')
    ax2.set_xlabel("Birth level")
    ax2.set_ylabel("Death level")
    ```
    """
    if mask is None:
        mask = np.ones_like(im, dtype=bool)
    if center is None:
        center = autocenter(im=im, mask=mask)

    g0 = Persistence(im).persistence
    birth_death = list()
    birth_death_indices = list()
    persistencies = list()
    candidates = list()
    for i, homclass in enumerate(g0):
        p_birth, bl, pers, p_death = homclass
        persistencies.append(pers)
        if pers <= bd_threshold:
            continue
        x, y = bl, bl - pers
        birth_death.append([x, y])
        birth_death_indices.append(i)
    for i, homclass in enumerate(g0):
        p_birth, bl, pers, p_death = homclass
        if pers <= prominence:
            continue
        y, x = p_birth
        candidates.append([x, y])
    if min_dist is not None:
        combos = combinations(candidates, 2)
        if type(min_dist) == int or type(min_dist) == float:
            points_to_remove = [
                point2
                for point1, point2 in combos
                if np.linalg.norm(np.array(point1) - np.array(point2)) < min_dist
            ]

        else:
            points_to_remove = [
                point2
                for point1, point2 in combos
                if abs(point1[0] - point2[0]) <= min_dist[0]
                and abs(point1[1] - point2[1]) <= min_dist[1]
            ]
        candidates = [point for point in candidates if point not in points_to_remove]

    candidates = np.array(candidates).reshape(-1, 2)
    peaks = np.array(
        sorted(candidates, key=lambda p: np.linalg.norm(p - center))
    ).reshape(-1, 2)
    birth_death = np.array(birth_death).reshape(-1, 2)

    # remove peaks that are within the masked area
    if mask.sum() != mask.shape[0] * mask.shape[1]:
        peaks = np.array([p for p in peaks if mask[p[1], p[0]]])
        birth_death = np.array(
            [bd for p, bd in zip(peaks, birth_death) if mask[p[1], p[0]]]
        )
        birth_death_indices = np.array(
            [bdi for p, bdi in zip(peaks, birth_death_indices) if mask[p[1], p[0]]]
        )
        persistencies = np.array(
            [pers for p, pers in zip(peaks, persistencies) if mask[p[1], p[0]]]
        )
    return peaks, birth_death, birth_death_indices, persistencies


class UnionFind:

    """Union-find data structure.

    Each unionFind instance X maintains a family of disjoint sets of
    hashable objects, supporting the following two methods:

    - X[item] returns a name for the set containing the given item.
      Each set is named by an arbitrarily-chosen one of its members; as
      long as the set remains unchanged it will keep the same name. If
      the item is not yet part of a set in X, a new singleton set is
      created for it.

    - X.union(item1, item2, ...) merges the sets containing each item
      into a single larger set.  If any item is not yet part of a set
      in X, it is added to X as one of the members of the merged set.
    """

    def __init__(self):
        """Create a new empty union-find structure."""
        self.weights = {}
        self.parents = {}

    def add(self, object, weight):
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = weight

    def __contains__(self, object):
        return object in self.parents

    def __getitem__(self, object):
        """Find and return the name of the set containing the object."""

        # check for previously unknown object
        if object not in self.parents:
            assert False
            self.parents[object] = object
            self.weights[object] = 1
            return object

        # find path of objects leading to the root
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root

    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure."""
        return iter(self.parents)

    def union(self, *objects):
        """Find the sets containing the objects and merge them all."""
        roots = [self[x] for x in objects]
        heaviest = max([(self.weights[r], r) for r in roots])[1]
        for r in roots:
            if r != heaviest:
                self.parents[r] = heaviest


class Persistence:
    def __init__(self, im):
        self.image = im
        self.calculate()

    def calculate(self):
        h, w = self.image.shape

        # Get indices orderd by value from high to low
        indices = [(i, j) for i in range(h) for j in range(w)]
        indices.sort(key=lambda p: self.get(p), reverse=True)

        # Maintains the growing sets
        self.uf = UnionFind()

        self._groups0 = {}

        # Process pixels from high to low
        for i, p in enumerate(indices):
            v = self.get(p)
            ni = [self.uf[q] for q in self.iter_neighbors(p, w, h) if q in self.uf]
            nc = sorted([(self.get_comp_birth(q), q) for q in set(ni)], reverse=True)

            if i == 0:
                self._groups0[p] = (v, v, None)

            self.uf.add(p, -i)

            if len(nc) > 0:
                oldp = nc[0][1]
                self.uf.union(oldp, p)

                # Merge all others with oldp
                for bl, q in nc[1:]:
                    if self.uf[q] not in self._groups0:
                        # print(i, ": Merge", uf[q], "with", oldp, "via", p)
                        self._groups0[self.uf[q]] = (bl, bl - v, p)
                    self.uf.union(oldp, q)

        self._groups0 = [
            (k, self._groups0[k][0], self._groups0[k][1], self._groups0[k][2])
            for k in self._groups0
        ]
        self._groups0.sort(key=lambda g: g[2], reverse=True)
        self.persistence = self._groups0

    def get_comp_birth(self, p):
        return self.get(self.uf[p])

    def get(self, p):
        return self.image[p[0]][p[1]]

    def iter_neighbors(self, p, w, h):
        y, x = p

        # 8-neighborship
        neigh = [(y + j, x + i) for i in [-1, 0, 1] for j in [-1, 0, 1]]
        # 4-neighborship
        # neigh = [(y-1, x), (y+1, x), (y, x-1), (y, x+1)]

        for j, i in neigh:
            if j < 0 or j >= h:
                continue
            if i < 0 or i >= w:
                continue
            if j == y and i == x:
                continue
            yield j, i
