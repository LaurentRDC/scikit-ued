# -*- coding: utf-8 -*-
import numpy as np

from npstreams import deprecated

def hypot(*args):
    """ Generalized np.hypot """
    return np.sqrt(np.sum(np.square(args)))

# TODO: test
#       Not exporting this function until tests are written
def calq(I, crystal, peak_indices, miller_indices):
    """
    Determine the scattering vector q corresponding to a diffraction pattern
    and a known crystal structure.

    Parameters
    ----------
    I : `~numpy.ndarray`, ndim 2
        Diffraction pattern. It is assumed that the diffraction
        pattern is defined on an equidistant grid.
    crystal : skued.Crystal instance
        Crystal that gave rise to the diffraction pattern ``I``.
    peak_indices : 2-tuple of ints
        Array index location of two diffraction peaks in the array ``I``. For best
        results, peaks should be well-separated.
    miller_indices : iterable of 3-tuples
        Indices associated with the peaks of ``peak_indices``.
        E.g. ``indices = [(2,2,0), (-3,0,2)]``

    Returns
    -------
    qx, qy : `~numpy.ndarray`, ndim 2
        Scattering vectors (x and y) associated with the intensity profile ``I``.
    
    Raises
    ------
    ValueError : if ``I`` is not a 2D diffraction pattern.
    """
    I = np.asarray(I)

    if I.ndim != 2:
        raise ValueError('Expected 2D diffraction pattern, but received shape {}'.format(I.shape))

    hkl1, hkl2 = miller_indices
    qx1, qy1, _ = crystal.scattering_vector(*hkl1)
    qx2, qy2, _ = crystal.scattering_vector(*hkl2)

    # calibration is done by fitting a line
    # Expecting that I is defined on an 
    # equally-spaced grid [0, 1, ..., I.size]
    peak_indices_x = [peak_index[0] for peak_index in peak_indices]
    peak_indices_y = [peak_index[1] for peak_index in peak_indices]

    slope_x, intercept_x = np.polyfit(np.asarray(peak_indices_x), np.asarray([qx1, qx2]), deg = 1)
    slope_y, intercept_y = np.polyfit(np.asarray(peak_indices_y), np.asarray([qy1, qy2]), deg = 1)

    x_range = slope_x * np.arange(0, I.shape[0]) + intercept_x
    y_range = slope_y * np.arange(0, I.shape[1]) + intercept_y
    return np.meshgrid(x_range, y_range, indexing = 'xy')

@deprecated('calibrate_scattvector has been renamed to powder_calq')
def calibrate_scattvector(*args, **kwargs):
    return powder_calq(*args, **kwargs)

def powder_calq(I, crystal, peak_indices, miller_indices):
    """
    Determine the scattering vector q corresponding to a polycrystalline diffraction pattern
    and a known crystal structure.

    For best results, multiple peaks (and corresponding Miller indices) should be provided; the
    absolute minimum is two.

    Parameters
    ----------
    I : `~numpy.ndarray`, ndim 1
        Polycristalline diffraction pattern. It is assumed that the diffraction
        pattern is defined on an equidistant grid.
    crystal : skued.Crystal instance
        Crystal that gave rise to the diffraction pattern ``I``.
    peak_indices : iterable of 2-tuple
        Array index location of two diffraction peaks in the array ``I``. For best
        results, peaks should be well-separated.
        E.g. ``peak_indices = [(1028, 123), (10, 891)]``
    miller_indices : iterable of 3-tuples
        Indices associated with the peaks of ``peak_indices``.
        E.g. ``indices = [(2,2,0), (-3,0,2)]``

    Returns
    -------
    q : `~numpy.ndarray`, ndim 1
        Scattering vectors associated with the intensity profile ``I``.
    
    Raises
    ------
    ValueError : if ``I`` is not a 1D diffraction pattern.
    ValueError : if the number of peak indices does not match the number of Miller indices.
    ValueError : if the number of peaks given is lower than two.
    """
    # I is not strictly required at this time. Only the shape of I is important
    # However, we might refine the peak positions in the future
    I = np.asarray(I)

    if I.ndim > 1:
        raise ValueError('Expected 1D diffraction intensity, but received shape {}'.format(I.shape))
    
    if len(peak_indices) != len(miller_indices):
        raise ValueError('Number of array indices {} does not match the \
                          number of Miller indices {}'.format(len(peak_indices), len(miller_indices)))
    
    if len(peak_indices) < 2:
        raise ValueError('Two peaks are required to calibrate, but received {}'.format(len(peak_indices)))
    
    # scattering vector length based on known structure
    qs = [hypot(*crystal.scattering_vector(*hkl)) for hkl in miller_indices]

    # calibration is done by fitting a line
    # Expecting that I is defined on an 
    # equally-spaced grid [0, 1, ..., I.size]
    slope, intercept = np.polyfit(np.asarray(peak_indices), np.asarray(qs), deg = 1)
    return slope * np.arange(0, I.size) + intercept