# -*- coding: utf-8 -*-


import numpy as np
from skimage.filters import gaussian

from crystals import Crystal
from skued import detector_scattvectors, electron_wavelength, powder_calq, powdersim

np.random.seed(23)


def test_powder_calq_simulation():
    """
    Test calibration from simulation, down to 1% error . Peaks (200) and (220) from monoclinic VO2 are used,
    """

    s = np.linspace(0.11, 0.8, 1024)
    q = 4 * np.pi * s
    c = Crystal.from_database("vo2-m1")
    I = powdersim(c, s)

    peak1 = (2, 0, 0)
    Gx1, Gy1, Gz1 = c.scattering_vector(peak1)
    q1 = np.sqrt(Gx1**2 + Gy1**2 + Gz1**2)
    arr_index1 = np.argmin(np.abs(q - q1))

    peak2 = (2, 2, 0)
    Gx2, Gy2, Gz2 = c.scattering_vector(peak2)
    q2 = np.sqrt(Gx2**2 + Gy2**2 + Gz2**2)
    arr_index2 = np.argmin(np.abs(q - q2))

    calibrated = powder_calq(
        I, c, peak_indices=(arr_index1, arr_index2), miller_indices=(peak1, peak2)
    )

    assert I.shape == calibrated.shape
    assert np.allclose(q, calibrated, rtol=0.01)


def test_powder_calq_simulation_3_peaks():
    """
    Test calibration from simulation, down to 1% error . Peaks (200) and (220) from monoclinic VO2 are used,
    """
    s = np.linspace(0.11, 0.8, 1024)
    q = 4 * np.pi * s
    c = Crystal.from_database("vo2-m1")
    I = powdersim(c, s)

    peak1 = (2, 0, 0)
    Gx1, Gy1, Gz1 = c.scattering_vector(peak1)
    q1 = np.sqrt(Gx1**2 + Gy1**2 + Gz1**2)
    arr_index1 = np.argmin(np.abs(q - q1))

    peak2 = (2, 2, 0)
    Gx2, Gy2, Gz2 = c.scattering_vector(peak2)
    q2 = np.sqrt(Gx2**2 + Gy2**2 + Gz2**2)
    arr_index2 = np.argmin(np.abs(q - q2))

    peak3 = (3, 0, -2)
    Gx2, Gy2, Gz2 = c.scattering_vector(peak3)
    q3 = np.sqrt(Gx2**2 + Gy2**2 + Gz2**2)
    arr_index3 = np.argmin(np.abs(q - q3))

    calibrated = powder_calq(
        I,
        c,
        peak_indices=(arr_index1, arr_index2, arr_index3),
        miller_indices=(peak1, peak2, peak3),
    )

    assert I.shape == calibrated.shape
    assert np.allclose(q, calibrated, rtol=0.01)


def test_detector_scattvectors_center():
    """Test that the placement of the center is working as intended."""
    qx, qy, qz = detector_scattvectors(
        keV=200,
        camera_length=1,
        shape=(512, 512),
        pixel_size=1e-6,
        center=(128, 128),
    )

    q_parallel = np.sqrt(qx**2 + qy**2)
    assert np.unravel_index(np.argmin(q_parallel), qx.shape) == (128, 128)


def test_detector_scattvectors_default_center():
    """Test that the placement of the center
    by default is in the center of the detector."""
    qx, qy, qz = detector_scattvectors(
        keV=200, camera_length=1, shape=(512, 512), pixel_size=1e-6, center=None
    )

    q_parallel = np.sqrt(qx**2 + qy**2)
    assert np.unravel_index(np.argmin(q_parallel), qx.shape) == (256, 256)


def test_detector_scattvectors_ewald_radius():
    """Test that the norm of total scattering vector norm is constant and equal to the
    Ewald sphere radius"""
    qx, qy, qz = detector_scattvectors(
        keV=200, camera_length=1, shape=(128, 128), pixel_size=1e-6, center=None
    )

    q_norm = np.sqrt(qx**2 + qy**2 + qz**2)
    ewald_sphere_radius = 2 * np.pi / electron_wavelength(keV=200)

    assert np.allclose(q_norm, ewald_sphere_radius)
