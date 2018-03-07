.. include:: references.txt

.. _api:

*************
Reference/API
*************

.. currentmodule:: skued

======================
Baseline-determination
======================

Please refer to the :ref:`tutorial on baseline-determination <baseline_tutorial>` for some examples.

.. autosummary::
    :toctree: functions/
    :nosignatures:

    baseline_dt
    baseline_dwt
    dtcwt
    idtcwt

====================
Time-series Analysis
====================

Time-series exploration and analysis.

Time-zero tracking
------------------

Measurement of time-shifts between physically-equivalent time traces:

.. autosummary::
    :toctree: functions/
    :nosignatures:

    time_shift
    time_shifts

Robust statistics
-----------------

.. autosummary::
    :toctree: functions/
    :nosignatures:

    mad

Non-uniform Fast Fourier Transform
----------------------------------

.. autosummary::
    :toctree: functions/
    :nosignatures:

    nfft
    nfftfreq

==============
Image Analysis
==============

Combine the routines below with
`npstreams`_ to process diffraction data in parallel. Refer 
to the :ref:`tutorial on image manipulation <image_analysis_tutorial>` 
for some examples.

Symmetry
--------

.. autosummary::
    :toctree: functions/
    :nosignatures:

    nfold

Polycrystalline diffraction
---------------------------

.. autosummary::
    :toctree: functions/
    :nosignatures:

    azimuthal_average
    powder_center
    calibrate_scattvector

Image alignment
---------------

.. autosummary::
    :toctree: functions/
    :nosignatures:

    align
    ialign
    diff_register
    shift_image
    itrack_peak

Correlations
------------

.. autosummary::
    :toctree: functions/
    :nosignatures:
    
    xcorr
    mnxc2

Image masking
-------------

.. autosummary::
    :toctree: functions/
    :nosignatures:

    mask_from_collection
    combine_masks
    mask_image

Image noise
-----------

.. autosummary::
    :toctree: functions/
    :nosignatures:

    snr_from_collection
    isnr
    triml
    trimr

==========
Simulation
==========

.. autosummary::
    :toctree: functions/
    :nosignatures:

    structure_factor
    affe
    powdersim
    electrostatic
    pelectrostatic
    bounded_reflections

============
Input/Output
============

General diffraction image I/O and plotting. Note that
for :func:`diffshow`, the packages PyQtGraph and PyQt5 must be
installed.

.. autosummary::
    :toctree: functions/
    :nosignatures:

    diffread
    diffshow

Merlin Image Binary (.mib) files:

.. autosummary::
    :toctree: functions/
    :nosignatures:

    mibheader
    mibread
    imibread

==============
Plot Utilities
==============

.. autosummary::
    :toctree: functions/
    :nosignatures:

    spectrum_colors
    rgb_sweep

===============
Array Utilities
===============

.. autosummary::
    :toctree: functions/
    :nosignatures:

    mirror
    repeated_array
    cart2polar
    polar2cart
    cart2spherical
    spherical2cart
    plane_mesh

===================
Electron Properties
===================

.. autosummary::
    :toctree: functions/
    :nosignatures:

    electron_wavelength
    electron_velocity
    interaction_parameter
    lorentz

=============
Voigt Profile
=============

.. autosummary::
    :toctree: functions/
    :nosignatures:

    gaussian
    lorentzian
    pseudo_voigt

=================
Affine Transforms
=================

.. autosummary::
    :toctree: functions/
    :nosignatures:

    affine_map
    transform
    change_of_basis
    change_basis_mesh
    is_basis
    is_rotation_matrix
    minimum_image_distance
    rotation_matrix
    translation_matrix
    translation_rotation_matrix