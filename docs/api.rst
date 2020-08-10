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
    available_dt_filters
    available_first_stage_filters

======================
Structure manipulation
======================

.. currentmodule:: crystals

Structure manipulation is provided by the `crystals <https://crystals.rtfd.io>`_ library. A summary of the documentation
is provided below.

.. autosummary::
    :toctree: classes/
    :nosignatures:

    Crystal
    Atom

Bases classes
-------------

The :class:`Lattice` class allows for manipulating lattice information separately from
atomic information.

.. autosummary::
    :toctree: classes/
    :nosignatures:

    Lattice
    LatticeSystem
    AtomicStructure

Parsers
-------

Structure parsers are used to build :class:`Crystal` instances, mostly through :class:`Crystal` class methods.

.. autosummary::
    :toctree: classes/
    :nosignatures:

    CIFParser
    CODParser
    PDBParser

.. currentmodule:: skued

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

    register_time_shift
    register_time_shifts

Robust statistics
-----------------

.. autosummary::
    :toctree: functions/
    :nosignatures:

    mad

Fitting
-------

.. autosummary::
    :toctree: functions/
    :nosignatures:

    exponential
    biexponential
    with_irf

Non-uniform Fast Fourier Transform
----------------------------------

.. autosummary::
    :toctree: functions/
    :nosignatures:

    nfft
    nfftfreq

Time-series Selections
----------------------

.. autosummary::
    :toctree: classes/
    :nosignatures:

    Selection
    RectSelection
    DiskSelection
    RingSelection
    RincArcSelection
    ArbitrarySelection

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
    reflection

Polycrystalline diffraction
---------------------------

.. autosummary::
    :toctree: functions/
    :nosignatures:

    azimuthal_average

Calibrations
------------

.. autosummary::
    :toctree: functions/
    :nosignatures:
    
    powder_calq

Image alignment
---------------

.. autosummary::
    :toctree: functions/
    :nosignatures:

    align
    ialign
    itrack_peak

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
    kinematicsim
    electrostatic
    pelectrostatic

=============
Visualization
=============

.. autosummary::
    :toctree: functions/
    :nosignatures:

    potential_map
    potential_synthesis
    patterson

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

.. note::

    You should prefer to use :func:`diffread` to the specific format 
    functions like :func:`mibread` and :func:`dmread`

Digital Micrograph DM3/DM4:

.. autosummary::
    :toctree: functions/
    :nosignatures:

    dmread

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
    spectrum_cmap
    rgb_sweep
    indices_to_text

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

============================
Thin Film Optical Properties
============================

.. autosummary::
    :toctree: functions/
    :nosignatures:

    film_optical_coefficients

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