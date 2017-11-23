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

==============
Image Analysis
==============

Combine the routines below with
`npstreams`_ to process diffraction data in parallel. Refer 
to the :ref:`tutorial on image manipulation <image_analysis_tutorial>` 
for some examples.

.. autosummary::
    :toctree: functions/
    :nosignatures:

    azimuthal_average
    powder_center
    align
    ialign
    diff_register
    shift_image
    nfold
    mnxc2
    mask_from_collection
    combine_masks
    mask_image
    snr_from_collection
    isnr
    triml
    trimr
    
=================
Crystal structure
=================
Handling crystal structure information is crucial for many data analysis and modelling tasks.
See the :ref:`Structure tutorial <structure_tutorial>` for some examples on how to use the following
classes.

.. autosummary::
    :toctree: classes/
    :nosignatures:
    
    Crystal
    Atom
    Lattice

================
Dataset Handling
================

.. autosummary::
    :toctree: classes/
    :nosignatures:

    AbstractRawDataset
    McGillRawDataset

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