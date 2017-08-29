.. include:: references.txt

.. _api:

*************
Reference/API
*************

.. currentmodule:: skued

Baseline-determination
======================

Please refer to the :ref:`tutorial on baseline-determination <baseline_tutorial>` for some examples.

.. autosummary::
    :toctree: functions/

    baseline_dt
    baseline_dwt

Dual-tree Complex Wavelet Transform
-----------------------------------
.. autosummary::
    :toctree: functions/

    dtcwt
    idtcwt

Image Analysis
==============

Combine the routines below with
`npstreams`_ to process diffraction data in parallel. Refer 
to the :ref:`tutorial on image manipulation <image_analysis_tutorial>` 
for some examples.

.. autosummary::
    :toctree: functions/

Routines for polycrystalline data
---------------------------------

.. autosummary::
    :toctree: functions/

    azimuthal_average
    powder_center

Image alignment
---------------

.. autosummary::
    :toctree: functions/

    align
    ialign
    diff_register
    shift_image
    nfold
    mnxc2

Pixel mask manipulations
------------------------

.. autosummary::
    :toctree: functions/

    mask_from_collection
    combine_masks
    mask_image

Image metrics
-------------

.. autosummary::
    :toctree: functions/

    snr_from_collection
    isnr
    
Crystal structure
=================
Handling crystal structure information is crucial for many data analysis and modelling tasks.
See the :ref:`Structure tutorial <structure_tutorial>` for some examples on how to use the following
classes.

.. autosummary::
    :toctree: classes/
    
    Crystal
    Atom
    Lattice

Simulation
==========

Polycrystalline diffraction
---------------------------

.. autosummary::
    :toctree: functions/

    powdersim

Electrostatic potential
-----------------------

.. autosummary::
    :toctree: functions/

    electrostatic
    pelectrostatic

Structure factor calculations
-----------------------------

.. autosummary::
    :toctree: functions/

    structure_factor
    structure_factor_miller

Miscellaneous
-------------

.. autosummary::
    :toctree: functions/

    bounded_reflections

Plot Utilities
==============

.. autosummary::
    :toctree: functions/

    spectrum_colors
    rgb_sweep

Array Utilities
===============

.. autosummary::
    :toctree: functions/

    mirror
    repeated_array
    cart2polar
    polar2cart
    plane_mesh

Electron Properties
===================

.. autosummary::
    :toctree: functions/

    electron_wavelength
    interaction_parameter
    lorentz

Voigt Profile
=============

.. autosummary::
    :toctree: functions/

    gaussian
    lorentzian
    pseudo_voigt

Affine Transforms
=================

.. autosummary::
    :toctree: functions/

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