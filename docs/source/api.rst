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

Combine the routines in the :mod:`skued.image` module with
`npstreams`_ to process diffraction data in parallel. Please refer 
to the :ref:`tutorial on image manipulation <image_analysis_tutorial>` 
for some examples.

.. autosummary::
    :toctree: functions/

    image.azimuthal_average
    image.powder_center
    image.align
    image.ialign
    image.diff_register
    image.shift_image
    image.nfold
    image.mnxc2

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

.. autosummary::
    :toctree: functions/

    skued.structure.symmetry_expansion

Simulation
==========

.. autosummary::
    :toctree: functions/

    skued.simulation.powdersim

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

Structure Parsing
=================

.. autosummary::
    :toctree: classes/

    structure.CIFParser
    structure.PDBParser

Quantities
==========

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