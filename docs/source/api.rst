.. _api:

*************
Reference/API
*************

Parallel Utilities
==================
.. autofunction:: skued.pmap

.. autofunction:: skued.preduce

Plot Utilities
==============
.. autofunction:: skued.spectrum_colors

.. autofunction:: skued.rgb_sweep

Array Utilities
===============
.. automodule:: skued.array_utils

Iteration/Generator Utilities
=============================
.. automodule:: skued.iter_utils

Quantities
==========
.. automodule:: skued.quantities

Voigt Profile
=============
.. automodule:: skued.voigt

Affine Transforms
=================
.. automodule:: skued.affine

Crystal structure
=================
Handling crystal structure information is crucial for many data analysis and modelling tasks.
See the :ref:`Structure tutorial <structure_tutorial>` for some examples on how to use the following
classes.

.. autoclass:: skued.Crystal
   :members:

.. autoclass:: skued.Atom
   :members:

.. autoclass:: skued.Lattice
   :members:

Simulation
==========
.. autofunction:: skued.simulation.powdersim

Baseline-determination
======================
.. autofunction:: skued.baseline.baseline_dt

.. autofunction:: skued.baseline.baseline_dwt

.. autofunction:: skued.baseline.dtcwt

Image Analysis
==============
.. automodule:: skued.image.powder

.. automodule:: skued.image.alignment

.. automodule:: skued.image.streaming

.. automodule:: skued.image.symmetry

.. automodule:: skued.image.correlation