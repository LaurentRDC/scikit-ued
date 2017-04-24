.. include:: ../references.txt

.. _baseline_tutorial:

************************
Baseline Tutorials
************************

Due to the high electron cross-section, background signals (or baseline) are
much more of a problem for electron diffraction than equivalent X-ray experiments.

Contents
========

* :ref:`dual_tree_baseline`

.. _dual_tree_baseline:

Iterative Baseline Determination using the Dual-Tree Complex Wavelet Transform
==============================================================================

In the case of polycrystalline samples (i.e. 1D diffraction signals),
the definite way of removing background signals is to use an iterative
approach based on the dual-tree complex wavelet transform.

First, we load an example dataset::
    
    import numpy as np
    import matplotlib.pyplot as plt

    # TODO: this

:ref:`Return to Top <baseline_tutorial>`