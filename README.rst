scikit-ued
==========

.. image:: https://img.shields.io/appveyor/ci/LaurentRDC/scikit-ued/master.svg
    :target: https://ci.appveyor.com/project/LaurentRDC/scikit-ued
    :alt: Windows Build Status
.. image:: https://readthedocs.org/projects/scikit-ued/badge/?version=master
    :target: http://scikit-ued.readthedocs.io
    :alt: Documentation Build Status
.. image:: https://img.shields.io/pypi/v/scikit-ued.svg
    :target: https://pypi.python.org/pypi/scikit-ued
    :alt: PyPI Version

Collection of algorithms and functions for ultrafast electron diffraction. It aims to be a fully-tested package
taking advantage of Python's most recent features.

Examples
--------

First and foremost, ``skued`` is concerned with the diffraction of crystals. Hence, building a ``Crystal`` object
from various sources is easy::

    from skued import Crystal

    tise2 = Crystal.from_cif('tise2.cif')
    graphite = Crystal.from_database('C')       # Internal database
    hemoglobin = Crystal.from_pdb('1gzx')       # Protein Data Bank
    vo2 = Crystal.from_cod(1521124)             # Crystallography Open Database

The ``Crystal`` object encodes multiple attributes and methods related to space-group, symmetry, and scattering.

Another important part of ultrafast electron diffraction is image processing. Image-alignment, 
streaming image operations, exploitation of symmetry and more are included.

Baseline-determination in polycrystalline data is a cornerstone of this package. Using the dual-tree complex 
wavelet transform, time-varying baselines can be extracted with high accuracy, even in cases where diffraction
peaks overlap over the entire dataset. 

For more examples, see our `tutorials <http://scikit-ued.readthedocs.io/en/latest/tutorials/index.html>`_.

API Reference
-------------

The `API Reference on readthedocs.io <http://scikit-ued.readthedocs.io>`_ provides API-level documentation, as 
well as tutorials.

Installation
------------

scikit-ued is available on PyPI; it can be installed with `pip <https://pip.pypa.io>`_.::

    python -m pip install scikit-ued

To install the latest development version from `Github <https://github.com/LaurentRDC/scikit-ued>`_::

    python -m pip install git+git://github.com/LaurentRDC/scikit-ued.git

Each version is tested against Python 3.5 and 3.6. If you are using a different version, tests can be run
using the standard library's `unittest` module.

After installing scikit-ued you can use it like any other Python module as ``skued``.

Citations
---------

If you are using the baseline-removal functionality of scikit-ued, please consider citing the following publication:

    .. [#] L. P. René de Cotret and B. J. Siwick, A general method for baseline-removal in ultrafast 
           electron powder diffraction data using the dual-tree complex wavelet transform, Struct. Dyn. 4 (2017) DOI: 10.1063/1.4972518.

Support / Report Issues
-----------------------

All support requests and issue reports should be
`filed on Github as an issue <https://github.com/LaurentRDC/scikit-ued/issues>`_.

License
-------

scikit-ued is made available under the MIT License. For more details, see `LICENSE.txt <https://github.com/LaurentRDC/scikit-ued/blob/master/LICENSE.txt>`_.
