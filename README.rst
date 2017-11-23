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

For examples, see our `tutorials <http://scikit-ued.readthedocs.io/en/latest/tutorials/index.html>`_.

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

Optional dependencies
---------------------

While it is not strictly required, the fourier transform routines from ``pyfftw`` will be preferred If
``pyfftw`` is installed.

Related projects
----------------

A graphical user interface for the dual-tree complex wavelet transform baseline-removal routine is available as a 
`separate package <https://pypi.python.org/pypi/dtgui>`_.

Streaming operations on NumPy arrays are available in the `npstreams package <https://pypi.python.org/pypi/npstreams>`_.

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
