scikit-ued
==========

.. image:: https://img.shields.io/appveyor/ci/LaurentRDC/scikit-ued/master.svg
    :target: https://ci.appveyor.com/project/LaurentRDC/scikit-ued
    :alt: Windows Build Status
.. image:: https://readthedocs.org/projects/scikit-ued/badge/?version=latest
    :target: http://scikit-ued.readthedocs.io
    :alt: Documentation Build Status
.. image:: https://img.shields.io/pypi/v/scikit-ued.svg
    :target: https://pypi.python.org/pypi/scikit-ued
    :alt: PyPI Version

Collection of algorithms and functions for ultrafast electron diffraction.

Getting Started with scikit-ued
-------------------------------

scikit-ued is available on PyPI can be installed with `pip <https://pip.pypa.io>`_.::

    $ python -m pip install scikit-ued

To install the latest development version from `Github <https://github.com/LaurentRDC/scikit-ued>`_::

    $ python -m pip install git+git://github.com/LaurentRDC/scikit-ued.git

Each version is tested against Python 3.5 and 3.6. If you are using a different version, tests can be run
using the standard library's `unittest` module.

After installing scikit-ued you can use it like any other Python module as :code:`skued`.

Citations
---------

If you are using the :code:`skued.baseline` subpackage, consider citing the following publication:

    .. [#] L. P. René de Cotret and B. J. Siwick, A general method for baseline-removal in ultrafast 
           electron powder diffraction data using the dual-tree complex wavelet transform, Struct. Dyn. 4 (2017) DOI: 10.1063/1.4972518.

API Reference
-------------

The `API Reference on readthedocs.io <http://scikit-ued.readthedocs.io>`_ provides API-level documentation.

Support / Report Issues
-----------------------

All support requests and issue reports should be
`filed on Github as an issue <https://github.com/LaurentRDC/scikit-ued/issues>`_.

License
-------

scikit-ued is made available under the MIT License. For more details, see `LICENSE.txt <https://github.com/LaurentRDC/scikit-ued/blob/master/LICENSE.txt>`_.
