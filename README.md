scikit-ued
==========

[![Documentation Build Status](https://readthedocs.org/projects/scikit-ued/badge/?version=master)](http://scikit-ued.readthedocs.io) [![PyPI Version](https://img.shields.io/pypi/v/scikit-ued.svg)](https://pypi.org/project/scikit-ued/) [![Conda-forge Version](https://img.shields.io/conda/vn/conda-forge/scikit-ued.svg)](https://anaconda.org/conda-forge/scikit-ued) [![DOI badge](https://img.shields.io/badge/DOI-10.1186%2Fs40679--018--0060--y-blue)](https://doi.org/10.1186/s40679-018-0060-y)

Collection of algorithms and functions for ultrafast electron diffraction. It aims to be a fully-tested package taking advantage of Python's most recent features.

For examples, see our [tutorials](https://scikit-ued.readthedocs.io/).

API Reference
-------------

The [API Reference on readthedocs.io](https://scikit-ued.readthedocs.io) provides API-level documentation, as well as tutorials.

Installation
------------

scikit-ued is available on PyPI; it can be installed with [pip](https://pip.pypa.io):

    python -m pip install scikit-ued

To also install optional dependencies required to view diffraction images interactively:

    python -m pip install scikit-ued[diffshow]

scikit-ued is also available on the conda-forge channel for the [conda package manager](https://conda.io/docs/):

    conda config --add channels conda-forge
    conda install scikit-ued

To install the latest development version from [Github](https://github.com/LaurentRDC/scikit-ued):

    python -m pip install git+https://github.com/LaurentRDC/scikit-ued.git

After installing scikit-ued you can use it like any other Python module
as `skued`.

Each version is tested against **Python 3.7+**. If you are using a
different version, tests can be run using the `pytest` package.

Optional dependencies
---------------------

For displaying diffraction images with interactive contrast using the
`skued.diffshow` function, PyQtGraph is required.

Contributing
------------

If you want to contribute to `scikit-ued`, take a look at [`CONTRIBUTING.md`](https://github.com/LaurentRDC/scikit-ued/blob/master/CONTRIBUTING.md).

Related projects
----------------

Streaming operations on NumPy arrays are available in the [npstreams package](https://pypi.org/pypi/npstreams).

Interactive exploration of ultrafast electron diffraction data with the [iris-ued package](https://pypi.org/project/iris-ued/).

Crystal structure manipulation (including symmetry-determination) with the [crystals package](https://pypi.org/project/crystals/). (Included
with scikit-ued)

A graphical user interface for the dual-tree complex wavelet transform
baseline-removal routine is available as a [separate package](https://pypi.org/pypi/dtgui).

Citations
---------

If you find this software useful, please consider citing the following
publication:

> L. P. René de Cotret, M. R. Otto, M. J. Stern. and B. J. Siwick, *An open-source software ecosystem for the interactive exploration of ultrafast electron scattering data*, Advanced Structural and Chemical Imaging 4:11 (2018) [DOI: 10.1186/s40679-018-0060-y.](https://ascimaging.springeropen.com/articles/10.1186/s40679-018-0060-y)

If you are using the baseline-removal functionality of scikit-ued,
please consider citing the following publication:

> L. P. René de Cotret and B. J. Siwick, *A general method for baseline-removal in ultrafast electron powder diffraction data using the dual-tree complex wavelet transform*, Struct. Dyn. 4 (2017) [DOI: 10.1063/1.4972518](https://doi.org/10.1063/1.4972518).


Support / Report Issues
-----------------------

All support requests and issue reports should be [filed on Github as an issue](https://github.com/LaurentRDC/scikit-ued/issues).

License
-------

scikit-ued is made available under the GPLv3 License. For more details,
see [LICENSE.txt](https://github.com/LaurentRDC/scikit-ued/blob/master/LICENSE.txt).
