.. include:: references.txt

.. _installation:

************
Installation
************

Requirements
============

Scikit-ued works on Linux, Mac OS X and Windows. It requires Python 3.7+. Packages requirements are `listed here <https://github.com/LaurentRDC/scikit-ued/blob/master/requirements.txt>`_.

Install scikit-ued
==================

You can install the latest stable version from PyPI::

    pip install scikit-ued

To also install optional dependencies required to view diffraction images interactively with :func:`diffshow`::

    python -m pip install scikit-ued[diffshow]

scikit-ued is also available on the conda-forge channel for the `conda package manager <https://conda.io/docs/>`_::

    conda config --add channels conda-forge
    conda install scikit-ued

You can install the latest **developer** version of scikit-ued by cloning the git
repository::

    git clone https://github.com/LaurentRDC/scikit-ued.git

...then installing the package with::

    cd scikit-ued
    python -m pip install .

Alternatively, if you want to use :func:`diffshow`::

    git clone https://github.com/LaurentRDC/scikit-ued.git
    cd scikit-ued
    python -m pip install .[diffshow]


Testing
=======

Testing requires `pytest`. If you want to check that all the tests are running correctly with your Python
configuration, type::

    python -m pytest --pyargs skued