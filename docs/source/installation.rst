.. include:: references.txt

.. _installation:

************
Installation
************

Requirements
============

Scikit-ued works on Linux, Mac OS X and Windows. It requires Python 3.6+. Packages requirements are `listed here <https://github.com/LaurentRDC/scikit-ued/blob/master/requirements.txt>`_.

Install scikit-ued
==================

You can install the latest stable version from PyPI::

    pip install scikit-ued

You can install the latest **developer** version of scikit-ued by cloning the git
repository::

    git clone https://github.com/LaurentRDC/scikit-ued.git

...then installing the package with::

    cd scikit-ued
    python setup.py install

Installation on Windows
-----------------------

Some of scikit-ued's dependencies require compilation. If you are experiencing problems installing scikit-ued on Windows, here are some potential solutions:

    * Install a C/C++ compiler. The easiest way to do so is to install the `Visual Studio Build Tools <https://www.visualstudio.com/downloads/?q=build+tools>`_. More information is available on the `Python Wiki <https://wiki.python.org/moin/WindowsCompilers>`_.
    * Download the wheels from scikit-ued's `wheelhouse <https://github.com/LaurentRDC/scikit-ued/tree/master/wheelhouse>`_. These are pre-compiled dependencies that will only work on Windows. To install a wheel, you can use pip: `pip install some-pkg.whl`.
    * Install the dependencies using the `conda package manager <https://conda.io/docs/>`_. Most notably, spglib and pycifrw are both available in the conda-forge channel.

To install packages available in the conda-forge channel::

    conda config --add channels conda-forge
    conda install spglib pycifrw (...)


Testing
=======

If you want to check that all the tests are running correctly with your Python
configuration, type::

    python setup.py test

This will only work if you have downloaded a copy of the source code.