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

Testing
=======

If you want to check that all the tests are running correctly with your Python
configuration, type::

    python setup.py test

This will only work if you have downloaded a copy of the source code.