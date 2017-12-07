.. include:: references.txt

.. _installation:

************
Installation
************

Requirements
============

.. note::

    Users are strongly recommended to manage these dependencies with the
    excellent `Intel Distribution for Python <https://software.intel.com/en-us/intel-distribution-for-python>`_
    which provides easy access to all of the above dependencies and more.

works on Linux, Mac OS X and Windows. It requires Python 3.5+ 
as well as the following packages:

* `numpy`_
* `scipy`_
* `scikit-image`_
* `pywavelets`_


Install scikit-ued
==================

You can install the latest developer version of scikit-ued by cloning the git
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
