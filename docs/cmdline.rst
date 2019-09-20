.. include:: references.txt

.. _cmdline:

**********************
Command-line utilities
**********************

.. currentmodule:: skued

Scikit-ued includes command-line utilities for repetitive tasks. To see the available 
commands on your system::

    python -m skued -h

or equivalently::

    python -m skued --help

Interactive image viewer
------------------------

One command-line utility is the interactive image viewer. To read and display
a diffraction pattern/image located at `path/to/file`::

    python -m skued diffshow path/to/file

Make sure that `PyQtGraph`, an optional dependency, is installed.