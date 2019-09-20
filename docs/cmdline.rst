.. include:: references.txt

.. _cmdline:

**********************
Command-line utilities
**********************

.. currentmodule:: skued

Scikit-ued includes command-line utilities for repetitive tasks. To see the available 
commands on your system::

    skued --help

Interactive image viewer
------------------------

One command-line utility is the interactive image viewer. To read and display
a diffraction pattern/image located at `path/to/file`::

    skued diffshow path/to/file

Make sure that `PyQtGraph`, an optional dependency, is installed.