.. include:: ../references.txt

.. _plotting_tutorial:

*****************
Plotting Tutorial
*****************

Time-resolved diffraction measurements can benefit from careful plotting, especially
for time-series of data. This tutorial will go over some of the functions that make this easier
in :mod:`skued`.

Contents
========

* :ref:`colors`

.. _colors:

Colors
======

Time-order
----------

Time-resolved data possesses obvious ordering, and as such we can use colors to
enhance visualization of such data.

My favourite way to plot time-series data is to make use of rainbow colors to indicate time:

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    from skued import spectrum_colors

    t = np.linspace(0, 10, num = 1000)
    y = np.exp(-t/0.4) + 10 *(1 - np.exp(-t/2)) + np.sin(2*np.pi*t)

    fig, ax = plt.subplots(1,1)
    ax.scatter(t, y, c = list(spectrum_colors(y.size)))

    ax.set_title('Some realistic dynamics')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Intensity (a.u.)')

    plt.show()

This functionality can be accessed through the :func:`spectrum_colors` generator::

    from skued import spectrum_colors

    colors = spectrum_colors(5)
    
    # Colors always go from purple to red (increasing wavelength)
    purple = next(colors)
    red = list(colors)[-1]

You can see some examples of uses of :func:`spectrum_colors` in the :ref:`baseline tutorial <baseline_tutorial>`.

:ref:`Return to Top <plotting_tutorial>`