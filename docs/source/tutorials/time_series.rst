.. include:: ../references.txt

.. _timeseries_tutorial:

********************
Time-series analysis
********************

Time-resolved diffraction measurements reduce down to a (large) series of time-series. Exploring those time-series
reveal the dynamics captured by experiments. This tutorial will go over some of the functions that help this data exploration.

Contents
========

* :ref:`fitting`

.. _fitting:

Fitting
=======

Time-resolved experiments are powerful tools because they can be used to distinguish separate processes by their characteristic time-scales.
Therefore, extracting time-constants from dynamics is an important tool.

:mod:`skued` exports some curves that make it much easier to fit to scattering time-series::

    import numpy as np
    import skued
    from scipy.optimize import curve_fit

    # Load data from file first
    # 2 x N array:
    #   first row is time-delay
    #   second row is diffracted intensity
    block = np.load('tseries1.npy')
    time, intensity = block[0, :], block[1, :]

    # Compute initial guesses for this curve (optional)
    initial_guesses = (0,                                   # time-zero
                       intensity.max() - intensity.min(),   # amplitude
                       1,                                   # time-constant
                       intensity.min())                     # offset
    
    params, pcov = curve_fit(skued.exponential, time, intensity, 
                             p0 = initial_guesses)
                                                         
    tzero, amplitude, tconst, offset = params
    best_fit_curve = skued.exponential(time, *params)
    # Equivalent: 
    #   best_fit_curve = skued.exponential(time, tzero, amplitude, tconst, offset)

We can plot the result:

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    import skued
    from scipy.optimize import curve_fit

    # Load data from file first
    # 2 x N array, first row is time-delay, second row is diffracted intensity
    block = np.load('tseries1.npy')
    time, intensity = block[0, :], block[1, :]

    # Compute initial guesses for this curve
    initial_guesses = (0,                                   # time-zero
                       intensity.max() - intensity.min(),   # amplitude
                       1,                                   # time-constant
                       intensity.min())                     # offset
    
    params, pcov = curve_fit(skued.exponential, time, intensity, 
                                                         p0 = initial_guesses)
                                                         
    tzero, amplitude, tconst, offset = params
    best_fit_curve = skued.exponential(time, *params)
    # Equivalent: 
    #   best_fit_curve = skued.exponential(time, tzero, amplitude, tconst, offset)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.axvline(x = tzero, linestyle = 'dashed', color = 'grey', label = 'Time-zero')
    ax.plot(time, intensity, '.k', label = 'Time-series')
    ax.plot(time, best_fit_curve, '-r', label = 'Best fit')
    ax.set_xlabel('Time-delay [ps]')
    ax.set_ylabel('Diffracted Intensity [a.u.]')
    ax.set_xlim([time.min(), 40])
    ax.legend()
    plt.show()

:ref:`Return to Top <timeseries_tutorial>`