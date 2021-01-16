.. include:: ../references.txt

.. _timeseries_tutorial:

.. currentmodule:: skued

********************
Time-series analysis
********************

Time-resolved diffraction measurements reduce down to a (large) series of time-series. Exploring those time-series
reveal the dynamics captured by experiments. This tutorial will go over some of the functions that help this data exploration.

Contents
========

* :ref:`fitting`
* :ref:`selection`

.. _fitting:

Fitting
=======

Time-resolved experiments are powerful tools because they can be used to distinguish separate processes by their characteristic time-scales.
Therefore, extracting time-constants from dynamics is an important tool.

:mod:`skued` exports some curves that make it much easier to fit to scattering time-series:

    >>> import numpy as np
    >>> import skued
    >>> from scipy.optimize import curve_fit
    >>> 
    >>> # Load data from file first
    >>> # 2 x N array:
    >>> #   first row is time-delay
    >>> #   second row is diffracted intensity
    >>> block = np.load('docs/tutorials/data/tseries1.npy')
    >>> time, intensity = block[0, :], block[1, :]
    >>> 
    >>> # Compute initial guesses for this curve (optional)
    >>> initial_guesses = (0,                                   # time-zero
    ...                    intensity.max() - intensity.min(),   # amplitude
    ...                    1,                                   # time-constant
    ...                    intensity.min())                     # offset
    >>>
    >>> params, pcov = curve_fit(skued.exponential, time, intensity, 
    ...                          p0 = initial_guesses)
    >>>                                                  
    >>> tzero, amplitude, tconst, offset = params
    >>> best_fit_curve = skued.exponential(time, *params)
    >>> # Equivalent: 
    >>> #   best_fit_curve = skued.exponential(time, tzero, amplitude, tconst, offset)

We can plot the result:

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    import skued
    from scipy.optimize import curve_fit

    # Load data from file first
    # 2 x N array, first row is time-delay, second row is diffracted intensity
    block = np.load('data/tseries1.npy')
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

Taking into account the IRF
---------------------------

No fitting of time-resolved data is complete without taking into account the instrument response function, 
or IRF. To help with this, `scikit-ued` provides the :func:`with_irf` decorator factory. It can be used as follows:

    >>> from skued import exponential, with_irf
    >>> 
    >>> # If the data is recorded in picoseconds, then the following
    >>> # applies a Gaussian IRF with a full-width at half-max of 150fs (0.15 picoseconds) 
    >>> @with_irf(0.150)
    ... def exponential_with_irf(time, *args, **kwargs):
    ...     return exponential(time, *args, **kwargs)
    
Let's see what :func:`with_irf` in action:

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    import skued

    times = np.linspace(-5, 15, 256)
    data = skued.exponential(times, 0, 1, 3)

    @skued.with_irf(4) # IRF is larger than expected time-constant
    def exponential_with_irf(time, *args, **kwargs):
        return skued.exponential(time, *args, **kwargs)
    conv = exponential_with_irf(times, 0, 1, 3)

    plt.figure()
    plt.plot(times, data, ".k", label='No IRF')
    plt.plot(times, conv, "-r", label='With IRF')
    plt.legend()
    plt.show()


:func:`with_irf` also works with data that is not defined on an even grid, like many ultrafast experiments:

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    import skued

    times = np.concatenate(
        (np.arange(-10, -2, step=1), np.arange(-2, 2, step=0.05), np.arange(2, 10, step=1))
    )
    data = skued.exponential(times, 0, 1, 1)
    conv = skued.with_irf(2)(skued.exponential)(times, 0, 1, 1)

    plt.figure()
    plt.plot(times, data, ".k", label='No IRF')
    plt.plot(times, conv, "-r", label='With IRF')
    plt.legend()
    plt.show()


.. _selection:

Selections
==========

In the context of ultrafast electron/x-ray scattering, time-series are assembled by selection a portion 
of scattering patterns for each time-delay. The :class:`Selection` class (and related subclasses) is 
the generalization of selecting a rectangular area of scattering patterns to arbitrary patterns, 
e.g. disks, torii, etc.

Instances can be used like boolean masks to select portions of scattering patterns. Consider an example
where we want to know the integrated intensity in a Bragg peak:

    >>> from skued import DiskSelection, diffread
    >>>
    >>> im = diffread(...) # doctest: +SKIP
    >>> 
    >>> bragg = DiskSelection(shape = im.shape, center=(1024, 1024), radius=30) # doctest: +SKIP
    >>> intensity = np.sum(im[bragg]) # doctest: +SKIP

Selections really shine when trying to extract non-standard shapes, e.g. rings. You can use them
in combination with `iris-ued` `DiffractionDatasets` to assemble specialized time-series.


:ref:`Return to Top <timeseries_tutorial>`