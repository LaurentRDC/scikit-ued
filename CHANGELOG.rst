Changelog
=========

Release 2.1.13
--------------

* Fixed some deprecation warnings due to `matplotlib` 3.6

Release 2.1.12
--------------

* Improved the masking functionality in `bragg_peaks_persistence`.

Release 2.1.11
--------------

* Added new Bragg peak determination functionality to combat new datasets with sub-optimal signal-to-noise ratio. Further added Brillouin zone determination based on Bragg peak locations (#42). 

Release 2.1.10
--------------

* Fixed an issue where masks needed to be provided in `bragg_peaks` (#41).

Release 2.1.9
-------------

* Updated bounds on `pyqtgraph`.

Release 2.1.8
-------------

* Added explicit testing with python 3.10.
* Cleaned up unused imports.
* Updated some modules which were using deprecated code.

Release 2.1.7
-------------

* This release brings no changes, and only fixes the conda-forge package.

Release 2.1.6
-------------

* Fixed an issue where :func:`gaussian` would trip on a full-width at half-maximum of 0.
* Fixed an issue where the first stage of the dual-tree complex wavelet transform was not shifted properly (#36).

Release 2.1.5
-------------

* Releases are now automatically performed using Github Actions
* It is now possible to install all the dependencies required to use :func:`diffshow` using the following installation option: ``pip install scikit-ued[diffshow]``.

Release 2.1.4
-------------

* Increased the reliability of :func:`bragg_peaks` to distinguish between noise and Bragg peaks.

Release 2.1.3
-------------

* Added the function :func:`bragg_peaks` to determine the location of single-crystal diffraction peaks in an image.
* Fixed deprecation warnings regarding NumPy's dtypes.

Release 2.1.2
-------------

* Improved :func:`autocenter` for diffraction patterns with large Ewald sphere walkoff.
* :func:`diffread` now supports NumPy's ``*.npy`` format.
* Speedup of all routines that use the Fast Fourier transform (:func:`autocenter`, :func:`align`, :func:`ialign`, :func:`itrack_peak`, and :func:`kinematicsim`) by 50%.

Release 2.1.1
-------------

* Added the :func:`autocenter` routine, to automatically find the center of diffraction patterns. This works for both single-crystal and polycrystalline patterns.
* `Support for Python 3.6 and NumPy<1.17 has been dropped <https://numpy.org/neps/nep-0029-deprecation_policy.html>`_

Release 2.1.0
-------------

This release brings major infrastructure improvements, which in turn have allowed to squash some bugs.

* Migration of continuous integration testing to GitHub Actions.
* Migration of test infrastructure to pytest.
* Tests are now included in source distributions.
* Code snippets in documentation are now tested for correctness.

User-facing changes:

* Fixed an issue where the detected shift in `skued.align` and `skued.ialign` might be partial (i.e. only shift in one direction).
* Fixed an issue with packaging data on Linux.
* The `fast` argument to `skued.align` and `skued.ialign` has been removed. It was previously-marked as deprecated.
* Added pre-emptive support for Python 3.10+ by removing deprecations.
* Increased the precision of the pseudo-voigt approximation in `skued.pseudo_voigt`.
* Fixed many issues regarding documentation being out-of-date.

Release 2.0.6
-------------

* `scikit-ued` is being re-licensed from the MIT license to the GPLv3 license.
* The `fast` argument to `skued.align` and `skued.ialign` has been deprecated. Its value has no effect anymore.
* Official support for Python 3.9.
* Removed explicit requirement for the `tifffile` package.

Release 2.0.5
-------------

* Added `skued.kinematicsim`, a simple function to compute electron diffraction patterns from 
  crystals structures in the kinematic approximation (i.e. thin samples).
* Added the `skued.RingArcSelection` area.
* Various documentation improvements and fixes.

Release 2.0.4
-------------

* Added support for `crystals.ElectronicStructure`. This requires `crystals` version 1.1.0 and up.
* Added the function `with_irf`, which allows to modify fitting functions to include the effects of instrument response.
* Various documentation fixes.

Release 2.0.3
-------------

* Added the `Selection.mpatch` method to draw patches on Matplotlib plots.
* Added the `spectrum_cmap` Matplotlib colormap, available under the name `"spectrum"`.
* Fixed an issue where diffracted intensities were not correctly scaled in `potential_map`. 

Release 2.0.2
-------------

* Added the :class:`Selection` class and :class:`RectSelection`, :class:`DiskSelection`, :class:`RingSelection`, and 
  :class:`ArbitrarySelection` to assemble time-series. This is a generalization of iris-ued's time-series rects.
* Added real-time pixel value and cursor position to ``skued.diffshow``.
* Added `indices_to_text`, a plotting utility function to render Miller indices to Mathjax/LaTeX-style text (Matplotlib-compatible).

Release 2.0.1
-------------

* ``skued.diffshow`` will temporarily switch PyQtGraph's image axis order to the row-major, which is a saner default.
* Added skued command-line utilities. Images can be shown (with interactive viewer) using ``skued diffshow [path]``.
  Crystal information can be determined using ``skued crystinfo [path]``.
* Fixed an issue where a typo in ``electron_velocity`` would raise an exception.

Release 2.0.0
-------------

Due to a conflict between scikit-image and scikit-ued conventions, some breaking changes are required. 
Image conventions will now follow that of scikit-image. Most importantly:

* Changed the convention on image masks to align with the scikit-image convention. Masks will be ``True`` for valid pixels, and ``False`` on invalid pixels.

We took the opportunity to make other breaking changes:

* Broke off the ``skued.structure`` package into its own library, ``crystals``.
* Removed `masked_register_translation` in favour of the new scikit-image implementation ported from scikit-ued.
* Removed `xcorr` and `mnxc` as these were the backbone of `masked_register_translation` and are no longer needed.
* Added aspherical electron form factor parametrization from Zheng et al. 2009.
* Removed ``diff_register`` in favor of an analog of scikit-image's `register_translation` and `masked_register_translation`. 
* Removed `powder_center` due to unpredictable performance. 
* Removed `calibrate_scattvector`, which was deprecated.
* Removed `time_shift` and `time_shifts`, which were deprecated.
* Removed `shift_image` in favor of `scipy.ndimage.shift`.
* `bounded_reflections` has been removed in favor of ``Crystal.bounded_reflections`` in the crystals library (version >= 0.6.4)

We have also added some features:

* Added the `patterson` function to calculate Patterson pair-pair distribution functions from polycrystalline diffraction patterns.
* Added the `detector_scattvectors` function to determine the wavevectors visible on a detector, in transmission,
  based on experimental geometry.

Release 1.0.1.1
---------------

* Added time-series fitting.

Release 1.0.1.0
---------------

* Added support for Gatan Digital Micrograph image formats DM3 and DM4

Release 1.0.0.0
---------------

* ``available_dt_filters`` and ``available_first_stage_filters`` have been added to list available baseline-removal filters.
