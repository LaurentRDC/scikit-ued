Changelog
=========

Release 2.0.6
-------------

* `scikit-ued` is being re-licensed from the MIT license to the GPLv3 license.
* The `fast` argument to `skued.align` and `skued.ialign` has been deprecated. Its value has no effect anymore.
* Official support for Python 3.9.

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
