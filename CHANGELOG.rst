Changelog
=========

Release 2.0.0 (Development)
---------------------------

* Added aspherical electron form factor parametrization from Zheng et al. 2009.
* Added an analog of scikit-image's `register_translation`, `masked_register_translation`. It will eventually replace `diff_register`.

* Removed `mnxc2` in favor of `mnxc`, the n-dimensional masked normalized cross-correlation.
* Removed `powder_center` due to unpredictable performance. Warnings will be issued on every function call. It will be removed in an upcoming release

* Removed `calibrate_scattvector`, which was deprecated.
* Removed `time_shift` and `time_shifts`, which were deprecated.

Release 1.0.1.0
---------------

* Added support for Gatan Digital Micrograph image formats DM3 and DM4
