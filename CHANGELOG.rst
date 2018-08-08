Changelog
=========

Release 1.0.2.0 (Development)
---------------------------

* Added control for mask overlap in `mnxc2`. `mnxc2` now passes the same tests as reference implementation.
* Added an analog of scikit-image's `register_translation`, `masked_register_translation`. It will eventually replace `diff_register`.

* Deprecated `powder_center` due to unpredictable performance. Warnings will be issued on every function call. 
  It will be removed in an upcoming release

* Removed `calibrate_scattvector`, which was deprecated.
* Removed `time_shift` and `time_shifts`, which were deprecated.

Release 1.0.1.0
---------------

* Added support for Gatan Digital Micrograph image formats DM3 and DM4
