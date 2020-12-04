# -*- coding: utf-8 -*-

from datetime import datetime

import numpy as np


def mibheader(filepath, hoffset=0):
    """
    Get an image header from a Merlin Image Binary file.

    Parameters
    ----------
    filepath : str
        Path to the Merlin Image Binary file.
    hoffset : int, optional
        Bytes-offset of the header. For use in multi-image files only.

    Returns
    -------
    header : dict
        Header information. The following keys are available:

        * ``'ID'`` : Detector identification
        * ``'seq_num'`` : Image sequence number
        * ``'offset'`` : Length of the header in bytes
        * ``'nchips'`` :  Number of readout circuits used to build the sensor.
        * ``'shape'`` : image shape tuple (width, height)
        * ``'dtype'`` : resulting data-type.
        * ``'timestamp'`` : epoch representation of the time-stamp, down to the micro-second

    Notes
    -----
    Merlin Image Binary files can be composed of multiple images; in this case, the
    file is composed of alternating image headers and binary data.
    """
    # First step : read small part of the header
    # to get the offset
    with open(filepath, "r+b") as img_file:
        pre_header = img_file.read(20)
        _, _, offset, *_ = pre_header.decode("ascii").split(",")
        img_file.seek(hoffset)
        header = img_file.read(int(offset))

    header_items = header.decode("ascii").split(",")

    # More header items could be implemented if need be
    # but this is good enough for now.
    (
        header_ID,
        seq_num,
        data_offset,
        nchips,
        size_x,
        size_y,
        dtype_str,
        _,
        _,
        timestamp,
        *_,
    ) = header_items

    # Parse the date
    timestamp = datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S.%f").timestamp()

    # Special case : dtype
    # For NumPy, U16 is unicode; unsigned 16-bit (2 bytes) integers is u2
    # Also, must be big-endian ('>')
    nbytes = int(dtype_str[1::]) // 8
    dtype = np.dtype(">u" + str(nbytes))

    return {
        "ID": header_ID,
        "seq_num": int(seq_num),
        "offset": int(data_offset),
        "nchips": int(nchips),
        "shape": (int(size_x), int(size_y)),
        "dtype": dtype,
        "timestamp": timestamp,
    }


def imibread(filepath):
    """
    Generator of images contained in a Merlin Image Binary file as NumPy arrays.

    Parameters
    ----------
    filepath : str
        Path to the Merlin Image Binary file.

    Yields
    ------
    im : `~numpy.ndarray`, ndim 2
        NumPy array of the Merlin Image Binary.

    See Also
    --------
    mibread : Extract images from a MIB file as a dense NumPy array
    """
    coffset = 0  # current image offset (header + data)

    with open(filepath, "r+b") as binary:

        # Move to end of first image and check if there is another one
        # by reading one more byte after the image data
        # in the case of a single-image file, the next byte is b''
        while binary.read(1):
            binary.seek(coffset)

            # Information for the current image's header
            # these should not change from image to image
            header = mibheader(filepath, hoffset=coffset)
            size_x, size_y = header["shape"]
            im_dtype = header["dtype"]

            binary.seek(coffset + header["offset"])
            arr = np.fromfile(binary, dtype=im_dtype, count=size_x * size_y)
            yield np.reshape(arr, newshape=(size_x, size_y))

            coffset += header["offset"] + arr.nbytes


def mibread(filepath):
    """
    Read a MIB (Merlin Image Binary) file into a NumPy array.

    Parameters
    ----------
    filepath : str
        Path to the Merlin Image Binary file.

    Returns
    -------
    im : `~numpy.ndarray`
        NumPy array of the Merlin Image Binary. In the case of multi-image files,
        images will be stacked along axis 2.

    See Also
    --------
    imibread : generate images contained in a MIB file.
    """
    # Squeeze the resulting array so that single image array have shape (x, y)
    # and not (x, y, 1)
    images = imibread(filepath)
    return np.squeeze(np.dstack(tuple(images)))
