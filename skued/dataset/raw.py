# -*- coding: utf-8 -*-
"""
Handling UED Datasets
=====================
"""
from abc import ABC, abstractmethod, abstractproperty
from cached_property import cached_property

from os.path import isdir, join
from skimage.io import imread

class ExperimentalParameter:
    """
    Descriptor to experimental parameters and metadata.

    Parameters
    ----------
    name : str
        Parameter name
    ptype : type or callable
        Parameter type, e.g. float, or callable, e.g. numpy.array. 
    default : object or None
        Default value of the parameter. If None, no default value is set. Hence, the
        default value can never be None.
    readonly : bool, optional
        If True, attribute can never be set. Default is False.
    """
    def __init__(self, name, ptype, default = None, readonly = False):
        self.name = name
        self.type = ptype
        self.default = default
        self._readonly = readonly
    
    def __get__(self, instance, cls):
        if instance is None:
            return self
        return getattr(instance, self.name, default = self.default)
    
    def __set__(self, instance, value):
        """ If the value cannot be cast to the expected type, a ValueError is raised. """
        if self._readonly:
            raise AttributeError

        try:
            value = self.type(value)
        except ValueError:
            raise TypeError('Experimental parameter {} expects values of  \
                              type {}, but received {}'.format(self.name, self.ptype, value))
        else:
             return setattr(instance, self.name, value)

class RawDatasetBase(ABC):
    """
    Abstract base class for ultrafast electron diffraction data set based on filesystem. 
    RawDatasetBase allows for enforced metadata types and values, as well as a standard interface.

    Parameters
    ----------
    directory : str
        Location of the folder containing the raw data.
    metadata : dict
        Experimental parameters and metadata.
    
    Raises
    ------
    ValueError : if directory does not exist.
    TypeError : if an item from the metadata has an unexpected type.
    """
    # Function used to read raw data, e.g. TIFF images.
    # This function should take a single non-optional argument: a filename
    # By default, this is Scikit-image's imread.
    image_load_func = imread

    # List of valid metadata below
    # Using the ExperimentalParameter allows for automatic registering
    # of the parameters as valid.
    # These attributes can be accessed using the usual property access
    date =            ExperimentalParameter('date', str, default = '', readonly = True)
    pump_wavelength = ExperimentalParameter('pump_wavelength', int, default = 800, readonly = True)
    fluence =         ExperimentalParameter('fluence', float, default = 0)
    time_zero_shift = ExperimentalParameter('time_zero_shift', float, default = 0)
    temperature =     ExperimentalParameter('temperature', float, default = 293)
    exposure =        ExperimentalParameter('exposure', float, default = 1)
    resolution =      ExperimentalParameter('resolution', tuple, default = (2048, 2048), readonly = True)
    time_points =     ExperimentalParameter('time_points', tuple, default = tuple(), readonly = True)
    scans =           ExperimentalParameter('scans', tuple, default = (1,))
    notes =           ExperimentalParameter('notes', str, default = '')

    def __init__(self, directory, metadata, **kwargs):
        if not isdir(directory):
            raise ValueError('Directory {} does not exist.'.format(directory))
        self.directory = directory

        # Only record data that is valid
        # i.e.
        for k, v in metadata.items():
            if k in self.valid_metadata:
                setattr(self, k, v)
    
    def __repr__(self):
        string = '< RawDataset object. '
        for k, v in self.metadata.items():
            string.join('\n {key}: {value} '.format(key = k, value = v))
        return string
    
    @cached_property
    def valid_metadata(self):
        """ Iterable of valid parameter names """
        return {parameter.name for parameter in self.__dict__ if isinstance(parameter, ExperimentalParameter)}

    @property
    def metadata(self):
        """ Experimental parameters and dataset metadata as a dictionary. """
        # This property could be generated from some metadata file
        return {k:getattr(self, k) for k in self.valid_metadata}
        
    @abstractproperty
    def pumpoff_background(self):
        """ Image or composite of the dark background """
        pass
    
    @abstractproperty
    def pumpon_background(self):
        """ Image or composite of the background during photoexcitation """
        pass
    
    @abstractmethod
    def raw_data_filename(self, timedelay, scan, **kwargs):
        """
        Filename of the raw image.

        Parameters
        ----------
        timdelay : float
            Acquisition time-delay
        scan : int
            Scan number.
        """
        pass
    
    def raw_data(self, timedelay, scan, **kwargs):
        """
        Returns an array of the image at a timedelay and scan.
        
        Parameters
        ----------
        timdelay : float
            Acquisition time-delay.
        scan : int
            Scan number.
        
        Returns
        -------
        arr : `~numpy.ndarray`, ndim 2
        
        Raises
        ------
        IOError : Filename is not associated with an image/does not exist.
        """ 
        return self.image_load_func(self.raw_data_filename(timedelay, scan, **kwargs))