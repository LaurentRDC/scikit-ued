# -*- coding: utf-8 -*-
"""
Handling UED Datasets
=====================
"""
from abc import ABC, ABCMeta, abstractmethod, abstractproperty
from contextlib import suppress

from cached_property import cached_property

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
    __slots__ = ('name', 'type', 'default', 'readonly')

    def __init__(self, name, ptype, default = None, readonly = False):
        self.name = name
        self.type = ptype
        self.default = default
        self.readonly = readonly
    
    def __get__(self, instance, cls):
        if instance is None:
            return self
        return instance.__dict__.get(self.name, self.default)
    
    def __set__(self, instance, value):
        """ If the value cannot be cast to the expected type, a ValueError is raised. """
        if self.readonly:
            raise AttributeError

        try:
            value = self.type(value)
        except ValueError:
            raise TypeError('Experimental parameter {} expects values of type \
                             {}, but received {}'.format(self.name, self.type, value))
        else:
             instance.__dict__[self.name] = value

class MetaRawDataset(ABCMeta):
    """
    Metaclass for AbstractRawDataset. 
    
    This metaclass allows to determine the valid metadata that has been defined using the
    ExperimentalParameter class descriptor as class variables. For example, the AbstractRawDataset
    class already has some built-in ExperimentalParameter descriptors (date, notes, etc.)
    """
    def __init__(cls, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if not hasattr(cls, 'valid_metadata'):
            cls.valid_metadata = set([])
        
        # valid metadata as defined on the local class
        # Only metadata defined via the ExperimentalParameter descriptor will appear in
        # instance.metadata
        local_valid_metadata = {name for name, parameter in cls.__dict__.items() 
                                if isinstance(parameter, ExperimentalParameter)}
        cls.valid_metadata = cls.valid_metadata.union(local_valid_metadata)

        # If available, also include valid metadata from superclasses
        with suppress(AttributeError):
            cls.valid_metadata = set.union(cls.valid_metadata, super().valid_metadata) 

class AbstractRawDataset(metaclass = MetaRawDataset):
    """
    Abstract base class for ultrafast electron diffraction data set based on filesystem. 
    RawDatasetBase allows for enforced metadata types and values, as well as a standard interface.

    Minimally, the following method must be implemented in subclasses:

        * raw_data
    
    Raises
    ------
    ValueError : if directory does not exist.
    TypeError : if an item from the metadata has an unexpected type.
    """

    # List of valid metadata below
    # Using the ExperimentalParameter allows for automatic registering
    # of the parameters as valid.
    # These attributes can be accessed using the usual property access
    pump_wavelength = ExperimentalParameter('pump_wavelength', int, default = 800, readonly = True)
    fluence =         ExperimentalParameter('fluence', float, default = 0)
    time_zero_shift = ExperimentalParameter('time_zero_shift', float, default = 0)
    temperature =     ExperimentalParameter('temperature', float, default = 293)
    exposure =        ExperimentalParameter('exposure', float, default = 1)
    resolution =      ExperimentalParameter('resolution', tuple, default = (2048, 2048), readonly = True)
    time_points =     ExperimentalParameter('time_points', tuple, default = tuple(), readonly = True)
    scans =           ExperimentalParameter('scans', tuple, default = (1,))
    notes =           ExperimentalParameter('notes', str, default = '')    
    
    def __repr__(self):
        string = '< RawDataset object. '
        for k, v in self.metadata.items():
            string.join('\n {key}: {value} '.format(key = k, value = v))
        string.join(' >')
        return string

    @property
    def metadata(self):
        """ Experimental parameters and dataset metadata as a dictionary. """
        # This property could be generated from some metadata file
        return {k:getattr(self, k) for k in self.valid_metadata}
    
    @abstractmethod
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
        pass
