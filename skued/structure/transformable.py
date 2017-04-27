class Transformable(object):
    """
    This class defines objects that can be transformed.

    Methods
    -------
    transform
    
    rotate

    translate
    """
    # We need to have slots in Transformable
    # so that subclasses can have their own __slots__
    # E,g, Atom has __slots__, but Lattice doesn't
    # This is only possible if Transformable has __slots__
    __slots__ = tuple()
    
    def transform(self, *args, **kwargs):
        pass
    
    def rotate(self, *args, **kwargs):
        pass

    def translate(self, *args, **kwargs):
        pass