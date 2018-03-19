# -*- coding: utf-8 -*-

import numpy as np
from itertools import chain

class Base:
    """ 
    Base class that overrides the builtin behavior:
    >>> object() == object()
    False

    This allows for transparent multiple inheritance of subclasses.
    """
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return True
        return NotImplemented

    def __hash__(self):
        return 0

class AtomicStructure(Base):
    """
    Base class for atomic structures. These structures can be made
    out of :class:`Atom` objects, or other AtomicStructure subclasses.

    The AtomicStructure class provides an abstraction over structure with and without
    substructure. Subclasses can be iterated over as an iterable of atoms. Order of iteration
    is not guaranteed.
    
    Hierarchical containership of AtomicStructures is implemented, as well as 
    containership of :class:`Atom` instances.

    Parameters
    ----------
    atoms : iterable, optional
        Atoms not attached to a substructure.
    substructures : iterable, optional
        Atomic substructures. For example, different atom chains could form
        a secondary structure in a protein.
    """

    def __init__(self, atoms = tuple(), substructures = tuple(), **kwargs):
        self.atoms = frozenset(atoms)
        self.substructures = frozenset(substructures)
        super().__init__(**kwargs)
    
    def __iter__(self):
        """ Yields :class:`Atom` instances from the structure and substructures 
        recursively. Order is not guaranteed. """
        yield from iter(self.atoms)
        yield from chain(*self.substructures)
    
    def itersorted(self, *, key = None, reverse = False):
        """ 
        Yields :class:`Atom` in sorted order. By default, atoms are sorted by element. 

        Parameters
        ----------
        key : callable or None, optional 
            Function of one argument that is used to extract a comparison key for an :class:`Atom` instance.
        reverse : bool, optional
            If True, elements are yielded in reverse order.
        
        Yields
        ------
        atm : skued.Atom 
        """
        if key is None:
            key = lambda atm: atm.element
        yield from sorted(iter(self), key = key, reverse = reverse)

    def __contains__(self, item):
        """ Check containership of :class:`Atom` instances or :class:`AtomicStructure` substructures recursively."""
        if isinstance(item, AtomicStructure):
            return (item in self.substructures)
        
        # Either the item is an orphan atom or 
        # it is in one of the substructures
        # Checking containership of sets is faster than iterating
        if item in self.atoms:
            return True
        else:
            return any((item in struct) for struct in self.substructures)
    
    def __len__(self):
        """ Number of :class:`Atom` instances present in the structure and substructures """
        return len(self.atoms) + sum(len(struct) for struct in self.substructures)
    
    def __hash__(self):
        return hash((self.atoms, self.substructures)) | super().__hash__()
    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (set(self.atoms) == set(other.atoms) and 
                    set(self.substructures) == set(other.substructures) and
                    super().__eq__(other))
        return NotImplemented

    def __repr__(self):
        """ Verbose string representation of this instance. """
        # AtomicStructure subclasses need not override this method
        # since the class name is dynamically determined
        rep = '< {clsname} object with following orphan atoms:'.format(clsname = self.__class__.__name__)

        # Note that repr(Atom(...)) includes these '< ... >'
        # We remove those for cleaner string representation
        for atm in self.itersorted():
            rep += '\n    ' + repr(atm).replace('<', '').replace('>', '').strip()

        if self.substructures:
            rep += 'and the following substructures:'
            for struct in self.substructures:
                rep += '\n' + repr(struct)
        
        return rep + ' >'

    def __array__(self, *args, **kwargs):
        """ Returns an array in which each row represents an :class:`Atom` instance. Order is not guaranteed. """
        arr = np.empty(shape = (len(self), 4), *args, **kwargs)
        for row, atm in enumerate(self):
            arr[row, 0] = atm.atomic_number
            arr[row, 1:] = atm.coords
        return arr