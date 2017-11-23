# -*- coding: utf-8 -*-
import os
from copy import deepcopy as copy
from functools import lru_cache
from glob import glob
from urllib.request import urlretrieve

import numpy as np
from spglib import (get_error_message, get_spacegroup_type, 
                    get_symmetry_dataset, find_primitive)

from . import Atom, CIFParser, Lattice, PDBParser
from .. import affine_map

CIF_ENTRIES = glob(os.path.join(os.path.dirname(__file__), 'cifs', '*.cif'))

def symmetry_expansion(atoms, symmetry_operators):
    """
    Generate a set of unique atoms from an asymmetric cell and symmetry operators.

    Parameters
    ----------
    atoms : iterable of Atom
        Assymetric unit cell atoms. It is assumed that the atomic 
        coordinates are in fractional form.
    symmetry_operators : iterable of array_like
        Symmetry operators that generate the full unit cell.
    
    Yields
    ------
    Atom
    """
    # TODO: provide ability to reduce to primitive, niggli_reduce, etc.
    #       using spglib?
    uniques = set([])
    symmetry_operators = tuple(map(affine_map, symmetry_operators))

    for atm in atoms:
        for sym_op in symmetry_operators:
            new = copy(atm)
            new.transform(sym_op)
            new.coords[:] = np.mod(new.coords, 1)
            uniques.add(new)
    yield from uniques

class Crystal(Lattice):
    """
    This object is the basis for inorganic crystals such as VO2, 
    and protein crystals such as bR. 

    In addition to constructing the ``Crystal`` object yourself, other constructors
    are also available (and preferred):
    
    * ``Crystal.from_cif``: create an instance from a CIF file;
    
    * ``Crystal.from_pdb``: create an instance from a Protein Data Bank entry;
    
    * ``Crystal.from_database``: create an instance from the internal database of CIF files;
    
    * ``Crystal.from_cod``: create an instance from a Crystallography Open Database entry.

    * ``Crystal.from_ase``: create an instance from an ``ase.Atoms`` instance.

    Parameters
    ----------
    unitcell : iterable of ``Atom``
        Unit cell atoms. It is assumed that the atoms are in fractional coordinates.
    lattice_vectors : iterable of array_like
        Lattice vectors.
    source : str or None, optional
        Provenance, e.g. filename.
    """

    builtins = frozenset(map(lambda fn: os.path.basename(fn).split('.')[0], CIF_ENTRIES))

    def __init__(self, unitcell, lattice_vectors, source = None, **kwargs):
        self.unitcell = frozenset(unitcell)
        self.source = source
        super().__init__(lattice_vectors, **kwargs)
    
    def __iter__(self):
        yield from iter(self.unitcell)
    
    def __len__(self):
        return len(self.unitcell)
    
    def __repr__(self):
        return '< Crystal object with unit cell of {} atoms >'.format(len(self))
    
    def __eq__(self, other):
        return (isinstance(other, self.__class__) 
                and (set(self) == set(other))
                and super().__eq__(other))
    
    def __array__(self):
        """ Returns an array in which each row represents a unit cell atom """
        arr = np.empty(shape = (len(self), 4), dtype = np.float)
        for row, atm in enumerate(self):
            arr[row, 0] = atm.atomic_number
            arr[row, 1:] = atm.coords
        return arr

    @classmethod
    @lru_cache(maxsize = len(builtins)) # saves a lot of time in tests
    def from_cif(cls, path):
        """
        Returns a Crystal object created from a CIF 1.0, 1.1 or 2.0 file.

        Parameters
        ----------
        path : str
            File path
        
        References
        ----------
        .. [#] Torbjorn Bjorkman, "CIF2Cell: Generating geometries for electronic structure programs", 
               Computer Physics Communications 182, 1183-1186 (2011). doi: 10.1016/j.cpc.2011.01.013
        """
        with CIFParser(filename = path) as parser:
            return cls(unitcell = symmetry_expansion(parser.atoms(), parser.symmetry_operators()),
                       lattice_vectors = parser.lattice_vectors(),
                       source = str(path))
    
    @classmethod
    def from_database(cls, name):
        """ 
        Returns a Crystal object create from the internal CIF database.

        Parameters
        ----------
        name : str
            Name of tne databse entry. Available items can be retrieved from `Crystal.builtins`
        """
        if name not in cls.builtins:
            raise ValueError('Entry {} is not available in the database. See \
                              Crystal.builtins for valid entries.'.format(name))
        
        path = os.path.join(os.path.dirname(__file__), 'cifs', name + '.cif')
        return cls.from_cif(path)
    
    @classmethod
    def from_cod(cls, num, revision = None, download_dir = 'cod_cache', overwrite = False):
        """ 
        Returns a Crystal object built from the Crystallography Open Database. 

        Parameters
        ----------
        num : int
            COD identification number.
        revision : int or None, optional
            Revision number. If None (default), the latest revision is used.
        download_dir : path-like object, optional
            Directory where to save the CIF file. Default is a local folder in the current directory
        overwrite : bool, optional
            Whether or not to overwrite files in cache if they exist. If no revision 
            number is provided, files will always be overwritten. 
        """
        if revision is None:
            overwrite = True
        
        if not os.path.isdir(download_dir):
            os.mkdir(download_dir)
        
        url = 'http://www.crystallography.net/cod/{}.cif'.format(num)

        if revision is not None:
            url = url + '@' + str(revision)
            base = '{iden}-{rev}.cif'.format(iden = num, rev = revision)
        else:
            base = '{}.cif'.format(num)
        path = os.path.join(download_dir, base)

        if (not os.path.isfile(path)) or overwrite:
            urlretrieve(url, path)
        
        return cls.from_cif(path)

    @classmethod
    def from_pdb(cls, ID, download_dir = 'pdb_cache', overwrite = False):
        """
        Returns a Crystal object created from a Protein DataBank entry.

        Parameters
        ----------
        ID : str
            Protein DataBank identification. The correct .pdb file will be downloaded,
            cached and parsed.
        download_dir : path-like object, optional
            Directory where to save the PDB file. Default is a local folder in the current directory
        overwrite : bool, optional
            Whether or not to overwrite files in cache if they exist. If no revision 
            number is provided, files will always be overwritten. 
        """
        parser = PDBParser(ID = ID, download_dir = download_dir)
        return cls(unitcell = symmetry_expansion(parser.atoms(), parser.symmetry_operators()),
                   lattice_vectors = parser.lattice_vectors(),
                   source = str(parser.file))
    
    @classmethod
    def from_ase(cls, atoms):
        """
        Returns a Crystal object created from an ASE Atoms object.
        
        Parameters
        ----------
        atoms : ase.Atoms
            Atoms group.
        """
        lattice_vectors = atoms.get_cell()
        
        return cls(unitcell = [Atom.from_ase(atm) for atm in atoms], 
                   lattice_vectors = lattice_vectors)

    def _spglib_cell(self):
        """ Returns an array in spglib's cell format. """
        arr = np.asarray(self)
        return np.array(self.lattice_vectors), arr[:, 1:], arr[:, 0]

    def primitive(self, symprec = 1e-2):
        """ 
        Returns a Crystal object in the primitive unit cell.
        
        Parameters
        ----------
        symprec : float, optional
            Symmetry-search distance tolerance in Cartesian coordinates [Angstroms].

        Returns
        -------
        primitive : Crystal
            Crystal with primitive cell. If primitive cell is the same size as
            the source Crystal, a reference to the source Crystal is returned.

        Raises
        ------
        RuntimeError
            If primitive cell could not be found.
        
        Notes
        -----
        Optional atomic properties (e.g magnetic moment) might be lost in the reduction.
        """
        search = find_primitive(self._spglib_cell(), 
                                symprec = symprec)
        if search is None:
            raise RuntimeError('Primitive cell could not be found.')

        lattice_vectors, scaled_positions, numbers = search
        if numbers.size == len(self):   # Then there's no point in creating a new crystal
            return self

        atoms = list()
        for coords, Z in zip(scaled_positions, numbers):
            atoms.append(Atom(int(Z), coords = coords))

        return Crystal(unitcell = atoms, 
                       lattice_vectors = lattice_vectors, 
                       source = self.source)

    def ase_atoms(self, **kwargs):
        """ 
        Create an ASE Atoms object from a Crystal. 
        
        Parameters
        ----------
        kwargs
            Keyword arguments are passed to ase.Atoms constructor.
        
        Returns
        -------
        atoms : ase.Atoms
            Group of atoms ready for ASE's routines.
        
        Raises
        ------
        ImportError
            If ASE is not installed
        """
        from ase import Atoms
        
        return Atoms(symbols = [atm.ase_atom(lattice = self) for atm in iter(self)],
                     cell = np.array(self.lattice_vectors), **kwargs)
    
    def spacegroup_info(self, symprec = 1e-2, angle_tolerance = -1.0):
        """ 
        Returns a dictionary containing space-group information. This information
        is computed from the crystal unit cell, and is not taken from records if available.
        
        Parameters
        ----------
        symprec : float, optional
            Symmetry-search distance tolerance in Cartesian coordinates [Angstroms].
        angle_tolerance: float, optional
            Symmetry-search tolerance in degrees. If the value is negative (default), 
            an internally optimized routine is used to judge symmetry.
        
        Returns
        -------
        info : dict or None
            Dictionary of space-group information. The following keys are available:

            * ``'international_symbol'``: International Tables of Crystallography space-group symbol (short);

            * ``'international_full'``: International Tables of Crystallography space-group full symbol;

            * ``'hall_symbol'`` : Hall symbol;

            * ``'pointgroup'`` : International Tables of Crystallography point-group;

            * ``'international_number'`` : International Tables of Crystallography space-group number (between 1 and 230);

            * ``'hall_number'`` : Hall number (between 1 and 531).

            If symmetry-determination has failed, None is returned.
        
        Raises
        ------
        RuntimeError
            If symmetry-determination has yielded an error.
        
        Notes
        -----
        Note that crystals generated from the Protein Data Bank are often incomplete; 
        in such cases the space-group information will be incorrect.
        """
        dataset = get_symmetry_dataset(cell = self._spglib_cell(),
                                       symprec = symprec, 
                                       angle_tolerance = angle_tolerance)

        if dataset: 
            info = dict()
            info.update( {'international_symbol': dataset['international'],
                          'hall_symbol': dataset['hall'],
                          'international_number': dataset['number'],
                          'hall_number': dataset['hall_number']} )
            
            spg_type = get_spacegroup_type(info['hall_number'])
            info.update( {'international_full': spg_type['international_full'],
                          'pointgroup': spg_type['pointgroup_international']} )

            return info

        err_msg = get_error_message()
        if err_msg:
            raise RuntimeError('Symmetry-determination has returned the following error: {}'.format(err_msg))
