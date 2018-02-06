# -*- coding: utf-8 -*-
"""
CIF 1.0, 1.1 and 2.0 parser based on cif2cell.

References
----------
.. [CIF2CELL] Torbjorn Bjorkman, "CIF2Cell: Generating geometries for electronic structure programs", 
                 Computer Physics Communications 182, 1183-1186 (2011) doi: 10.1016/j.cpc.2011.01.013
""" 
import warnings
from contextlib import suppress
from functools import lru_cache
from re import sub
from string import digits, punctuation

import numpy as np
from CifFile import ReadCif, get_number_with_esd
from numpy.linalg import inv

from . import Atom, Lattice, ParseError, frac_coords
from .. import affine_map, transform
from .spg_data import HM2Hall, Number2Hall, SymOpsHall


def sym_ops(equiv_site):
    """ Parse a symmetry operator from an equivalent-site representation 
    
    Parameters
    ----------
    equiv_site : str or iterable of strings
        Either comma-separated string e.g. "+y, +x, -z + 1/2" or an
        iterable of the comma-separated values, e.g. ["+y", "+x", "-z + 1/2"] 
    
    Returns
    -------
    sym_ops : ndarray, shape (4,4)
        Symmetry operator as a 4x4 affine transformation on the FRACTIONAL
        coordinates.
    """
    symmetry = np.zeros( (3,3) )
    translation = np.zeros( (3,) )

    if isinstance(equiv_site, str):
        equiv_site = equiv_site.split(',')
        
    equiv_site = tuple(map(lambda s: s.strip().lower(), equiv_site))
    for j in range(3):
        xyz = equiv_site[j].replace('+',' +').replace('-',' -').split()
        for i in xyz:
            if i.strip("+-") == 'x':
                symmetry[0, j] = float(i.strip('x')+"1")
            elif i.strip("+-") == 'y':
                symmetry[1, j] = float(i.strip('y')+"1")
            elif i.strip("+-") == 'z':
                symmetry[2, j] = float(i.strip('z')+"1")
            
            if i.strip("+-xyz") != "":
                translation[j] = eval(i)
    
    symmetry[:] = np.transpose(symmetry)
    
    # Combination of transform and translation into a single matrix
    # is done in a 4x4 affine transform
    symmetry_operation = affine_map(symmetry)
    symmetry_operation[:3,3] = translation
    return symmetry_operation

class CIFParser:
    """
    Collection of methods that parses CIF files based on cif2cell. The preferred method
    of using this object is as a context manager.

    Parameters
    ----------
    filename : str or path-like
        Location of the CIF file.
    
    References
    ----------
    .. [#] Torbjorn Bjorkman, "CIF2Cell: Generating geometries for electronic structure programs", 
           Computer Physics Communications 182, 1183-1186 (2011) doi: 10.1016/j.cpc.2011.01.013
    """
    def __init__(self, filename, **kwargs):
        # ReadCIF would get confused between local files and URLs
        # Therefore, more clear to pass an open file
        self._handle = open(filename, mode = 'r')
        self.file = ReadCif(self._handle, **kwargs)
    
    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        self._handle.close()

    @property
    def _first_block(self):
        return self.file[self.file.keys()[0]]
    
    @lru_cache(maxsize = 1)
    def hall_symbol(self):
        """ Returns the Hall symbol """
        block = self._first_block

        hall_symbol = block.get('_symmetry_space_group_name_Hall') or block.get('_space_group_name_Hall')

        # In some rare cases, the given hall symbol in the file isn't standard,
        # otherwise it would be a key in SymOpsHall
        # Then, it is preferable to infer the conventional hall symbol from other info
        if (hall_symbol is None) or (hall_symbol not in SymOpsHall):
            h_m_symbol = block.get('_symmetry_space_group_name_H-M') or block.get('_space_group_name_H-M_alt')
            
            if h_m_symbol is not None:
                h_m_symbol = sub('\s+', '', h_m_symbol)
                with suppress(KeyError):    # Symbol could be meaningless, e.g. h_m_symbol = '?' (True story)
                    hall_symbol =  HM2Hall[h_m_symbol]
        
        # Again, if hall_symbol is still missing OR invalid
        if (hall_symbol is None) or (hall_symbol not in SymOpsHall):
            table_number = block.get('_symmetry_Int_Tables_number') or block.get('_space_group_IT_number')
                
            if table_number is not None:
                hall_symbol = Number2Hall[table_number]
                
        if hall_symbol is None:
            raise ParseError('Hall symbol could not be inferred')
    
        if hall_symbol[0] == "-":
            hall_symbol = "-" + hall_symbol[1].upper() + hall_symbol[2:].lower()
        else:
            hall_symbol = hall_symbol[0].upper() + hall_symbol[1:].lower()

        return hall_symbol

    @lru_cache(maxsize = 1)
    def lattice_parameters(self):
        """ 
        Returns the lattice parameters associated to a CIF structure.

        Returns
        ----------
        a, b, c : float
            Lengths of lattice vectors [Angstroms]
        alpha, beta, gamma : float
            Angles of lattice vectors [degrees]. 
        """
        block = self._first_block

        try:
            a, _ = get_number_with_esd(block["_cell_length_a"])
            b, _ = get_number_with_esd(block["_cell_length_b"])
            c, _ = get_number_with_esd(block["_cell_length_c"])
            alpha, _ = get_number_with_esd(block["_cell_angle_alpha"])
            beta, _ = get_number_with_esd(block["_cell_angle_beta"])
            gamma, _ = get_number_with_esd(block["_cell_angle_gamma"])
        except:
            raise ParseError('Lattice vectors could not be determined.')

        return a, b, c, alpha, beta, gamma

    @lru_cache(maxsize = 1)
    def lattice_vectors(self):
        """ 
        Returns the lattice vectors associated to a CIF structure.
        
        Returns
        -------
        lv : list of ndarrays, shape (3,)
        """
        return Lattice.from_parameters(*self.lattice_parameters()).lattice_vectors
    
    def symmetry_operators(self):
        """
        Returns the symmetry operators that map the fractional atomic positions in a
        CIF file to the crystal *conventional* unit cell.

        Yields
        ------
        sym_ops : ndarray, shape (4,4)
            Transformation matrices. Since translations and rotation are combined,
            the transformation matrices are 4x4.
        """
        block = self._first_block

        equivalent_sites_str = None
        for tag in ['_symmetry_equiv_pos_as_xyz','_space_group_symop_operation_xyz']:
            with suppress(KeyError):
                equivalent_sites_str = block.GetLoop(tag).get(tag)

        # P1 space group only has a single equivalent site
        if isinstance(equivalent_sites_str, str):
            equivalent_sites_str = [equivalent_sites_str]

        with suppress(ParseError):
            if not equivalent_sites_str:
                equivalent_sites_str = SymOpsHall[self.hall_symbol()]
            elif len(equivalent_sites_str) != len(SymOpsHall[self.hall_symbol()]):
                warnings.warn('The number of equivalent sites is not in line with the database. The file might be incomplete')

        yield from map(sym_ops, equivalent_sites_str)
    
    def atoms(self):
        """
        Asymmetric unit cell. Combine with CIFParser.symmetry_operators() for a full unit cell.

        Yields
        ------
        atoms : skued.Atom instance
        """
        block = self._first_block

        try:
            tmpdata = block.GetLoop('_atom_site_fract_x')
            cartesian = False
        except:
            try:
                tmpdata = block.GetLoop('_atom_site_Cartn_x')
                cartesian = True
            except:
                raise ParseError('Atomic positions could not be found or inferred.')
            
            t11 = block.get('_atom_sites_Cartn_tran_matrix_11')
            t12 = block.get('_atom_sites_Cartn_tran_matrix_12')
            t13 = block.get('_atom_sites_Cartn_tran_matrix_13')
            t21 = block.get('_atom_sites_Cartn_tran_matrix_21')
            t22 = block.get('_atom_sites_Cartn_tran_matrix_22')
            t23 = block.get('_atom_sites_Cartn_tran_matrix_23')
            t31 = block.get('_atom_sites_Cartn_tran_matrix_13')
            t32 = block.get('_atom_sites_Cartn_tran_matrix_23')
            t33 = block.get('_atom_sites_Cartn_tran_matrix_33')
            cart_trans_matrix_inv = np.array([[float(t11),float(t12),float(t13)],
                                              [float(t21),float(t22),float(t23)],
                                              [float(t31),float(t32),float(t33)]])
            cart_trans_matrix = inv(cart_trans_matrix_inv)

            if not all([t11, t12, t13, t21, t22, t23, t31, t32, t33]):
                raise ParseError('Cartesian coordinates in CIF but no transformation matrix given')
        
        if cartesian:
            xs = tmpdata.get('_atom_site_Cartn_x')
            ys = tmpdata.get('_atom_site_Cartn_y')
            zs = tmpdata.get('_atom_site_Cartn_z')
        else:
            xs = tmpdata.get('_atom_site_fract_x')
            ys = tmpdata.get('_atom_site_fract_y')
            zs = tmpdata.get('_atom_site_fract_z')
        # TODO: handle wildcards like '?', '.' in xs, ys, zs

        elements = tmpdata.get('_atom_site_type_symbol')
        if not elements:
            elements = tmpdata.get('_atom_site_label')
            if not elements:
                raise ParseError('Atom symbols could not be found or inferred.')
        elements = map(lambda s: s.strip(punctuation + digits).title(), elements)
        
        lv = self.lattice_vectors()
        for e, x, y, z in zip(elements, xs, ys, zs):
            coords = np.array([get_number_with_esd(x)[0], 
                               get_number_with_esd(y)[0], 
                               get_number_with_esd(z)[0]])
            
            # We normalize atom position to be within the unit cell
            # Therefore we need the fractional coordinates
            if cartesian:
                coords = transform(cart_trans_matrix, coords)
                coords[:] = frac_coords(coords, lv)
            
            yield Atom(element = e, coords = np.mod(coords, 1))
