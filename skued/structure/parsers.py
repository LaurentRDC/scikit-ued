# -*- coding: utf-8 -*-
"""
CIF 1.0, 1.1 and 2.0 parser based on cif2cell.

References
----------
.. [CIF2CELL] Torbjorn Bjorkman, "CIF2Cell: Generating geometries for electronic structure programs", 
                 Computer Physics Communications 182, 1183-1186 (2011) doi: 10.1016/j.cpc.2011.01.013
""" 
import gzip
import os
import warnings
from abc import abstractmethod
from contextlib import AbstractContextManager, suppress
from functools import lru_cache
from re import sub
from string import digits, punctuation
from urllib.request import urlretrieve

import numpy as np
from CifFile import ReadCif, get_number_with_esd
from numpy.linalg import inv

from . import Atom, Lattice, ParseError, frac_coords
from .. import affine_map, transform
from .spg_data import HM2Hall, Number2Hall, SymOpsHall


class AbstractStructureParser(AbstractContextManager):
    """
    Abstract base class for structure parsers. The preferred method
    of using this object is as a context manager.

    Parameters
    ----------
    filename : str or path-like
        Location of the CIF file.
    """
    @abstractmethod
    def __init__(self, *args, **kwargs):
        pass
    
    @abstractmethod
    def lattice_vectors(self):
        """ 
        Returns the lattice vectors associated to a structure.
        
        Returns
        -------
        lv : iterable of ndarrays, shape (3,)
        """
        pass
    
    @abstractmethod
    def symmetry_operators(self):
        """
        Returns the symmetry operators that map the fractional atomic positions in a
        structure to the crystal *conventional* unit cell.

        Yields
        ------
        sym_ops : ndarray, shape (4,4)
            Transformation matrices. Since translations and rotation are combined,
            the transformation matrices are 4x4.
        """
        pass
    
    @abstractmethod
    def atoms(self):
        """
        Asymmetric unit cell. Combine with CIFParser.symmetry_operators() for a full unit cell.

        Yields
        ------
        atoms : skued.Atom instance
        """
        pass


class PDBParser(AbstractStructureParser):
    """
    Collection of methods that parses PDB files. This object should be used as a context manager.
    
    Parameters
    ----------
    ID : str
        Protein DataBank identification. The correct .pdb file will be downloaded,
        cached and parsed.
    download_dir : path-like object
        Directory where to save the PDB file. Default is a local folder in the current directory
    overwrite : bool, optional
        Whether or not to overwrite files in cache if they exist. If no revision 
        number is provided, files will always be overwritten. 
    """

    def __init__(self, ID, download_dir = 'pdb_cache', overwrite = False):
        filename = self.retrieve_pdb_file(pdb_code = ID, 
                                          download_dir = download_dir, 
                                          overwrite = overwrite)
        self._handle = open(filename, 'r')
    
    def __exit__(self, *args, **kwargs):
        self._handle.close()

    @staticmethod
    def retrieve_pdb_file(pdb_code, download_dir = None, server = 'ftp://ftp.wwpdb.org', overwrite = False):
        """ 
        Retrieves a PDB structure file from the PDB server and
        stores it in a local file tree.

        Parameters
        ----------
        pdf_code : str, len 4
            PDB ID code
        download_dir : path-like object
            Directory where to save the PDB file. Default is a local folder in the current directory
        server : str, optional
            Address of the FTP server from which to download the PDB file. Default is the main server.
        overwrite : bool, optional
            If True, existing PDB file with the same structure will be overwritten. Default is False.

        Returns
        -------
        file : str 
            Pointer to the downloaded file
        """
        # Get the compressed PDB structure
        code = pdb_code.lower()
        archive_fn = "pdb{}.ent.gz".format(code)
        pdb_dir = "divided"
        url = (server + '/pub/pdb/data/structures/{}/pdb/{}/{}'.format(pdb_dir, code[1:3], archive_fn))

        # Where does the final PDB file get saved?
        if download_dir is None:
            path = os.path.join(os.getcwd(), code[1:3])
        else:
            path = download_dir
        if not os.access(path, os.F_OK):
            os.makedirs(path)

        filename = os.path.join(path, archive_fn)
        final_file = os.path.join(path, "pdb{}.ent".format(code))  # (decompressed)

        # Skip download if the file already exists
        if (not overwrite) and (os.path.exists(final_file)):
                return final_file

        urlretrieve(url, filename)

        # Uncompress the archive, delete when done
        # Can't use context manager with gzip.open until Python 2.7
        with gzip.open(filename, 'rb') as gz:
            with open(final_file, 'wb') as out:
                out.writelines(gz)
        os.remove(filename)

        return final_file
    
    @property
    def filename(self):
        return self._handle.name

    @lru_cache(maxsize = 1)
    def lattice_vectors(self):
        """ 
        Returns the lattice vectors associated to a PDB structure.
        
        Returns
        -------
        lv : list of ndarrays, shape (3,)
        
        Raises
        ------
        ParseError
            If the file does not contain a CRYST1 tag.
        """
        self._handle.seek(0)

        for line in filter(lambda l: l.startswith('CRYST1'), self._handle):
            # characters are described in the PDB file content guide.
            a, b, c = float(line[6:15]), float(line[15:24]), float(line[24:33])
            alpha, beta, gamma = float(line[33:40]), float(line[40:47]), float(line[47:54])
            break
        else:
            raise ParseError('No CRYST1 line found')

        return Lattice.from_parameters(a, b, c, alpha, beta, gamma).lattice_vectors
    
    def atoms(self):
        """
        Returns a list of atoms associated with a PDB structure. These atoms form the asymmetric unit cell.

        Yields
        ------
        atom: skued.structure.Atom
        """
        self._handle.seek(0)

        atoms = list()
        for line in filter(lambda l: l.startswith( ('ATOM', 'HEMATM') ), self._handle):
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            element = str(line[76:78]).replace(' ','')
            yield Atom(element = element, coords = frac_coords(np.array([x,y,z]), self.lattice_vectors()))
    
    def symmetry_operators(self):
        """
        Returns the symmetry operators that map the atomic positions in a
        PDB file to the crystal unit cell.

        Yields
        ------
        sym_ops : `~numpy.ndarray`, shape (4,4)
            Transformation matrices. Since translations and rotation are combined,
            the transformation matrices are 4x4.
        """
        self._handle.seek(0)

        # This is clunky af
        sym_ops = dict()
        for line in filter(lambda l: l.startswith('REMARK 290') and ('SMTRY' in l), self._handle):
            
            op_num = line[22:23]
            if op_num not in sym_ops:
                sym_ops[op_num] = {'rotation': list(), 'translation': list()}
                
            r1, r2, r3, t = np.fromstring(line[23:], dtype = np.float, count = 4, sep = ' ')
            sym_ops[op_num]['rotation'].append([r1,r2,r3])
            sym_ops[op_num]['translation'].append(t)

        if not sym_ops:
            raise ParseError('No symmetry could be parsed from file {}'.format(self._handle.filename))
        
        for op in sym_ops.values():
            mat = np.eye(4, dtype = np.float)
            mat[:3,:3] = np.array(op['rotation'])
            mat[:3,3] = np.array(op['translation'])

            yield mat


class CIFParser(AbstractStructureParser):
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

    def __exit__(self, *args, **kwargs):
        self._handle.close()

    @staticmethod
    def sym_ops_from_equiv(equiv_site):
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
    
    @property
    def structure_block(self):
        """ Retrieve which CIF block has the appropriate structural information """
        blocks = (self.file[key] for key in self.file.keys())
        for block in blocks:
            try:
                a, _ = get_number_with_esd(block["_cell_length_a"])
            except:
                continue
            else:
                return block
    
    @lru_cache(maxsize = 1)
    def hall_symbol(self):
        """ Returns the Hall symbol """
        block = self.structure_block

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
        block = self.structure_block

        try:
            a, _ = get_number_with_esd(block["_cell_length_a"])
            b, _ = get_number_with_esd(block["_cell_length_b"])
            c, _ = get_number_with_esd(block["_cell_length_c"])
            alpha, _ = get_number_with_esd(block["_cell_angle_alpha"])
            beta, _ = get_number_with_esd(block["_cell_angle_beta"])
            gamma, _ = get_number_with_esd(block["_cell_angle_gamma"])
        except:
            raise ParseError('Lattice vectors could not be determined.')
        else:
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
        block = self.structure_block

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

        yield from map(self.sym_ops_from_equiv, equivalent_sites_str)
    
    def atoms(self):
        """
        Asymmetric unit cell. Combine with CIFParser.symmetry_operators() for a full unit cell.

        Yields
        ------
        atoms : skued.Atom instance
        """
        block = self.structure_block

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

class CODParser(CIFParser):
    """
    Collection of methods that parses CIF files retrieved from the Crystallography Open Database. 
    The preferred method of using this object is as a context manager.

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
    def __init__(self, num, revision = None, download_dir = 'cod_cache', overwrite = False, **kwargs):
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
        
        return super().__init__(filename = path, **kwargs)