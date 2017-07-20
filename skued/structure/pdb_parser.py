# -*- coding: utf-8 -*-
import gzip
import os
from functools import lru_cache
from urllib.request import urlretrieve

import numpy as np

from . import Atom, frac_coords, lattice_vectors_from_parameters, ParseError

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

class PDBParser(object):
    """
    Collection of methods that parses PDB files. Insipred from BioPython.PDB module.
    
    Attributes
    ----------
    ID : str, len 4
        Unique PDB ID of the parser's target structure.
    file : str
        Absolute path to the PDB file associated with the parser's target structure.

    Methods
    -------
    lattice_vectors
        Lattice vectors of the crystal structure.

    Atoms
        Yields atoms making up the asymmetric unit cell

    symmetry_operators
        Symmetry operators that relate the asymmetric unit cell and the unit cell. The symmetry operators
        act on all ATM and HETATM in the file to produce the unit cel.
    """

    def __init__(self, ID, download_dir = 'pdb_cache'):
        self.file = retrieve_pdb_file(pdb_code = ID, download_dir = download_dir)

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
        with open(self.file) as pdb_file:
            for line in pdb_file:
                if line.startswith('CRYST1'):
                    # characters are described in the PDB file content guide.
                    a, b, c = float(line[6:15]), float(line[15:24]), float(line[24:33])
                    alpha, beta, gamma = float(line[33:40]), float(line[40:47]), float(line[47:54])
                    break
            else:
                raise ParseError('No CRYST1 line found')

        return lattice_vectors_from_parameters(a, b, c, alpha, beta, gamma)
    
    def atoms(self):
        """
        Returns a list of atoms associated with a PDB structure. These atoms form the asymmetric unit cell.

        Yields
        ------
        atom: skued.structure.Atom
        """
        atoms = list()
        with open(self.file) as pdb_file:
            for line in filter(lambda l: l.startswith( ('ATOM', 'HEMATM') ), pdb_file):
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
        # This is clunky af
        sym_ops = dict()
        with open(self.file) as pdb_file:
            for line in filter(lambda l: l.startswith('REMARK 290') and ('SMTRY' in l), pdb_file):
                
                op_num = line[22:23]
                if op_num not in sym_ops:
                    sym_ops[op_num] = {'rotation': list(), 'translation': list()}
                    
                r1, r2, r3, t = np.fromstring(line[23:], dtype = np.float, count = 4, sep = ' ')
                sym_ops[op_num]['rotation'].append([r1,r2,r3])
                sym_ops[op_num]['translation'].append(t)

        if not sym_ops:
            raise ParseError('No symmetry could be parsed from file {}'.format(self.file))
        
        for op in sym_ops.values():
            mat = np.eye(4, dtype = np.float)
            mat[:3,:3] = np.array(op['rotation'])
            mat[:3,3] = np.array(op['translation'])

            yield mat
