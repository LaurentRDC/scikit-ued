# Copyright 2010 Torbjorn Bjorkman
# This file is part of cif2cell
#
# cif2cell is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cif2cell is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cif2cell.  If not, see <http://www.gnu.org/licenses/>.
#
from distutils.core import setup
from glob import glob
from subprocess import call
from os import remove
import sys
import shutil
import tarfile

# Check for PyCifRW and attempt to install bundled version if not available.
pycifrwmodulelist = []
try:
    import CifFile
except:
    try:
        print "PyCifRW not found. Installing relevant modules."
        pycifrwmodulelist = ['CifFile','StarFile','YappsStarParser_1_0','YappsStarParser_1_1',
                             'YappsStarParser_DDLm','yapps3_compiled_rt']
        tar = tarfile.open('PyCifRW-3.3.tar.gz')
        tar.extractall()
        for mod in pycifrwmodulelist:
            shutil.copy("PyCifRW-3.3/"+mod+".py","./")
        shutil.rmtree('PyCifRW-3.3')
    except:
        print "PyCifRW could neither be found nor installed. Please install PyCIFRW manually, "+\
              "either from the\ncopy provided in this package or from: "+\
              "https://sourceforge.net/projects/pycifrw.berlios/\n"+\
              "(or maybe you have it installed and just need to adjust your PYTHONPATH)."
        sys.exit(1)

# Set up documentation
docfiles = ['docs/cif2cell.pdf']

# Get list of the cif example files
ciffiles = glob('cifs/*.cif')
periodiccifs = glob('cifs/periodic_table/*.cif')+['cifs/periodic_table/README']

# python modules
modulelist = pycifrwmodulelist+['utils','uctools','spacegroupdata','elementdata','ESPInterfaces']

setup(name='cif2cell',
      version='1.2.10',
      description='Construct a unit cell from CIF data',
      long_description='A command-line tool to generate the geometrical setup for various electronic structure codes from a CIF format file.',
      author='Torbjorn Bjorkman',
      author_email='torbjornb@gmail.com',
      url='http://cif2cell.sourceforge.net/',
      py_modules=modulelist,
      scripts=['cif2cell'],
      requires=['CifFile'],
      data_files=[('lib/cif2cell', ['LICENSE','HOWTOCITE']),
                  ('lib/cif2cell/sample_cifs', ciffiles),
                  ('lib/cif2cell/sample_cifs/periodic_table', periodiccifs),
                  ('lib/cif2cell/docs',docfiles)],
      license='GNU General Public License version 3'
      )

# cleanup
for mod in pycifrwmodulelist:
    remove("./"+mod+".py")
