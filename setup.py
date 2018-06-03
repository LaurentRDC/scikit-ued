# -*- coding: utf-8 -*-
import re
from glob import glob
from itertools import chain
from pathlib import Path
from unittest import TestLoader

from setuptools import find_packages, setup

#from Cython.Build import cythonize

# To upload to pypi.org:
#   >>> python setup.py sdist
#   >>> twine upload dist/scikit-ued-x.x.x.tar.gz

PACKAGE_NAME    = 'scikit-ued'
DESCRIPTION     = 'Collection of algorithms and functions for ultrafast electron diffraction'
URL             = 'http://scikit-ued.readthedocs.io'
DOWNLOAD_URL    = 'http://github.com/LaurentRDC/scikit-ued'
AUTHOR          = 'Laurent P. René de Cotret'
AUTHOR_EMAIL    = 'laurent.renedecotret@mail.mcgill.ca'
BASE_PACKAGE    = 'skued'

WAVELET_FILES   = chain.from_iterable([glob('skued\\baseline\\data\\*.npy'), 
                                       glob('skued\\baseline\\data\\*.npz')])

CIF_FILES       = chain.from_iterable([glob('skued\\structure\\cifs\\*.cif')])

base_path = Path(__file__).parent
with open(base_path / BASE_PACKAGE / '__init__.py') as f:
    module_content = f.read()
    VERSION = re.compile(r'.*__version__ = \'(.*?)\'', re.S).match(module_content).group(1)
    LICENSE = re.compile(r'.*__license__ = \'(.*?)\'', re.S).match(module_content).group(1)

with open('README.rst') as f:
    README = f.read()

with open('requirements.txt') as f:
    REQUIREMENTS = [line for line in f.read().split('\n') if len(line.strip())]

exclude = {'exclude': ['external*', 'docs', '*cache']}
PACKAGES = [BASE_PACKAGE + '.' + x for x in find_packages(str(base_path / BASE_PACKAGE), **exclude)]
if BASE_PACKAGE not in PACKAGES:
    PACKAGES.append(BASE_PACKAGE)

def skued_test_suite():
    return TestLoader().discover('.')

if __name__ == '__main__':
    setup(
        name = PACKAGE_NAME,
        description = DESCRIPTION,
        long_description = README,
        license = LICENSE,
        url = URL,
        download_url = DOWNLOAD_URL,
        version = VERSION,
        author = AUTHOR,
        author_email = AUTHOR_EMAIL,
        maintainer = AUTHOR,
        maintainer_email = AUTHOR_EMAIL,
        install_requires = REQUIREMENTS,
        keywords = 'ultrafast electron diffraction',
        project_urls = {
            'Documentation' : 'http://scikit-ued.readthedocs.io/en/master/',
            'Source'        : 'https://github.com/LaurentRDC/scikit-ued',
        },
        python_requires = '>=3.6',
        packages = PACKAGES,
        data_files = [('skued\\baseline\\data', WAVELET_FILES),
                      ('skued\\structure\\cifs', CIF_FILES)],
        include_package_data = True,
        zip_safe = False,
        test_suite = 'setup.skued_test_suite', 

        # list of possible classifiers:
        #  https://pypi.python.org/pypi?%3Aaction=list_classifiers
        classifiers = ['Environment :: Console',
                       'Development Status :: 5 - Production/Stable',
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: MIT License',
                       'Natural Language :: English',
                       'Operating System :: OS Independent',
                       'Programming Language :: Python',
                       'Programming Language :: Python :: 3.6',
                       'Topic :: Scientific/Engineering',
                       'Topic :: Scientific/Engineering :: Physics',]
    )
