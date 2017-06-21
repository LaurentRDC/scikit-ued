# -*- coding: utf-8 -*-
from glob import glob
from itertools import chain
import os
import re
from setuptools import setup, find_packages
#from Cython.Build import cythonize

# How to submit a package on PyPi:
# http://peterdowns.com/posts/first-time-with-pypi.html
#
# Test upload:
#   > python setup.py sdist upload -r pypitest
#
# Real upload:
#   > python setup.py sdist upload -r pypi

BASE_PACKAGE = 'skued'

wavelets = chain.from_iterable([glob('skued\\baseline\\data\\*.npy'), 
                                glob('skued\\baseline\\data\\*.npz')])

base_path = os.path.dirname(__file__)
with open(os.path.join(base_path, 'skued', '__init__.py')) as f:
    module_content = f.read()
    VERSION = re.compile(r'.*__version__ = \'(.*?)\'', re.S).match(module_content).group(1)
    LICENSE = re.compile(r'.*__license__ = \'(.*?)\'', re.S).match(module_content).group(1)


with open('README.rst') as f:
    readme = f.read()

with open('requirements.txt') as f:
    requirements = [line for line in f.read().split('\n') if len(line.strip())]

exclude = {'exclude': ['cif2cell*', 'pdb_cache', 'docs']}
packages = [BASE_PACKAGE + '.' + x for x in find_packages(os.path.join(base_path, BASE_PACKAGE), **exclude)]
if BASE_PACKAGE not in packages:
    packages.append(BASE_PACKAGE)

if __name__ == '__main__':
    setup(
        name = 'scikit-ued',
        description = 'Collection of algorithms and functions for ultrafast electron diffraction',
        long_description = readme,
        license = LICENSE,
        url = 'http://scikit-ued.readthedocs.io',
        download_url = 'http://github.com/LaurentRDC/scikit-ued',
        version = VERSION,
        author = 'Laurent P. René de Cotret',
        author_email = 'laurent.renedecotret@mail.mcgill.ca',
        maintainer = 'Laurent P. René de Cotret',
        maintainer_email = 'laurent.renedecotret@mail.mcgill.ca',
        install_requires = requirements,
        keywords = ['skued'],
        packages = packages,
        data_files = [('skued\\baseline\\data', wavelets)],
        include_package_data = True,
        zip_safe = False,
        classifiers = ['Environment :: Console',
                       'Intended Audience :: Science/Research',
                       'Topic :: Scientific/Engineering',
                       'Topic :: Scientific/Engineering :: Physics',
                       'License :: OSI Approved :: MIT License',
                       'Natural Language :: English',
                       'Operating System :: OS Independent',
                       'Programming Language :: Python',
                       'Programming Language :: Python :: 3.5',
                       'Programming Language :: Python :: 3.6']
    )
