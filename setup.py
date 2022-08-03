# -*- coding: utf-8 -*-
import re
from glob import glob
from itertools import chain
from pathlib import Path

# import numpy
# from Cython.Build import cythonize
from setuptools import find_packages, setup

PACKAGE_NAME = "scikit-ued"
DESCRIPTION = "Collection of algorithms and functions for ultrafast electron scattering"
URL = "http://scikit-ued.readthedocs.io"
DOWNLOAD_URL = "http://github.com/LaurentRDC/scikit-ued"
AUTHOR = "Laurent P. René de Cotret"
AUTHOR_EMAIL = "laurent.renedecotret@mail.mcgill.ca"
BASE_PACKAGE = "skued"

base_path = Path(__file__).parent
with open(base_path / BASE_PACKAGE / "__init__.py") as f:
    module_content = f.read()
    VERSION = (
        re.compile(r".*__version__ = \"(.*?)\"", re.S).match(module_content).group(1)
    )
    LICENSE = (
        re.compile(r".*__license__ = \"(.*?)\"", re.S).match(module_content).group(1)
    )

with open("README.md") as f:
    README = f.read()

with open("requirements.txt") as f:
    REQUIREMENTS = [line for line in f.read().split("\n") if len(line.strip())]

exclude = {"exclude": ["external*", "docs", "*cache"]}
PACKAGES = [
    BASE_PACKAGE + "." + x
    for x in find_packages(str(base_path / BASE_PACKAGE), **exclude)
]
if BASE_PACKAGE not in PACKAGES:
    PACKAGES.append(BASE_PACKAGE)


if __name__ == "__main__":
    setup(
        name=PACKAGE_NAME,
        description=DESCRIPTION,
        long_description=README,
        long_description_content_type="text/markdown",
        license=LICENSE,
        url=URL,
        download_url=DOWNLOAD_URL,
        version=VERSION,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        maintainer=AUTHOR,
        maintainer_email=AUTHOR_EMAIL,
        install_requires=REQUIREMENTS,
        extras_require={"diffshow": ["pyqtgraph>=0.12,<1", "PyQt5"]},
        keywords="ultrafast electron scattering",
        project_urls={
            "Documentation": "https://scikit-ued.readthedocs.io/",
            "Source": "https://github.com/LaurentRDC/scikit-ued",
        },
        python_requires=">=3.7",
        packages=PACKAGES,
        entry_points={"console_scripts": ["skued = skued.__main__:main"]},
        include_package_data=True,
        zip_safe=False,
        #        include_dirs = [numpy.get_include()],
        #        ext_modules = cythonize("skued/*/**.pyx",
        #                                 compiler_directives = {'language_level':3,
        #                                                        'boundscheck': False}),
        # list of possible classifiers:
        #  https://pypi.python.org/pypi?%3Aaction=list_classifiers
        classifiers=[
            "Environment :: Console",
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Physics",
        ],
    )
