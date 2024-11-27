# setup.py
import sysconfig
import sys
import os
import subprocess
import re

# Import setuptools at the beginning
import setuptools
from setuptools import setup, find_packages
from distutils.version import LooseVersion
from Cython.Distutils import build_ext

# Ensure dependencies are installed before setup
try:
    import numpy
    import Cython
    import pysam
except ImportError as e:
    missing_package = str(e).split("'")[1]
    raise ImportError(f"{missing_package} must be installed before running setup.py")

# Enforce Python 3 requirement
if sys.version_info < (3, 6):
    raise SystemExit("Python 3.6 or later is required to install this package.")

# Minimum setuptools version requirement
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    raise ImportError("Setuptools version >=1.1 is required")

# External dependency check
external_dependencies = [("wigToBigWig", "UCSC tools", 255), ("bedtools", "bedtools", 0)]
for tool, toolkit, expected in external_dependencies:
    retcode = subprocess.call(tool, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if retcode != expected:
        print(f"WARNING: Dependency check for {toolkit} ({tool}) failed with error code {retcode}")

# Adjust packages and directories
cgat_packages = find_packages(include=["cgat", "cgat.*"], exclude=['tests'])
cgat_package_dirs = {'cgat': 'cgat'}

# Build setup configuration
setup(
    name='cgat',
    version='0.7.4.1',
    description='cgat : the Computational Genomics Analysis Toolkit',
    author='Andreas Heger',
    author_email='andreas.heger@gmail.com',
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description='cgat : the Computational Genomics Analysis Toolkit',
    python_requires=">=3.6",
    packages=cgat_packages,
    package_dir=cgat_package_dirs,
    include_package_data=True,
    entry_points={'console_scripts': ['cgat = cgat.cgat:main']},
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
    test_suite="tests",
)
