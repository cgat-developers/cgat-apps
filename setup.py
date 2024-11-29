# setup.py
import sysconfig
import sys
import os
import subprocess
import re

# Import setuptools at the beginning
import setuptools
from setuptools import setup, find_packages, Extension
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

# Cython extensions and paths
conda_includes = [os.path.dirname(sysconfig.get_paths()["include"])]
conda_libdirs = [os.path.dirname(sysconfig.get_paths()["stdlib"])]
pysam_libraries = pysam.get_libraries()
pysam_libdirs = list(set(os.path.dirname(x) for x in pysam_libraries)) + conda_libdirs
pysam_libs = ["hts"] + [os.path.basename(x)[3:-3] for x in pysam_libraries]
pysam_dirname = os.path.dirname(pysam.__file__)
extra_link_args_pysam = [f'-Wl,-rpath,{x}' for x in pysam_libdirs + conda_libdirs]

extensions = [
    Extension(
        'cgat.Components',
        ['cgat/Components/Components.pyx', 'cgat/Components/connected_components.cpp'],
        include_dirs=[os.path.join('cgat', 'Components')] + conda_includes,
        language="c++",
    ),
    Extension(
        "cgat.NCL.cnestedlist",
        ["cgat/NCL/cnestedlist.pyx", "cgat/NCL/intervaldb.c"],
        language="c",
    ),
    Extension(
        "cgat.GeneModelAnalysis",
        ["cgat/GeneModelAnalysis.pyx"],
        include_dirs=conda_includes + pysam.get_include() + [numpy.get_include()],
        define_macros=pysam.get_defines(),
        language="c",
    ),
    Extension(
        "cgat.BamTools.bamtools",
        ["cgat/BamTools/bamtools.pyx"],
        include_dirs=conda_includes + pysam.get_include() + [numpy.get_include()],
        library_dirs=pysam_libdirs,
        libraries=pysam_libs,
        define_macros=pysam.get_defines(),
        language="c",
        extra_link_args=extra_link_args_pysam,
    ),
    Extension(
        "cgat.BamTools.geneprofile",
        ["cgat/BamTools/geneprofile.pyx"],
        include_dirs=conda_includes + pysam.get_include() + [numpy.get_include()],
        library_dirs=pysam_libdirs,
        libraries=pysam_libs,
        define_macros=pysam.get_defines(),
        language="c",
        extra_link_args=extra_link_args_pysam,
    ),
    Extension(
        "cgat.BamTools.peakshape",
        ["cgat/BamTools/peakshape.pyx"],
        include_dirs=conda_includes + pysam.get_include() + [numpy.get_include()],
        library_dirs=pysam_libdirs,
        libraries=pysam_libs,
        define_macros=pysam.get_defines(),
        language="c",
        extra_link_args=extra_link_args_pysam,
    ),
    Extension(
        "cgat.VCFTools",
        ["cgat/VCFTools/vcftools.pyx"],
        include_dirs=conda_includes + pysam.get_include() + [numpy.get_include()],
        library_dirs=pysam_libdirs,
        libraries=pysam_libs,
        define_macros=pysam.get_defines(),
        language="c",
        extra_link_args=extra_link_args_pysam,
    ),
    Extension(
        "cgat.FastqTools",
        ["cgat/FastqTools/fastqtools.pyx"],
        include_dirs=conda_includes + pysam.get_include() + [numpy.get_include()],
        library_dirs=pysam_libdirs,
        libraries=pysam_libs,
        define_macros=pysam.get_defines(),
        language="c",
        extra_link_args=extra_link_args_pysam,
    ),
]

# Build setup configuration
setup(
    packages=cgat_packages,
    package_dir=cgat_package_dirs,
    include_package_data=True,
    ext_modules=extensions,
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'cgat = cgat.cgat.main',
            ],
        },
    test_suite="tests",
)
