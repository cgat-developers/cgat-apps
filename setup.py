import sysconfig
import sys
import os
import subprocess
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from Cython.Distutils import build_ext as cython_build_ext
import numpy
import pysam

# Enforce Python version
if sys.version_info < (3, 6):
    raise SystemExit("Python 3.6 or later is required to install this package.")

# Package and directory settings
cgat_packages = find_packages(include=["cgat", "cgat.*"], exclude=["tests"])
cgat_package_dirs = {"cgat": "cgat"}

# Paths for build dependencies
conda_includes = [os.path.dirname(sysconfig.get_paths()["include"])]
conda_libdirs = [os.path.dirname(sysconfig.get_paths()["stdlib"])]
pysam_libdirs = list(set(os.path.dirname(lib) for lib in pysam.get_libraries())) + conda_libdirs
pysam_libs = ["hts"] + [os.path.basename(lib)[3:-3] for lib in pysam.get_libraries()]
extra_link_args_pysam = [f"-Wl,-rpath,{path}" for path in pysam_libdirs + conda_libdirs]

# Define Cython extensions
extensions = [
    Extension(
        "cgat.Components",
        ["cgat/Components/Components.pyx", "cgat/Components/connected_components.cpp"],
        include_dirs=["cgat/Components"] + conda_includes,
        language="c++",
    ),
    Extension(
        "cgat.NCL.cnestedlist",
        ["cgat/NCL/cnestedlist.pyx", "cgat/NCL/intervaldb.c"],
        include_dirs=["cgat/NCL"] + conda_includes,
        language="c",
    ),
    Extension(
        "cgat.GeneModelAnalysis",
        ["cgat/GeneModelAnalysis.pyx"],
        include_dirs=conda_includes + pysam.get_include() + [numpy.get_include()],
        define_macros=pysam.get_defines(),
        language="c",
    ),
    # Additional extensions for BamTools, VCFTools, and FastqTools
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
        "cgat.VCFTools",
        ["cgat/VCFTools/vcftools.pyx"],
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
    name="cgat",
    version="0.7.7",
    description="Computational Genomics Analysis Toolkit",
    author="Adam Cribbs",
    author_email="adam.cribbs@ndorms.ox.ac.uk",
    license="MIT",
    packages=cgat_packages,
    package_dir=cgat_package_dirs,
    include_package_data=True,
    package_data={
        "cgat.Components": ["*.h"],
    },
    ext_modules=extensions,
    cmdclass={"build_ext": cython_build_ext},
    entry_points={
        "console_scripts": [
            "cgat = cgat.cgat:main",
        ]
    },
    zip_safe=False,
    install_requires=[
        "Cython>=0.29.35",
        "numpy",
        "pysam",
    ],
)
