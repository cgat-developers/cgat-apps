import sysconfig
import sys
import os
import subprocess
import re
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

# Minimum setuptools version requirement
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    raise ImportError("Setuptools version >=1.1 is required")

# Define version and other package information
sys.path.insert(0, "cgat")
import version
version = version.__version__

IS_OSX = sys.platform == 'darwin'

# External dependency check
external_dependencies = [("wigToBigWig", "UCSC tools", 255), ("bedtools", "bedtools", 0)]
for tool, toolkit, expected in external_dependencies:
    retcode = subprocess.call(tool, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if retcode != expected:
        print(f"WARNING: Dependency check for {toolkit} ({tool}) failed with error code {retcode}")

# Python version check
if sys.version_info < (2, 7):
    raise SystemExit("Python 2.7 or later is required")

# Dependency and package configuration
install_requires = []
dependency_links = []
with open('requires.txt') as req_file:
    for line in req_file:
        if line.startswith("https"):
            install_requires.append(line.strip())
        else:
            REPO_REQUIREMENT = re.compile(r'^-e (?P<link>(?P<vcs>git|svn|hg|bzr).+#egg=(?P<package>.+)-(?P<version>\d(?:\.\d)*))$')
            HTTPS_REQUIREMENT = re.compile(r'^-e (?P<link>.*).+#(?P<package>.+)-(?P<version>\d(?:\.\d)*)$')
            match = REPO_REQUIREMENT.match(line)
            if match:
                install_requires.append(f"{match.group('package')}=={match.group('version')}")
                dependency_links.append(match.group('link'))
            else:
                match = HTTPS_REQUIREMENT.match(line)
                if match:
                    install_requires.append(f"{match.group('package')}=={match.group('version')}")
                    dependency_links.append(match.group('link'))
                else:
                    install_requires.append(line.strip())

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
    name='cgat',
    version=version,
    description='cgat : the Computational Genomics Analysis Toolkit',
    author='Andreas Heger',
    author_email='andreas.heger@gmail.com',
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description='cgat : the Computational Genomics Analysis Toolkit',
    classifiers=[_f for _f in """
        Development Status :: 3 - Alpha
        Intended Audience :: Science/Research
        Intended Audience :: Developers
        License :: OSI Approved
        Programming Language :: Python
        Topic :: Software Development
        Topic :: Scientific/Engineering
        Operating System :: POSIX
        Operating System :: Unix
        Operating System :: MacOS
    """.splitlines() if _f],
    url="http://www.cgat.org/cgat/Tools/",
    packages=cgat_packages,
    package_dir=cgat_package_dirs,
    include_package_data=True,
    entry_points={'console_scripts': ['cgat = cgat.cgat:main']},
    install_requires=install_requires,
    dependency_links=dependency_links,
    ext_modules=extensions,
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
    test_suite="tests",
)
