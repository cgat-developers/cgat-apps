.. _CGATInstallation:

=========================
Installation instructions
=========================

The section below describes how to install the cgat-apps. We distinguish between two different installation
types: production and development. The former refers to the released version of our tools, which can be installed
using pip or conda, and is the recommended installation for users. The latter refers to the installation
of the development version of the apps, which can be used to make changes to our code base.

Please note that we can not test our code on all systems and configurations out there so please bear with us.

Quick installation
==================

Install using Conda
-------------------

Conda Installation
------------------
The preferred method to install CGAT Apps is using the installation script, which uses
mamba_, the fast C implementation of conda_.

To install cgat-apps using mamba::

    mamba install -c conda-forge -c bioconda cgat-apps

Developers: try the installation script
---------------------------------------

Here are the steps::

        # Install conda and mamba according the the documentation. Next
        # install the conda packages for cgat-apps to work
        conda env create -f conda/environments/cgat-apps.yml

        # enable the conda environment
        conda activate cgat-a

        # Install the development version of cgat-apps
        python setup.py develop

        # finally, please run the cgat command-line tool to check the installation:
        cgat --help

The installation script will put everything under the specified location. The aim is to provide a portable
installation that does not interfere with the existing software. As a result, you will get a conda environment
working with CGAT Apps which can be enabled on demand according to your needs.


Install using pip
-----------------

You can also use pip_ to install the CGAT scripts. To go down this route, please type::

   pip install cgat

However, cgat-apps depends on numerous other python packages which themselves might require
manual intervention.


.. _conda: https://conda.io
.. _pip: https://pypi.python.org/pypi/CGAT
.. _here: https://doi.org/10.1093/bioinformatics/btt756
.. _mamba: https://github.com/mamba-org/mamba
