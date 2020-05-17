.. image:: https://travis-ci.org/cgat-developers/cgat-apps.svg?branch=master
    :target: https://travis-ci.org/cgat-developers/cgat-apps

=========
CGAT Apps
=========

CGAT Apps is a collection of scripts to analyse high-throughput sequencing data.

After installation, use the ``cgat`` command to see how to use them.

We are attempting to improve our documentation. However, our current documentation
can be found `here <https://www.cgat.org/downloads/public/cgat/documentation/>`_

For questions, please open a discussion on the GitHub 
`issue page <https://github.com/cgat-developers/cgat-apps/issues>`_.

Installation
============

End users: install using Conda
------------------------------
The preferred method to install CGAT Apps is using the installation script, which uses
`Conda <https://conda.io>`_.

To install cgat-apps using conda::
    
    conda install -c conda-forge -c bioconda cgat-apps

Developers: try the installation script
---------------------------------------
Alternatively, installation can be performed using our install script that will
automate the process for you.

Here are the steps::

        # download installation script:
        curl -O https://raw.githubusercontent.com/cgat-developers/cgat-apps/master/install.sh

        # see help:
        bash install.sh

        # install the development version (recommended, no production version yet):
        bash install.sh --devel [--location </full/path/to/folder/without/trailing/slash>]

        # the code is downloaded in zip format by default. If you want to get a git clone, use:
        --git # for an HTTPS clone
        --git-ssh # for a SSH clone (you need to be a cgat-developer contributor on GitHub to do this)

        # enable the conda environment as requested by the installation script
        # NB: you probably want to automate this by adding the instructions below to your .bashrc
        source </full/path/to/folder/without/trailing/slash>/conda-install/etc/profile.d/conda.sh
        conda activate base
        conda activate cgat-a

        # finally, please run the cgatflow command-line tool to check the installation:
        cgat --help

The installation script will put everything under the specified location. The aim is to provide a portable
installation that does not interfere with the existing software. As a result, you will get a conda environment
working with CGAT Apps which can be enabled on demand according to your needs.

Usage
=====

Run the ``cgat --help`` command to see what scripts are available and how to use them.
For example, to strip sequence and quality information from a bam_ file, type::

   cgat bam2bam --strip=sequence < in.bam > out.bam

.. _cgat-core: https://github.com/cgat-developers/cgat-core
.. _bam: http://en.wikipedia.org/wiki/SAMtools
