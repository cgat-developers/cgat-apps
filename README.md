
<p align="left">
<a href='https://cgat-apps.readthedocs.io/en/latest/?badge=latest'>
  <img src='https://readthedocs.org/projects/cgat-apps/badge/?version=latest' alt='Documentation Status' /></a>
  <a href='https://github.com/cgat-developers/cgat-apps/actions/workflows/cgatapps_python.yml'>
    <img src='https://github.com/cgat-developers/cgat-apps/actions/workflows/cgatapps_python.yml/badge.svg' alt='Github Action status' /></a>
  <a href="https://twitter.com/cgat_oxford?lang=en", alt="Twitter followers">
    <img src="https://img.shields.io/twitter/url/http/shields.io.svg?style=social&logo=twitter" /></a>
  <a href="https://twitter.com/cgat_oxford?lang=en", alt="Twitter followers">
	<img src="https://img.shields.io/twitter/url/http/shields.io.svg?style=social&logo=twitter" /></a>
</p>

CGAT Apps
=========

CGAT Apps is a collection of scripts to aid the analysis of high-throughput sequencing data.

After installation, use the ``cgat`` command to see how to use them.

Documentation
=============

We are attempting to improve our documentation. However, our current documentation
can be found [here](https://cgat-apps.readthedocs.io/en/latest/)

For questions, please open a discussion on the GitHub
[issue page](https://github.com/cgat-developers/cgat-apps/issues)

Installation
============

Conda Installation
------------------
The preferred method to install CGAT Apps is using the installation script, which uses
[mamba](https://github.com/mamba-org/mamba), the fast C implementation of [Conda](https://conda.io).

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

Usage
=====

Run the ``cgat --help`` command to see what scripts are available and how to use them.
For example, to strip sequence and quality information from a bam_ file, type::

   cgat bam2bam --strip=sequence < in.bam > out.bam
