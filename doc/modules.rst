.. _modules:

=======
Modules
=======

This section documents the modules used in CGAT scripts.

CGAT generic toolboxes
======================

These are the modules that every script or module should use.

.. toctree::
   :maxdepth: 1

   modules/experiment.rst
   modules/iotools.rst

Genomics
========

File formats
------------

Modules for parsing and working for data in specific formats.

.. toctree::
   :maxdepth: 1

   modules/AGP.rst
   modules/BamTools.rst
   modules/Bed.rst
   modules/Biomart.rst
   modules/Blat.rst
   modules/CBioPortal.rst
   modules/FastaIterator.rst
   modules/Fastq.rst
   modules/GFF3.rst
   modules/GTF.rst
   modules/IndexedFasta.rst
   modules/IndexedGenome.rst
   modules/Sra.rst
   modules/VCF.rst

Algorithms
----------

.. toctree::
   :maxdepth: 1

   modules/AString.rst
   modules/Counts.rst
   modules/Expression.rst
   modules/Genomics.rst
   modules/Intervals.rst
   modules/Motifs.rst
   modules/SequenceProperties.rst
   modules/Variants.rst


Wrappers
--------

These modules wrap tools and provide routines for parsing their
output.

.. toctree::
   :maxdepth: 1

   modules/WrapperCodeML.rst
   modules/Bioprospector.rst
   modules/MAST.rst
   modules/IGV.rst
   modules/Masker.rst

Data processing
===============

Math and Stats
--------------

.. toctree::
   :maxdepth: 1

   modules/Histogram.rst
   modules/Histogram2D.rst
   modules/Stats.rst
   modules/MatrixTools.rst

Toolboxes
----------

Toolboxes for generic problems.

.. toctree::
   :maxdepth: 1

   modules/CSV.rst
   modules/Iterators.rst
   modules/Database.rst
   modules/SetTools.rst
   modules/Tree.rst
   modules/TreeTools.rst


CGAT infrastructure
===================

Below is a list of modules that are involved in maintainig the
CGAT infrastructure such as logging, dependency tracking, etc.

.. toctree::
   :maxdepth: 1

   modules/Requirements.rst
   modules/Logfile.rst
   modules/CSV2DB.rst

Other
-----

.. toctree::
   :maxdepth: 1

   modules/RLE.rst
   modules/SVGdraw.rst
   modules/RateEstimation.rst

Unsorted
---------

Modules not sorted into categories.

.. toctree::
   :maxdepth: 2
