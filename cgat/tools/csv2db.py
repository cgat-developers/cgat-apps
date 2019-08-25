'''csv2db.py - upload table to database
====================================

:Tags: Python

Purpose
-------

create a table from a csv separated file and load data into it.

This module supports backends for postgres and sqlite3. Column types are
auto-detected.

Read a table from stdin and create an sqlite3 database. By default,
the database will reside in a file called csvdb and in a table csv.

.. todo::

   Use file import where appropriate to speed up loading. Currently, this is
   not always the case.

Usage
-----

Example::

   python csv2db.py -b sqlite < stdin 

Type::

   python csv2db.py --help

for command line help.

Command line options
--------------------

'''

import sys
from cgatcore.csv2db import main

import csv
csv.field_size_limit(sys.maxsize)


if __name__ == "__main__":
    sys.exit(main())
