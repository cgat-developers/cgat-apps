'''table2table.py - operate on tables
==================================

:Tags: Python

Purpose
-------

This script implements a few methods for manipulating tables.

Methods working on all tables:
++++++++++++++++++++++++++++++

transpose
   transpose a table

split-fields
   Split muliple-value fields in each row at ``--separator``. Output
   multiple rows with all combinations.

group
   Group values by column

join-column
   Join rows in a table by columns

expand-table
   If a field in a row contains multiple values,
   the row is expanded into multiple rows such
   that all values have space.

flatten-table
   Output a table as row/column/value tuples.

as-column
   Output table as a single column. Colums in the original table are
   appended and output.

collapse-table
   Collapse a table of two columns with row names in the first
   column. Outputs a table with multiple columns for each row name.

Methods for numerical columns
+++++++++++++++++++++++++++++

Some methods make only sense for columns containing numerical values.
If a table contains both numerical and non-numerical data, the
numerical columns can be specified by the ``--columns`` option.

normalize-by-value
   divide all cells in a table by a value

multiply-by-value
   multiply all cells in a table by a value

lower-bound
   replace all cells with a value of less than lower bound with the lower
   bound.

upper-bound
   replace all cells with a value of more than upper bound with the upper
   bound.

normalize-by-table
   divide each cell in a table with the corresponding entry in a secondary
   table.

normalize-by-max
   divide table columns by maximum per column

kullback-leibler
   compute kullback-leibler divergence between two columns. Compute
   both D(a||b), D(b||a) and (D(a||b) + D(b||a)) / 2

rank
   substitute cells with their ranks in a column

fdr
   compute an FDR over selected columns. Replaces the columns
   with the qvalues.

Usage
-----

Example::

   python table2table.py --help

Type::

   python table2table.py --help

for command line help.

Command line options
--------------------

'''
import sys
import cgatcore.table
from cgatcore.table import main
from cgat import Stats

# monkeypatch: import Stats in cgatcore.Table's namespace.
cgatcore.table.Stats = Stats


if __name__ == "__main__":
    sys.exit(main())
