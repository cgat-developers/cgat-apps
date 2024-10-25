'''test_style - test coding style confirmation of CGAT code
===========================================================

:Author: Adam Cribbs
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script runs pep8 on all scripts in the CGAT
code collection.

This script is best run within nosetests::

   nosetests tests/test_style.py

'''
import pep8
import glob
import os

# DIRECTORIES to examine for python modules/scripts
EXPRESSIONS = (
    ('tests', 'tests/*.py'),
    ('scripts', 'cgat/tools/*.py'),
    ('cgat', 'cgat/*.py')
)

# Codes to ignore in the pep8 BaseReport
IGNORE = set(('E101', 'E201', 'E202', 'E122', 'E265', 'E501', 'E502', 'E731',
              'W191', 'W291', 'W293', 'W391', 'W503', 'W601', 'W602', 'E701', 'E713', 'E702'))


def check_style(filename):
    '''check style of filename.'''

    p = pep8.StyleGuide(quiet=True)
    report = p.check_files([filename])

    # Exclude non-violation counters from the result
    take = [y for x, y in report.counters.items() if x not in IGNORE and x not in ('files', 'directories', 'physical lines', 'logical lines')]
    found = ['%s:%i' % (x, y) for x, y in report.counters.items() if x not in IGNORE and x not in ('files', 'directories', 'physical lines', 'logical lines')]
    total = sum(take)

    assert total == 0, f'pep8 style violations in {filename}: {", ".join(found)}'


def test_style():
    '''test style of scripts'''

    for label, expression in EXPRESSIONS:
        files = glob.glob(expression)
        files.sort()

        for f in files:
            if os.path.isdir(f):
                continue
            check_style(os.path.abspath(f))
