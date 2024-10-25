'''test_import - test importing all modules and pipelines
=========================================================

:Author: Adam Cribbs
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script attempts to import all the python libraries and
pipeline scripts in the CGAT code collection.

Importing a script/module is a pre-requisite for building
documentation with sphinx. A script/module that can not be imported
will fail within sphinx.

This script is best run within nosetests::

   nosetests tests/test_import.py


'''
import os
import glob
import traceback
import importlib.util

# DIRECTORIES to examine for python modules/scripts
EXPRESSIONS = (
    ('tests', 'tests/*.py'),
    ('tools', 'cgat/tools/*.py'),
    ('cgat', 'cgat/*.py'))

# Scripts to exclude as they fail imports.
EXCLUDE = (
    'pipeline_intervals', 'PipelinePeakcalling', 'IndexedFasta',
    'pipeline_peakcalling', 'bam2transcriptContribution', 'beds2counts',
    'fasta2bed', 'bed2table', 'table2table', 'combine_tables')


def check_import(filename, outfile):
    prefix, suffix = os.path.splitext(filename)
    dirname, basename = os.path.split(prefix)

    if basename in EXCLUDE:
        return

    if os.path.exists(prefix + ".pyc"):
        try:
            os.remove(prefix + ".pyc")
        except OSError:
            pass

    # ignore scripts with pyximport for now
    blob = open(filename).read()
    if "import pyximport" in blob:
        return

    try:
        # Use importlib to load the module dynamically
        spec = importlib.util.spec_from_file_location(basename, filename)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

    except ImportError as msg:
        outfile.write(f"FAIL {basename}\n{msg}\n")
        outfile.flush()
        traceback.print_exc(file=outfile)
        assert False, f'{basename} scripts/modules - ImportError: {msg}'
    except Exception as msg:
        outfile.write(f"FAIL {basename}\n{msg}\n")
        outfile.flush()
        traceback.print_exc(file=outfile)
        assert False, f'{basename} scripts/modules - Exception: {msg}'

    assert True


def test_imports():
    '''test importing modules and scripts'''

    with open('test_import.log', 'a') as outfile:
        for label, expression in EXPRESSIONS:
            files = glob.glob(expression)
            files.sort()

            for f in files:
                if os.path.isdir(f):
                    continue
                check_import(os.path.abspath(f), outfile)
