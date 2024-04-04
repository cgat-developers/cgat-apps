'''
test_commandline - Tests coding style conformity of CGAT code collection.
==========================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------
This script tests the command line usage of all scripts in the CGAT code collection.
It's recommended to run this within nosetests for efficiency and coverage:

   nosetests tests/test_commandline.py --nocapture

Before running these tests, ensure to execute:

   python setup.py develop

to make all package scripts available for import and testing.
'''

import glob
import os
import importlib
import yaml
import re
import sys
import copy
import argparse

from nose.tools import ok_
import cgatcore.experiment as E
import cgatcore.iotools as iotools
import TestUtils

# Preserve the original E.Start function for later restoration
ORIGINAL_START = E.start

# Placeholders for parser object and tested script's module
PARSER = None
TESTED_MODULE = None

# Directories to examine for Python modules/scripts
EXPRESSIONS = (
    ('tools', 'cgat/tools/*.py'),
)

# Files to exclude from checks
EXCLUDE = [
    "__init__.py",
    "version.py",
    "cgat.py",
    "gtf2table.py",   # Fails with pysam include issue
    "bed2table.py",   # Fails with pysam include issue
    "fasta2bed.py",   # Fails due to pybedtools rebuild requirements
]

# Filename for the black/white list of options
FILENAME_OPTIONLIST = "tests/option_list.tsv"


class DummyError(Exception):
    """Custom exception for controlling test flow."""
    pass


def filter_files(files):
    '''Filter list of files according to filters set in the configuration file tests/_test_commandline.yml'''
    testing_dir = TestUtils.get_tests_directory()
    config_file = os.path.join(testing_dir, "_test_commandline.yml")

    if os.path.exists(config_file):
        with open(config_file) as cf:
            config = yaml.safe_load(cf)
            if config and "restrict" in config:
                values = config["restrict"]
                if "manifest" in values:
                    scriptdirs = [x.strip() for x in open("MANIFEST.in")
                                  if x.startswith("include CGAT/tools") and x.endswith(".py\n")]
                    take = set(re.sub("include\s*", "", x) for x in scriptdirs)
                    files = [x for x in files if x in take]
                if "regex" in values:
                    rx = re.compile(values["regex"])
                    files = list(filter(rx.search, files))
    return files


def LocalStart(parser, *args, **kwargs):
    '''Stub for E.start - captures the parser for inspection.'''
    global PARSER
    kwargs.update({'return_parser': True})
    PARSER = ORIGINAL_START(parser, **kwargs)
    raise DummyError()


def load_script(script_name):
    '''Attempts to import a script as a module for inspection.'''
    script_path = os.path.splitext(script_name)[0]
    script_dir, script_base = os.path.split(script_path)
    module_name = ".".join(filter(None, [script_dir.replace(os.sep, '.'), script_base]))

    # Remove compiled files to ensure fresh import
    compiled_script = script_path + ".pyc"
    if os.path.exists(compiled_script):
        os.remove(compiled_script)

    try:
        module = importlib.import_module(module_name)
    except ImportError as e:
        sys.stderr.write(f'ImportError for {module_name}: {e}\n')
        return None, None

    return module, module_name


def test_cmdline():
    '''Test command line interfaces of scripts for style and conformity.'''
    global ORIGINAL_START

    # Load option actions from list
    option_actions = iotools.read_map(
        iotools.open_file(FILENAME_OPTIONLIST),
        columns=(0, 1),
        has_header=True
    )

    # Compile list of scripts to test
    files = [f for label, expr in EXPRESSIONS for f in glob.glob(expr)]
    files = filter_files(files)

    # Prioritise the current directory for module lookup
    sys.path.insert(0, ".")

    for script in files:
        if os.path.isdir(script) or os.path.basename(script) in EXCLUDE:
            continue

        script_name = os.path.abspath(script)
        module, module_name = load_script(script)
        if not module:
            yield fail_, f"Module {script_name} could not be imported."
            continue

        # Replace the start function to capture parser
        E.start = LocalStart

        # Attempt to run script's main function to access its parser
        try:
            module.main(argv=["--help"])
        except DummyError:
            # Expected flow interruption by LocalStart
            pass
        except Exception as e:
            yield fail_, f"Error invoking main of {script_name}: {e}"
            continue

        if PARSER:
            for action in PARSER._actions:  # Iterate through the actions stored in the parser
                if isinstance(action, argparse._HelpAction):  # Skip help actions
                    continue
                opt_strings = action.option_strings  # Get the list of CLI flags
                if not opt_strings:  # This skips positional arguments
                    continue
                for opt_string in opt_strings:
                    if opt_string.startswith("--"):
                        opt_string = opt_string[2:]
                    yield check_option, opt_string, script_name, option_actions

        # Reset module to avoid conflicts
        if module_name in sys.modules:
            del sys.modules[module_name]


def check_option(option, script_name, option_actions):
    print(f"Checking option: {option} in script: {script_name}")  # Diagnostic print
    assert option in option_actions, f"Option {option} in script {script_name} is unknown or not allowed."
    assert option_actions[option] == "ok", f"Option {option} in script {script_name} is not allowed."


def fail_(msg):
    '''Generate a failing test with the provided message.'''
    ok_(False, msg)

# Reset E.start to its original function after testing
E.start = ORIGINAL_START
