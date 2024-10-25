
"""test_scripts.py
==================

test script for CGAT scripts

This script builds test cases from directories in the :file:`tests`
subdirectories. Test data for scripts are contained in a directory
with the name of the script and described in a :file:`tests.yaml`
within that directory.

To permit the parallelization of running tests, tests can be run in
chunks (tasks). The script will look for the following environment variables:

CGAT_TASK_ID
   The starting number of this task starting at 1

CGAT_TASK_STEPSIZE
   The number of tests to run within a chunk

"""
# tests/test_scripts.py

import subprocess
import tempfile
import os
import shutil
import re
import glob
import gzip
import yaml
import time
import hashlib
import platform
import pytest

import TestUtils  # Ensure this module is compatible with pytest

TRAVIS = os.environ.get("TRAVIS", False) == "true"
JENKINS = os.environ.get("USER", "") == "jenkins"
PYTHON_VERSION = platform.python_version()
SUBDIRS = ("gpipe", "optic")

# Directory where tools are located, relative to root
TOOLS_DIR = os.path.join("cgat", "tools")

# Setup logging
LOGFILE = open("test_scripts.log", "a")
DEBUG = os.environ.get("CGAT_DEBUG", False)


def compute_checksum(filename):
    '''Return MD5 checksum of file.'''
    with open(filename, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()


def _read(fn):
    if fn.endswith(".gz"):
        with gzip.open(fn) as inf:
            data = inf.read()
    else:
        with open(fn, "rb") as inf:
            data = inf.read()

    try:
        data = data.decode("ascii")
    except UnicodeDecodeError:
        return data

    data = [x for x in data.splitlines() if not x.startswith("#")]

    return data


def _check_script(test_name, script, stdin,
                  options, outputs,
                  references,
                  working_dir,
                  current_dir):
    '''Check script.'''
    tmpdir = tempfile.mkdtemp()

    t1 = time.time()

    stdout = os.path.join(tmpdir, 'stdout')
    if stdin:
        if stdin.endswith(".gz"):
            # zcat on osX requires .Z suffix
            stdin_cmd = f'gunzip < {os.path.abspath(working_dir)}/{stdin} |'
        else:
            stdin_cmd = f'cat {os.path.abspath(working_dir)}/{stdin} |'
    else:
        stdin_cmd = ""

    if options:
        options = re.sub(r"%TMP%|<TMP>", tmpdir, options)
        options = re.sub(r"%DIR%|<DIR>", os.path.abspath(working_dir), options)
    else:
        options = ""

    options = options.replace("\n", "")
    short_name = os.path.basename(script)[:-3]

    # Use /bin/bash to enable "<( )" syntax in shells
    statement = (f"/bin/bash -c '{stdin_cmd} cgat {short_name} {options} > {stdout}'")

    process = subprocess.Popen(statement,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               cwd=tmpdir)

    if DEBUG:
        print(f"tmpdir={tmpdir}", end=" ")

    process_stdout, process_stderr = process.communicate()

    fail = False
    msg = ""

    if process.returncode != 0:
        fail = True
        msg = f"error in statement: {statement}; stderr={process_stderr.decode()}"

    # For version tests, do not compare output
    if test_name == "version":
        pass
    elif not fail:
        # Compare line by line, ignoring comments
        for output, reference in zip(outputs, references):
            if output == "stdout":
                output_path = stdout
            elif output.startswith("<DIR>/") or output.startswith("%DIR%/"):
                output_path = os.path.join(working_dir, output[6:])
            else:
                output_path = os.path.join(tmpdir, output)

            if not os.path.exists(output_path):
                fail = True
                msg = f"output file '{output_path}' does not exist: {statement}"
                break

            reference_path = os.path.join(working_dir, reference)
            if not os.path.exists(reference_path):
                fail = True
                msg = f"reference file '{reference_path}' does not exist ({tmpdir}): {statement}"
                break

            a = _read(output_path)
            b = _read(reference_path)
            if a != b:
                fail = True
                msg = (f"files {output_path} and {reference_path} are not the same\n"
                       f"{statement}\nmd5: output={compute_checksum(output_path)}, "
                       f"reference={compute_checksum(reference_path)}")

                diffs = []
                for aa, bb in zip(a, b):
                    if aa != bb:
                        diffs.append((aa, bb))
                        if len(diffs) > 10:
                            break

                diff_str = "\n--\n".join(["\n".join(map(str, x)) for x in diffs])
                msg += f"first 10 differences: {diff_str}"
                break

    t2 = time.time()
    LOGFILE.write(f"{script}\t{test_name}\t{t2 - t1}\n")
    LOGFILE.flush()

    # Preserve coverage information
    coverage_files = glob.glob(os.path.join(tmpdir, ".coverage*"))
    for f in coverage_files:
        shutil.move(os.path.abspath(f),
                    os.path.join(current_dir, os.path.basename(f)))

    if not DEBUG:
        shutil.rmtree(tmpdir)

    assert not fail, msg


def check_script(test_name, script, stdin,
                 options, outputs,
                 references,
                 working_dir,
                 current_dir):
    '''Check script and assert no failure.'''
    _check_script(test_name, script, stdin,
                  options, outputs,
                  references,
                  working_dir,
                  current_dir)


def check_main(script):
    '''Test if a script can be imported and has a main function.'''
    # Substitute gpipe and other subdirectories.
    for s in SUBDIRS:
        script = re.sub(f"{s}_", f"{s}/", script)

    # Check for text match
    with open(script) as inf:
        code = inf.read()
    assert "def main(" in code or "import main" in code, "No main function"


def collect_test_cases():
    '''Collect all test cases for parametrization.'''
    current_dir = os.getcwd()
    testing_dir = TestUtils.get_tests_directory()
    scripts_dir = os.path.join(os.path.dirname(testing_dir), TOOLS_DIR)
    test_dirs = glob.glob(os.path.join(testing_dir, "*.py"))

    # The config file
    config_file = os.path.join(testing_dir, "_test_scripts.yml")

    if os.path.exists(config_file):
        with open(config_file) as f:
            config = yaml.safe_load(f)
        if config and "restrict" in config and config["restrict"]:
            values = config["restrict"]
            if "glob" in values:
                test_dirs = glob.glob(os.path.join(testing_dir, values["glob"]))
            if "manifest" in values:
                # Take scripts defined in the MANIFEST.in file
                with open("MANIFEST.in") as f:
                    manifest_lines = f.readlines()
                test_dirs = [re.sub(r"include\s*cgat/tools/", "tests/",
                                    x.strip()) for x in manifest_lines
                             if x.startswith("include cgat/tools") and x.endswith(".py")]
            if "regex" in values:
                rx = re.compile(values["regex"])
                test_dirs = list(filter(rx.search, test_dirs))

    # Ignore those which don't exist as tests
    test_dirs = [x for x in test_dirs if os.path.exists(x)]

    # Ignore non-directories
    test_dirs = [x for x in test_dirs if os.path.isdir(x)]

    test_dirs.sort()

    # Restrict tests run according to chunk parameters
    starting_test_number = os.getenv('CGAT_TASK_ID', None)
    test_increment = os.getenv('CGAT_TASK_STEPSIZE', None)

    try:
        starting_test_number, test_increment = \
            (int(starting_test_number) - 1,
             int(test_increment))
        test_dirs = test_dirs[starting_test_number:
                              starting_test_number + test_increment]
    except (TypeError, ValueError):
        pass

    test_cases = []

    for test_script in test_dirs:
        script_name = os.path.basename(test_script)
        script_path = os.path.abspath(os.path.join(scripts_dir, script_name))

        # Add main function check
        test_cases.append({
            "func": check_main,  # Direct function reference
            "args": (script_path,),
            "description": os.path.join(script_name, "def_main")
        })

        tests_yaml = os.path.join(test_script, "tests.yaml")
        if not os.path.exists(tests_yaml):
            continue

        with open(tests_yaml) as f:
            script_tests = yaml.safe_load(f)

        for test, values in sorted(script_tests.items()):
            if "skip_python" in values:
                versions = [x.strip() for x in str(values["skip_python"]).split(",")]
                if any(PYTHON_VERSION.startswith(x) for x in versions):
                    continue
            if "skip_travis" in values and TRAVIS:
                continue
            if "skip_jenkins" in values and JENKINS:
                continue

            # Handle scripts in subdirectories
            if "_" in script_name:
                parts = script_name.split("_")
                potential_path = os.path.join("tools", parts[0], "_".join(parts[1:]))
                if os.path.exists(potential_path):
                    script_path = os.path.abspath(os.path.join(scripts_dir, script_name))
                else:
                    script_path = os.path.abspath(os.path.join(scripts_dir, script_name))

            test_cases.append({
                "func": check_script,  # Direct function reference
                "args": (
                    test,
                    script_path,
                    values.get('stdin', None),
                    values['options'],
                    values['outputs'],
                    values['references'],
                    test_script,
                    current_dir
                ),
                "description": os.path.join(script_name, test)
            })

    return test_cases


# Collect all test cases
TEST_CASES = collect_test_cases()


@pytest.mark.parametrize("test_case", TEST_CASES, ids=[tc["description"] for tc in TEST_CASES])
def test_scripts(test_case):
    '''Run parametrized script tests.'''
    func = test_case["func"]  # Directly retrieve the function object
    args = test_case["args"]
    func(*args)
