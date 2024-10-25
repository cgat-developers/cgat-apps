'''
cgat_get_options.py - build a sorted list of all options used in scripts
========================================================================

:Author:
:Tags: Python

Purpose
-------

Go through all scripts in the cgat code collection and collect
options used in the scripts.

This script expects to be executed at the root of the
cgat code repository.


Usage
-----

.. Example use case

Example::

   python cgat_get_options.py

Type::

   python cgat_get_options.py --help

for command line help.

Command line options
--------------------

'''
import sys
import os
import glob
import importlib.util  # Use importlib instead of imp
import collections
import pandas
import cgatcore.experiment as E
import cgatcore.iotools as iotools

ORIGINAL_START = None

PARSER = None

EXPRESSIONS = (
    ('scripts', 'scripts/*.py'),)

EXCLUDE = ("__init__.py",
           "cgat.py",)


class DummyError(Exception):
    pass


def LocalStart(parser, *args, **kwargs):
    '''stub for E.start - set return_parser argument to true'''
    global PARSER
    PARSER = ORIGINAL_START(parser,
                            return_parser=True,
                            **kwargs
                            )
    raise DummyError()


def collectOptionsFromScript(script_name):
    '''collect options used in script *script_name*.'''

    prefix, suffix = os.path.splitext(script_name)

    dirname = os.path.dirname(script_name)
    basename = os.path.basename(script_name)[:-3]

    if os.path.exists(prefix + ".pyc"):
        os.remove(prefix + ".pyc")

    with iotools.open_file(script_name) as inf:
        if "getopt" in inf.read():
            E.warn("script %s uses getopt directly" % script_name)
            return []

    try:
        # Using importlib to load the module dynamically
        spec = importlib.util.spec_from_file_location(basename, script_name)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    except ImportError as msg:
        E.warn('could not import %s - skipped: %s' % (basename, msg))
        return []

    E.start = LocalStart

    try:
        module.main(argv=["--help"])
    except AttributeError:
        E.warn("no main method in %s" % script_name)
        return []
    except SystemExit:
        E.warn("script exits - possibly does not use E.start()")
        return []
    except DummyError:
        pass

    result = []
    for option in PARSER.option_list:
        if option.dest is None:
            continue

        optstring = option.get_opt_string()
        if optstring.startswith("--"):
            optstring = optstring[2:]
        result.append(optstring)

    return result


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--inplace", dest="inplace", action="store_true",
        help="update option list in place. New options will"
        "be added to the list given by --options-tsv-file. "
        "Options will only be added, not removed ")

    parser.add_argument(
        "--options-tsv-file", dest="tsv_file", type=str,
        help="existing table with options. Will be updated if "
        "--in-place is set [default]")

    parser.set_defaults(
        inplace=False,
        tsv_file=None)

    args = E.start(parser, argv=argv)

    old_options = None
    if args.tsv_file:
        if not os.path.exists(args.tsv_file):
            raise OSError(
                "filename %s not found, see --options-tsv-file" %
                args.tsv_file)
        old_options = pandas.read_csv(
            iotools.open_file(args.tsv_file),
            sep="\t",
            index_col=0,
        )
        old_options = old_options.fillna("")

    global ORIGINAL_START
    ORIGINAL_START = E.start

    all_options = collections.defaultdict(list)

    for label, expression in EXPRESSIONS:

        files = glob.glob(expression)
        files.sort()

        for f in files:

            E.debug("processing %s" % f)
            if os.path.isdir(f):
                continue
            if os.path.basename(f) in EXCLUDE:
                continue
            collected_options = collectOptionsFromScript(os.path.abspath(f))
            for o in collected_options:
                all_options[o].append(f)

    for x in old_options.index:
        if x not in all_options:
            all_options[x].append("--")

    if args.inplace:
        outfile = iotools.open_file(args.tsv_file, "w")
        E.info("updating file '%s'" % args.tsv_file)
    else:
        outfile = args.stdout

    outfile.write("option\taction\tcomment\talternative\tfiles\n")
    for o, v in sorted(all_options.items()):
        try:
            action, comment, alternative, ff = old_options.xs(o)
        except KeyError:
            action, comment, alternative, ff = "", "", "", ""

        if comment == "nan":
            comment = ""
        if alternative == "nan":
            alternative = ""

        outfile.write("\t".join((list(map(
            str, (o, action, comment, alternative, ",".join(v)))))) + "\n")

    if outfile != args.stdout:
        outfile.close()

    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
