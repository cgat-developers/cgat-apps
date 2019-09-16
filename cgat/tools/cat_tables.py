"""
cat_tables.py - concatenate tables
==================================

:Tags: Python

Purpose
-------

concatenate tables. Headers of subsequent files are ignored.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Command line options
--------------------

"""

import sys
import fileinput

import cgatcore.experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.set_defaults(
    )

    # add common options (-h/--help, ...) and parse command line
    (args, unknown) = E.start(parser,
                              argv=argv,
                              unknowns=True)

    if len(unknown) == 0 or (len(unknown) == 1 and unknown[0] == "-"):
        infile = args.stdin
    else:
        infile = fileinput.FileInput(args)

    # do sth
    ninput, nskipped, noutput = 0, 0, 0

    header = False

    for line in infile:
        ninput += 1
        if line.startswith("#"):
            pass
        elif not header:
            header = line
        elif line == header:
            nskipped += 1
            continue

        args.stdout.write(line)
        noutput += 1

    E.info("ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped))

    # write footer and output benchmark information.
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
