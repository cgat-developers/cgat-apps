"""
bed2graph.py - compute the overlap graph between two bed files
==============================================================

:Tags: Python

Purpose
-------

This script ouputs a list of the names of all overlapping intervals
between two bed files.

Usage
-----

Type::

   python bed2graph.py A.bed.gz B.bed.gz > graph.out

for command line help.

Command line options
--------------------

"""

import sys
import cgatcore.experiment as E
import cgatcore.iotools as iotools
import cgat.Bed as Bed


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("-o", "--output-section", dest="output", type=str,
                        choices=("full", "name"),
                        help="output either ``full`` overlapping entries, only the ``name``s.")

    parser.set_defaults(
        output="full",
    )

    # add common options (-h/--help, ...) and parse command line
    (args, unknown) = E.start(parser,
                              argv=argv,
                              unknowns=True)

    if len(unknown) != 2:
        raise ValueError("two arguments required")

    if unknown[0] == "-":
        infile1 = args.stdin
    else:
        infile1 = iotools.open_file(unknown[0], "r")

    infile2 = iotools.open_file(unknown[1], "r")

    idx = Bed.readAndIndex(infile2, with_values=True)

    output = args.output
    outfile = args.stdout

    if output == "name":
        outfile.write("name1\tname2\n")
        outf = lambda x: x.fields[0]
    else:
        outf = str

    for bed in Bed.iterator(infile1):
        try:
            overlaps = idx[bed.contig].find(bed.start, bed.end)
        except (KeyError, IndexError):
            # ignore missing contig and zero length intervals
            continue

        for o in overlaps:
            outfile.write("\t".join((outf(bed), outf(o[2]))) + "\n")

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
