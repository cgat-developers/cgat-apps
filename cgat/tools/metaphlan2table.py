'''
script_template.py
=============================================

:Tags: Python

Purpose
-------

Convert the output of a metaphlan analysis to a preferred table format


Usage
-----

Example::

   python metaphlan2table.py --help

Type::

   python metaphlan2table.py --help

for command line help.

Documentation
-------------

Code
----

'''

import sys
import argparse
import cgat.Metaphlan as Metaphlan
import cgatcore.experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("-t", "--sequence-type", dest="type", type=str,
                        choices=("read_map", "rel_ab"), help="type of file to be parsed to a table")

    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    assert args.type, "must specify infile type"
    if args.type == "read_map":
        args.stdout.write(
            "seq_id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
        for entry in Metaphlan.read_map_iterator(sys.stdin):
            args.stdout.write("\t".join(
                [entry.seq_id, entry.kingdom, entry.phylum, entry.c_lass, entry.order, entry.family, entry.genus, entry.species]) + "\n")

    elif args.type == "rel_ab":
        args.stdout.write("taxon_level\ttaxon\trel_abundance\n")
        for entry in Metaphlan.relative_abundance_iterator(sys.stdin):
            args.stdout.write(
                "\t".join([entry.taxon_level, entry.taxon, entry.abundance]) + "\n")

    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
