'''
csv2csv.py - operate on tables
==============================

:Tags: Python

Purpose
-------

operate on tables.

Usage
-----

Example::

   python csv2csv.py --help

Type::

   python csv2csv.py --help

for command line help.

Command line options
--------------------

'''
import sys
import csv
import cgatcore.experiment as E
import cgatcore.iotools as iotools


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument(
        "-s", "--method=sort --sort-order", dest="sort", type=str,
        help="fields to take (in sorted order).")

    (args) = E.start(parser, add_csv_options=True)

    reader = csv.DictReader(E.stdin, dialect=args.csv_dialect)

    if args.sort:
        fields = args.sort.split(",")
    else:
        fields = None

    writer = csv.DictWriter(E.stdout,
                            fields,
                            dialect=args.csv_dialect,
                            lineterminator=args.csv_lineterminator,
                            extrasaction='ignore')

    E.stdout.write("\t".join(fields) + "\n")

    for row in reader:
        row = iotools.convertDictionary(row)
        writer.writerow(row)

    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
