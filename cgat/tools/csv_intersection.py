'''
csv_intersection.py - intersect two tables
======================================================

:Tags: Python

Purpose
-------

.. todo::

   describe purpose of the script.

Usage
-----

Example::

   python csv_intersection.py --help

Type::

   python csv_intersection.py --help

for command line help.

Command line options
--------------------

'''
import sys
import cgatcore.experiment as E
import cgatcore.iotools as iotools
from cgatcore.csvutils import readTable
import csv
import hashlib


class UniqueBuffer:
    mKeys = {}

    def __init__(self, outfile):
        self.mOutfile = outfile

    def write(self, out):
        key = hashlib.md5(out).digest()
        if key not in self.mKeys:
            self.mKeys[key] = True
            self.mOutfile.write(out)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("-u", "--unique", dest="unique", action="store_true",
                        help="output rows are uniq.")

    parser.set_defaults(
        remove=False,
        unique=False,
    )

    (args, unknown) = E.start(parser, add_csv_options=True, unknowns=True)

    if len(args) != 2:
        raise ValueError("please specify two files to join")

    args.filename1, args.filename2 = unknown

    table1 = readTable(iotools.open_file(args.filename1, "r"))
    table2 = readTable(iotools.open_file(args.filename2, "r"))

    if args.unique:
        outfile = UniqueBuffer(sys.stdout)
    else:
        outfile = args.stdout

    # build new field list
    new_fields = []

    for x in args.join_fields1:
        new_fields.append(x)

    for x in fields1:
        if x not in args.join_fields1:
            new_fields.append(x)
        if x not in args.join_fields2:
            new_fields.append(x)

        writer = csv.DictWriter(outfile,
                                fields,
                                dialect=args.csv_dialect,
                                lineterminator=args.csv_lineterminator,
                                extrasaction='ignore')

    if len(lines) > 0:

        old_fields = lines[0][:-1].split("\t")

        if args.remove:
            fields = []
            for x in old_fields:
                if x not in input_fields:
                    fields.append(x)
        else:
            fields = input_fields

        reader = csv.DictReader(lines,
                                dialect=args.csv_dialect)

        print("\t".join(fields))

        first_row = True
        for row in reader:
            row = iotools.convertDictionary(row)
            writer.writerow(row)

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
