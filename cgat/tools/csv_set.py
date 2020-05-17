'''
csv_set.py - set operations on a table
======================================

:Tags: Python

Purpose
-------

.. todo::

   describe purpose of the script.

Usage
-----

Example::

   python csv_set.py --help

Type::

   python csv_set.py --help

for command line help.

Command line options
--------------------

'''
import sys

import cgatcore.experiment as E
from cgatcore.csvutils import readTable
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

    parser.add_argument("-1", "--join-fields1", dest="join_fields1", type=str,
                        help="join fields in first table.")
    parser.add_argument("-2", "--join-fields2", dest="join_fields2", type=str,
                        help="join fields in second table.")
    parser.add_argument("-m", "--method", dest="method", type=str,
                        help="set operation to perform.", choices=("intersection", "rest", "union"))

    parser.set_defaults(
        remove=False,
        unique=False,
        join_fields1=None,
        join_fields2=None,
        method="intersection",
    )

    (args, unknown) = E.start(parser,
                              add_csv_options=True,
                              unknowns=True)

    if len(unknown) != 2:
        raise ValueError("please specify two files to join")

    if not args.join_fields1 or not args.join_fields2:
        raise ValueError("please specify at least one join field per table")

    args.join_fields1 = args.join_fields1.split(",")
    args.join_fields2 = args.join_fields2.split(",")

    args.filename1, args.filename2 = unknown

    fields1, table1 = readTable(open(args.filename1, "r"))
    fields2, table2 = readTable(open(args.filename2, "r"))

    if args.unique:
        outfile = UniqueBuffer(sys.stdout)
    else:
        outfile = args.stdout

    nfields1 = []
    for x in range(len(fields1)):
        if fields1[x] in args.join_fields1:
            nfields1.append(x)
    nfields2 = []
    for x in range(len(fields2)):
        if fields2[x] in args.join_fields2:
            nfields2.append(x)

    # calculate row indices: double keys are not taken care of here
    keys = {}
    for row1 in table1:
        v = [row1[x] for x in nfields1]
        key = hashlib.md5("".join(v)).digest()
        keys[key] = row1

    if args.method == "intersection":
        # build new field list
        take = list(range(len(fields1)))
        c = len(take)
        for x in fields2:
            if x not in args.join_fields2:
                take.append(c)
            c += 1

        t = fields1 + fields2

        new_fields = [t[x] for x in take]

        print("\t".join(new_fields))

        for row2 in table2:
            v = [row2[x] for x in nfields2]
            key = hashlib.md5("".join(v)).digest()
            if key in keys:
                new_row = keys[key] + row2
                outfile.write(
                    "\t".join([new_row[x] for x in take]) + "\n")

    elif args.method == "rest":

        new_fields = fields2
        print("\t".join(new_fields))

        for row2 in table2:
            v = [row2[x] for x in nfields2]
            key = hashlib.md5("".join(v)).digest()
            if key not in keys:
                outfile.write("\t".join(row2) + "\n")

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
