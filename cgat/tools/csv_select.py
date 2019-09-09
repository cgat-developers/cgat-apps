'''
csv_select.py - select rows from a table
========================================

:Tags: Python

Purpose
-------

extract rows from a csv-formatted table.

The select statement is a one-line, for example::

   csv_select.py "int(r['mC-foetal-sal-R4']) > 0" < in > out

Note the required variable name r for denoting field names. Please
also be aware than numeric values need to be converted first before
testing.

Usage
-----

Type::

   python csv_select.py --help

for command line help.

Command line options
--------------------

'''
import sys
import csv
import _csv
import cgatcore.experiment as E
from cgatcore.csvutils import CommentStripper, DictReaderLarge


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("-r", "--remove", dest="remove", action="store_true",
                        help="remove specified columns, keep all others.")

    parser.add_argument("-u", "--unique", dest="unique", action="store_true",
                        help="output rows are uniq.")

    parser.add_argument("-l", "--large", dest="large", action="store_true",
                        help="large columns. Do not use native python csv module [default=%default].")

    parser.add_argument("-f", "--filename-fields", dest="filename_fields", type=str,
                        help="filename with field information.")

    parser.set_defaults(
        remove=False,
        unique=False,
        filename_fields=None,
    )

    (args, unknown) = E.start(parser,
                              add_csv_options=True,
                              quiet=True,
                              unknowns=True)

    statement = " ".join(unknown)

    if args.large:
        reader = DictReaderLarge(CommentStripper(sys.stdin),
                                 dialect=args.csv_dialect)
    else:
        reader = csv.DictReader(CommentStripper(sys.stdin),
                                dialect=args.csv_dialect)

    exec("f = lambda r: %s" % statement, globals())
    counter = E.Counter()
    writer = csv.DictWriter(args.stdout,
                            reader.fieldnames,
                            dialect=args.csv_dialect,
                            lineterminator=args.csv_lineterminator)

    writer.writerow(dict((fn, fn) for fn in reader.fieldnames))
    while 1:
        counter.input += 1
        try:
            row = next(reader)
        except _csv.Error as msg:
            args.stderr.write("# error while parsing: %s\n" % (msg))
            counter.errors += 1
            continue
        except StopIteration:
            break

        if not row:
            break

        if f(row):
            writer.writerow(row)
            counter.output += 1
        else:
            counter.filtered += 1

    E.info("%s" % counter)

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
