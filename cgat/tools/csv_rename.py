'''
csv_rename.py - rename columns in a table
=========================================

:Tags: Python

Purpose
-------

rename columns in a csv file

Usage
-----

Example::

   csv_rename.py gene=id < stdin

Type::

   python csv_rename.py --help

for command line help.

Command line options
--------------------

'''
import sys
import cgatcore.experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("-r", "--remove", dest="remove", action="store_true",
                        help="remove specified columns, keep all others.")

    parser.add_argument("-u", "--unique", dest="unique", action="store_true",
                        help="output rows are uniq.")

    parser.add_argument("-f", "--filename-fields", dest="filename_fields", type=str,
                        help="filename with field information.")

    parser.set_defaults(
        filename_fields=None,
    )

    (args) = E.start(parser,
                     add_csv_options=True)
    mapper = {}
    for x in args:
        a, b = x.split("=")
        mapper[a.strip()] = b.strip()

    while 1:
        line = args.stdin.readline()

        if not line:
            E.stop()
            sys.exit(0)

        if line[0] == "#":
            args.stdout.write(line)
            continue

        break

    header = []
    nreplaced = 0
    for x in line[:-1].split():
        if x in mapper:
            nreplaced += 1
            header.append(mapper[x])
        else:
            header.append(x)

    args.stdout.write("\t".join(header) + "\n")
    nlines = 0
    for line in args.stdin:
        nlines += 1
        args.stdout.write(line)

    if args.loglevel >= 1:
        ninput = len(header)
        noutput = ninput
        args.stdout.write("# ninput=%i, noutput=%i, nreplaced=%i, nlines=%i\n" % (
            ninput, noutput, nreplaced, nlines))

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
