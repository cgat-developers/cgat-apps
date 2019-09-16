'''
csvs2csv.py - join tables
=========================

:Tags: Python

Purpose
-------

This script reads several tab-separated tables and joins them.

.. note::
   working with multiple columns per table and sorting is
   not implemented correctly and likely to fail.

Usage
-----

Example::

   python combine_tables.py --help

Type::

   python combine_tables.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import os
import glob

import cgatcore.experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument(
        "-t", "--no-titles", dest="titles", action="store_false",
        help="no titles in input.")

    parser.add_argument(
        "-i", "--skip-titles", dest="skip_titles", action="store_true",
        help="skip output of titles.")

    parser.add_argument(
        "-m", "--missing-value", dest="missing_value", type=str,
        help="entry to use for missing values.")

    parser.add_argument("--header-names", dest="headers", type=str,
                        help="add headers for files.")

    parser.add_argument(
        "-c", "--columns", dest="columns", type=str,
        help="columns to use for joining. Multiple columns can be specified "
        "as a comma-separated list.")

    parser.add_argument(
        "-g", "--glob", dest="glob", type=str,
        help="wildcard expression for table names.")

    parser.add_argument(
        "-s", "--sort-order", dest="sort", type=str,
        help="sort by column titles alphabetical|numeric|list of columns.")

    parser.add_argument(
        "-e", "--merge-overlapping", dest="merge", action="store_true",
        help="simply merge tables without matching up rows. ")

    parser.add_argument(
        "--sort-keys", dest="sort_keys", type=str,
        choices=("numeric", "alphabetic"),
        help="sort key columns by value.")

    parser.add_argument(
        "--keep-empty", dest="ignore_empty", action="store_false",
        help="keep empty tables. The default is to ignore them.")

    parser.add_argument(
        "--add-file-prefix", dest="add_file_prefix", action="store_true",
        help="add file prefix to columns headers in multi-column tables ")

    parser.add_argument(
        "--regex-filename", dest="regex_filename", type=str,
        help="pattern to apply to filename to build prefix ")

    parser.set_defaults(
        titles=True,
        skip_titles=False,
        missing_value="na",
        headers=None,
        sort=None,
        glob=None,
        columns="1",
        sort_keys=False,
        merge=False,
        ignore_empty=True,
        add_file_prefix=False,
        regex_filename="(.*)"
    )

    (args, unknowns) = E.start(parser,
                               unknowns=True)

    if args.headers:
        if "," in args.headers:
            args.headers = args.headers.split(",")
        else:
            args.headers = re.split("\s+", args.headers.strip())

    if args.sort and args.sort not in ("numeric", "alphabetic"):
        if "," in args.sort:
            args.sort = args.sort.split(",")
        else:
            args.sort = re.split("\s+", args.sort)

    if args.merge:
        args.columns = []
    else:
        args.columns = [int(x) - 1 for x in args.columns.split(",")]

    args.filenames = []

    if args.glob:
        args.filenames += glob.glob(args.glob)

    args.filenames += unknown

    if len(args.filenames) < 1:
        print(USAGE, "no tables specified/found.")
        sys.exit(1)

    if args.loglevel >= 1:
        args.stdlog.write("# combining %i tables.\n" %
                          len(args.filenames))
        sys.stdout.flush()
        if len(args.filenames) == 1:
            for line in iotools.open_file(args.filenames[0]):
                args.stdout.write(line)
            E.stop()
            sys.exit(0)

    if args.headers and args.headers[0] != "auto" and \
       len(args.headers) != len(args.filenames):
        raise "number of provided headers (%i) is not equal to number filenames (%i)." %\
              (len(args.headers), len(args.filenames))

    tables = []
    keys = {}
    sorted_keys = []
    sizes = {}
    if args.merge:
        titles = ["count"]
    else:
        titles = []

    for filename in args.filenames:

        prefix = os.path.basename(filename)

        if os.path.exists(filename):
            file = iotools.open_file(filename, "r")
            lines = [x for x in file if x[0] != "#"]

        else:
            lines = []

        if len(lines) == 0 and args.ignore_empty:
            continue

        table = {}
        sizes = {}
        max_size = 0
        ncolumns = 0

        if args.titles:
            data = lines[0][:-1].split("\t")
            if not titles:
                key = "-".join([data[x] for x in args.columns])
                titles = [key]
            for x in range(len(data)):
                if x in args.columns:
                    continue
                ncolumns += 1
                if args.add_file_prefix:
                    p = re.search(args.regex_filename, prefix).groups()[0]
                    titles.append("%s_%s" % (p, data[x]))
                else:
                    titles.append(data[x])

            del lines[0]
        else:
            ncolumns = 1

        n = 0
        for line in lines:
            data = line[:-1].split("\t")
            row_keys = [data[x] for x in args.columns]
            if args.sort_keys:
                if args.sort_keys == "numeric":
                    row_keys.sort(lambda x, y: cmp(float(x), float(y)))
                else:
                    row_keys.sort()
            if args.merge:
                key = n
            else:
                key = "-".join(row_keys)

            if key not in keys:
                sorted_keys.append(key)
                keys[key] = 1
                sizes[key] = 0

            max_size = max(len(data) - len(args.columns), max_size)
            table[key] = [data[x]
                          for x in [x for x in range(0, len(data)) if x not in args.columns]]
            n += 1

        # enter columns of "na" for empty tables.
        if max_size == 0:
            max_size = ncolumns

        tables.append((max_size, table))

    if len(tables) == len(titles) - 1:

        if args.headers:
            headers = ["bin"]
            if args.headers[0] == 'auto':
                for t in range(len(tables)):
                    headers.append(os.path.basename(args.filenames[t]))
                    headers += [""] * (tables[t][0] - 1)

            else:
                for t in range(len(tables)):
                    headers.append(args.headers[t])
                    headers += [""] * (tables[t][0] - 1)

            # use headers as titles, if headers is given and skip-titles is
            # turned on
            if args.titles and args.skip_titles:
                titles = headers
            else:
                # otherwise: print the headers out right away
                sys.stdout.write("\t".join(headers) + "\n")

        order = list(range(0, len(tables) + 1))

        if args.titles:

            if args.sort:
                sort_order = []

                if args.sort == "numeric":
                    t = list(zip(list(map(int, titles[1:])), list(range(1, len(titles) + 1))))
                    t.sort()

                    for tt in t:
                        sort_order.append(titles[tt[1]])

                elif args.sort == "alphabetical":
                    t = list(zip(titles[1:], list(range(1, len(titles) + 1))))
                    t.sort()

                    for tt in t:
                        sort_order.append(titles[tt[1]])
                else:
                    sort_order = args.sort

                map_title2pos = {}
                for x in range(1, len(titles)):
                    map_title2pos[titles[x]] = x

                order = [0, ]
                for x in sort_order:
                    if x in map_title2pos:
                        order.append(map_title2pos[x])

            else:
                order = list(range(0, len(titles)))

            sys.stdout.write(
                "\t".join([titles[order[x]] for x in range(len(titles))]))
            sys.stdout.write("\n")

        if args.sort_keys:
            if args.sort_keys:
                if args.sort_keys == "numeric":
                    sorted_keys.sort(lambda x, y: cmp(float(x), float(y)))
                else:
                    sorted_keys.sort()

        for key in sorted_keys:

            sys.stdout.write("%s" % key)

            for x in order[1:]:
                max_size, table = tables[x - 1]
                c = 0
                if key in table:
                    sys.stdout.write("\t")
                    sys.stdout.write("\t".join(table[key]))
                    c = len(table[key])

                assert(max_size == 1)

                sys.stdout.write(
                    "\t%s" % args.missing_value * (max_size - c))

            sys.stdout.write("\n")

    else:

        # for multi-column table, just write
        if args.titles:
            sys.stdout.write(
                "\t".join([titles[x] for x in range(len(titles))]))
            sys.stdout.write("\n")

        for key in sorted_keys:

            sys.stdout.write("%s" % key)

            for x in range(len(tables)):

                max_size, table = tables[x]
                c = 0
                if key in table:
                    sys.stdout.write("\t")
                    sys.stdout.write("\t".join(table[key]))
                    c = len(table[key])

                sys.stdout.write(
                    "\t%s" % args.missing_value * (max_size - c))

            sys.stdout.write("\n")

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
