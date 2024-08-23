'''combine_tables.py - join tables
==================================

:Tags: Python

Purpose
-------

This script reads several tab-separated tables and joins them into a
single one.

Usage
-----

The option ``--header-names`` sets the column titles explicitely. Add
``--skip-titles`` if you want to avoid echoing the original title in
the input files.


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
import collections
import pandas

import cgatcore.iotools as iotools
import cgatcore.experiment as E


def read_tables(filenames, *args, **kwargs):

    tables = []
    for filename in filenames:
        try:
            table = pandas.read_csv(filename, *args, **kwargs)
        except pandas.errors.EmptyDataError:
            E.warn("file '{}' is empty".format(filename))
            continue
        except pandas.errors.ParserError as ex:
            E.warn("file '{}' has parsing error: {}".format(filename, ex))
            continue
        if len(table) == 0:
            E.warn("table '{}' is empty".format(filename))
            continue
        tables.append((table, filename))
    return zip(*tables)


def concatenate_tables(filenames,
                       regex_filename="(\S+)",
                       separator="\t",
                       headers=None,
                       missing_value=None,
                       cat=None):

    '''concatenate tables.'''

    rx = re.compile(regex_filename)

    tables, filenames = read_tables(filenames, sep=separator)

    if headers is None or headers == "auto":
        row_headers = [
            [y for y in rx.search(x).groups()] for x in filenames]
    else:
        row_headers = [headers]

    if cat is None:
        if len(row_headers) == 1:
            row_head_titles = ["filename"]
        else:
            row_head_titles = ["pattern" + str(x) for x in range(len(row_headers))]
    else:
        row_head_titles = [x.strip() for x in cat.split(",")]
        if len(row_headers[0]) != len(row_head_titles):
            raise ValueError(
                "row header (%i) has different number of fields in "
                "regular expression than supplied by the --cat option (%i)" %
                (len(row_headers[0]), len(row_head_titles)))

    # avoid MultiIndex if only single level
    if len(row_head_titles) == 1:
        row_head_titles = row_head_titles
        row_headers = [x[0] for x in row_headers]

    df = pandas.concat(tables, axis=0, keys=row_headers, names=row_head_titles) \
               .reset_index() \
               .drop(["level_1"], axis=1)
    return df


def main(argv=sys.argv):

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("-t", "--no-titles",
                        dest="input_has_titles",
                        action="store_false",
                        help="no titles in input .")

    parser.add_argument("--ignore-titles",
                        dest="ignore_titles",
                        action="store_true",
                        help="ignore titles in input ")

    parser.add_argument("-i", "--skip-titles",
                        dest="skip_titles",
                        action="store_true",
                        help="skip output of titles.")

    parser.add_argument("-m", "--missing-value",
                        dest="missing_value",
                        type=str,
                        help="entry to use for missing values.")

    parser.add_argument("--header-names", dest="headers", type=str,
                        help="add headers for files as a ,-separated "
                        "list .")

    parser.add_argument("-c", "--columns", dest="columns", type=str,
                        help="columns to use for joining. Multiple columns "
                        "can be specified as a comma-separated list ")

    parser.add_argument("-k", "--take",
                        dest="take",
                        type=str,
                        action="append",
                        help="columns to take. If not set, all columns "
                        "except for "
                        "the join columns are taken ")

    parser.add_argument("-g", "--glob", dest="glob", type=str,
                        help="wildcard expression for table names.")

    parser.add_argument(
        "-s", "--sort-order", dest="sort", type=str,
        help="sort by column titles in particular given order: "
        "alphabetical|numeric|list of columns.")

    parser.add_argument(
        "-e", "--merge-overlapping", dest="merge", action="store_true",
        help="simply merge tables without matching up "
        "rows.")

    parser.add_argument("-a", "--cat", dest="cat", type=str,
                        help="simply concatenate tables. Adds an "
                        "additional column called X with the filename ")

    parser.add_argument("--sort-keys", dest="sort_keys", type=str,
                        choices=("numeric", "alphabetic"),
                        help="sort key columns by value.")

    parser.add_argument("--keep-empty", dest="ignore_empty",
                        action="store_false",
                        help="keep empty tables. The default is "
                        "to ignore them.")

    parser.add_argument("--ignore-empty",
                        dest="ignore_empty",
                        action="store_true",
                        help="ignore empty tables - this is "
                        "the default .")

    parser.add_argument("--add-file-prefix",
                        dest="add_file_prefix",
                        action="store_true",
                        help="add file prefix to "
                        "columns headers. Suitable for multi-column"
                        "tables")

    parser.add_argument("--use-file-prefix",
                        dest="use_file_prefix",
                        action="store_true",
                        help="use file prefix as column headers. "
                        "Suitable for two-column tables ")

    parser.add_argument("--prefixes", dest="prefixes", type=str,
                        help="list of prefixes to use. "
                        ", separated list of prefixes. "
                        "The number of prefixes need to correspond to the "
                        "number of input files")

    parser.add_argument("--regex-filename", dest="regex_filename",
                        type=str,
                        help="pattern to apply to filename to "
                        "build prefix")

    parser.add_argument("--regex-start",
                        dest="regex_start",
                        type=str,
                        help="regular expression to start "
                        "collecting table in a file")

    parser.add_argument("--regex-end",
                        dest="regex_end",
                        type=str,
                        help="regular expression to end collecting "
                        "table in a file")

    parser.add_argument("--sep",
                        dest="separator",
                        type=str,
                        help="table separator to use. The default is to use tabs. ")

    parser.add_argument("--test", dest="test",
                        type=int,
                        help="test combining tables with "
                        "first X rows")

    parser.set_defaults(
        input_has_titles=True,
        skip_titles=False,
        missing_value=None,
        headers=None,
        sort=None,
        glob=None,
        columns="1",
        sort_keys=False,
        merge=False,
        ignore_empty=True,
        regex_start=None,
        regex_end=None,
        add_file_prefix=False,
        use_file_prefix=False,
        cat=None,
        take=[],
        regex_filename="(.*)",
        prefixes=None,
        test=0,
        separator="\t",
    )

    (args, unknown) = E.start(parser,
                              argv=argv,
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
        raise ValueError("no tables found.")

    E.info("combining %i tables" % len(args.filenames))

    # Remove the if statement and call concatenate_tables directly
    table = concatenate_tables(args.filenames,
                               regex_filename=args.regex_filename,
                               separator=args.separator,
                               headers=args.headers,
                               missing_value=args.missing_value,
                               cat=args.cat)

    # Ensure the table object is not None before attempting to write it
    if table is not None:
        table.to_csv(args.stdout, sep=args.separator, index=False)
    else:
        E.warn("No tables were concatenated.")

    E.stop()


if __name__ == '__main__':
    sys.exit(main(sys.argv))
