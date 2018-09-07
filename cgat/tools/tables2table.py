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

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--no-titles",
                      dest="input_has_titles",
                      action="store_false",
                      help="no titles in input [%default].")

    parser.add_option("--ignore-titles",
                      dest="ignore_titles",
                      action="store_true",
                      help="ignore titles in input [%default]")

    parser.add_option("-i", "--skip-titles",
                      dest="skip_titles",
                      action="store_true",
                      help="skip output of titles.")

    parser.add_option("-m", "--missing-value",
                      dest="missing_value",
                      type="string",
                      help="entry to use for missing values.")

    parser.add_option("--header-names", dest="headers", type="string",
                      help="add headers for files as a ,-separated "
                      "list [%default].")

    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to use for joining. Multiple columns "
                      "can be specified as a comma-separated list "
                      "[default=%default].")

    parser.add_option("-k", "--take",
                      dest="take",
                      type="string",
                      action="append",
                      help="columns to take. If not set, all columns "
                      "except for "
                      "the join columns are taken [%default]")

    parser.add_option("-g", "--glob", dest="glob", type="string",
                      help="wildcard expression for table names.")

    parser.add_option(
        "-s", "--sort-order", dest="sort", type="string",
        help="sort by column titles in particular given order: "
        "alphabetical|numeric|list of columns.")

    parser.add_option(
        "-e", "--merge-overlapping", dest="merge", action="store_true",
        help="simply merge tables without matching up "
        "rows. [default=%default].")

    parser.add_option("-a", "--cat", dest="cat", type="string",
                      help="simply concatenate tables. Adds an "
                      "additional column called X with the filename "
                      " [default=%default].")

    parser.add_option("--sort-keys", dest="sort_keys", type="choice",
                      choices=("numeric", "alphabetic"),
                      help="sort key columns by value.")

    parser.add_option("--keep-empty", dest="ignore_empty",
                      action="store_false",
                      help="keep empty tables. The default is "
                      "to ignore them.")

    parser.add_option("--ignore-empty",
                      dest="ignore_empty",
                      action="store_true",
                      help="ignore empty tables - this is "
                      "the default [%default].")

    parser.add_option("--add-file-prefix",
                      dest="add_file_prefix",
                      action="store_true",
                      help="add file prefix to "
                      "columns headers. Suitable for multi-column"
                      "tables [default=%default]")

    parser.add_option("--use-file-prefix",
                      dest="use_file_prefix",
                      action="store_true",
                      help="use file prefix as column headers. "
                      "Suitable for two-column tables "
                      "[default=%default]")

    parser.add_option("--prefixes", dest="prefixes", type="string",
                      help="list of prefixes to use. "
                      ", separated list of prefixes. "
                      "The number of prefixes need to correspond to the "
                      "number of input files [default=%default]")

    parser.add_option("--regex-filename", dest="regex_filename",
                      type="string",
                      help="pattern to apply to filename to "
                      "build prefix [default=%default]")

    parser.add_option("--regex-start",
                      dest="regex_start",
                      type="string",
                      help="regular expression to start "
                      "collecting table in a file [default=%default]")

    parser.add_option("--regex-end",
                      dest="regex_end",
                      type="string",
                      help="regular expression to end collecting "
                      "table in a file [default=%default]")

    parser.add_option("--sep",
                      dest="separator",
                      type="string",
                      help="table separator to use. The default is to use tabs. "
                      "[default=%default]")

    parser.add_option("--test", dest="test",
                      type="int",
                      help="test combining tables with "
                      "first X rows [default=%default]")

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

    (options, args) = E.start(parser, argv=argv)

    if options.headers:
        if "," in options.headers:
            options.headers = options.headers.split(",")
        else:
            options.headers = re.split("\s+", options.headers.strip())

    if options.sort and options.sort not in ("numeric", "alphabetic"):
        if "," in options.sort:
            options.sort = options.sort.split(",")
        else:
            options.sort = re.split("\s+", options.sort)

    if options.merge:
        options.columns = []
    else:
        options.columns = [int(x) - 1 for x in options.columns.split(",")]

    options.filenames = []

    if options.glob:
        options.filenames += glob.glob(options.glob)

    options.filenames += args

    if len(options.filenames) < 1:
        raise ValueError("no tables found.")

    E.info("combining %i tables" % len(options.filenames))

    if options.cat:
        table = concatenate_tables(options.filenames,
                                   regex_filename=options.regex_filename,
                                   separator=options.separator,
                                   headers=options.headers,
                                   missing_value=options.missing_value,
                                   cat=options.cat)

    table.to_csv(options.stdout, sep=options.separator, index=False)
    E.stop()


if __name__ == '__main__':
    sys.exit(main(sys.argv))
