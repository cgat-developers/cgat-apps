"""compute statistics on a table
================================

This tool outputs summary metrics for a tab-separated table.

The table is read into memory, so there is a limit to the size of
table that can be analyzed.

The tool outputs four columns:

metric
    the name of the metric
count
    number of entities
percent
    percent value of metric
info
    additional information, typically the denominator that the
    percent value is computed from

"""

import sys
import re
import pandas
import cgatcore.experiment as E
import cgatcore.iotools as iotools


def compute_table_summary(table):

    nrows = len(table)
    ncolumns = len(table.columns)
    ncells = nrows * ncolumns
    yield "rows", nrows, nrows, "rows"
    yield "columns", ncolumns, ncolumns, "columns"
    ncells = nrows * ncolumns
    yield "cells", ncells, ncells, "cells"
    yield "rows_with_na", nrows - len(table.dropna(axis=0)), nrows, "rows"
    yield "columns_with_na", ncolumns - len(table.dropna(axis=1).columns), ncolumns, "columns"
    yield "cells_with_na", sum(table.isnull().values.ravel()), ncells, "cells"

    for column in table.columns:
        nna = table[column].dropna()
        yield "column_na", nrows - len(nna), nrows, column
        yield "column_unique", len(nna.unique()), None, column


def main(argv=None):

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument(
        "-d", "--delimiter", dest="delimiter", type=str,
        help="delimiter to separate columns ")

    parser.add_argument(
        "-m", "--method", dest="methods", type=str,
        action="append",
        choices=["row-describe", "column-describe"],
        help="additional methods to apply ")

    parser.set_defaults(
        delimiter="\t",
        methods=[],
    )

    (args) = E.start(parser,
                     argv=argv,
                     add_output_options=True)

    if not args.methods:
        args.methods = ["summary"]

    table = pandas.read_csv(args.stdin, args.delimiter)

    args.stdout.write("metric\tcount\tpercent\tinfo\n")

    for method in args.methods:
        label = re.sub("-", "_", method)
        if method == "summary":
            for category, count, denominator, info in compute_table_summary(table):
                args.stdout.write("\t".join(map(str, (
                    category,
                    count,
                    iotools.pretty_percent(count, denominator, na=""),
                    info))) + "\n")
        elif method == "column-describe":
            df = table.describe().T.stack()
            with E.open_output_file(label) as outf:
                outf.write("label\tcategory\tvalue\n")
                df.to_csv(outf, sep="\t")
        elif method == "row-describe":
            df = table.T.describe().stack()
            with E.open_output_file(label) as outf:
                outf.write("label\tcategory\tvalue\n")
                df.to_csv(outf, sep="\t")

    E.stop()


if __name__ == "__main__":
    sys.exit(main())
