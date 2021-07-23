'''
data2histogram.py - histogram data in a table
=============================================

:Tags: Python

Purpose
-------

This script computes histograms over one or more
columns of a table.

Usage
-----

Example::

   python data2histogram.py --help

Type::

   python data2histogram.py --help

for command line help.

Command line options
--------------------

'''
import sys
import cgatcore.experiment as E
import cgat.Histogram as Histogram
import numpy


def main(argv=None):

    if not argv:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("-r", "--range", dest="range", type=str,
                        help="range to calculate histogram for.")
    parser.add_argument("-b", "--bin-size", dest="bin_size", type=str,
                        help="bin size.")
    parser.add_argument("-i", "--titles", dest="titles", action="store_true",
                        help="use supplied column titles.")
    parser.add_argument("--no-null", dest="nonull", action="store_true",
                        help="do not output null values")
    parser.add_argument("--no-titles", dest="titles", action="store_false",
                        help="no column titles given.")
    parser.add_argument("-c", "--columns", dest="columns", type=str,
                        help="columns to take for calculating histograms.")
    parser.add_argument("--min-data", dest="min_data", type=int,
                        help="minimum amount of data required, if less data, then the histogram will be empty [default=%default].")
    parser.add_argument("--min-value", dest="min_value", type=float,
                        help="minimum value for histogram.")
    parser.add_argument("--max-value", dest="max_value", type=float,
                        help="maximum value for histogram.")
    parser.add_argument("--no-empty-bins", dest="no_empty_bins", action="store_true",
                        help="do not display empty bins.")
    parser.add_argument("--with-empty-bins", dest="no_empty_bins", action="store_false",
                        help="display empty bins.")
    parser.add_argument("--normalize", dest="normalize", action="store_true",
                        help="normalize histogram.")
    parser.add_argument("--cumulative", dest="cumulative", action="store_true",
                        help="calculate cumulative histogram.")
    parser.add_argument("--reverse-cumulative", dest="reverse_cumulative", action="store_true",
                        help="calculate reverse cumulative histogram.")
    parser.add_argument("--header-names", dest="headers", type=str,
                        help="use the following headers.")
    parser.add_argument("--ignore-out-of-range", dest="ignore_out_of_range", action="store_true",
                        help="ignore values that are out of range (as opposed to truncating them to range border.")
    parser.add_argument("--missing-value", dest="missing_value", type=str,
                        help="entry for missing values .")
    parser.add_argument("--use-dynamic-bins", dest="dynamic_bins", action="store_true",
                        help="each value constitutes its own bin.")
    parser.add_argument("--on-the-fly", dest="on_the_fly", action="store_true",
                        help="on the fly computation of histograms. Requires setting of min-value, max-value and bin_size.")

    parser.set_defaults(
        bin_size=None,
        range=None,
        titles=True,
        columns="all",
        append=(),
        no_empty_bins=True,
        min_value=None,
        max_value=None,
        normalize=False,
        cumulative=False,
        reverse_cumulative=False,
        nonull=None,
        ignore_out_of_range=False,
        min_data=1,
        headers=None,
        missing_value="na",
        dynamic_bins=False,
        on_the_fly=False,
        bin_format="%.2f",
        value_format="%6.4f",
    )

    (args) = E.start(parser)

    if args.columns != "all":
        args.columns = [int(x) - 1 for x in args.columns.split(",")]

    if args.range:
        args.min_value, args.max_value = list(map(
            float, args.range.split(",")))

    if args.headers:
        args.headers = args.headers.split(",")

    if args.on_the_fly:
        if args.min_value is None or args.max_value is None or \
           args.bin_size is None:
            raise ValueError("please supply columns, min-value, max-value and "
                             "bin-size for on-the-fly computation.")

        # try to glean titles from table:
        if args.titles:
            while 1:
                line = sys.stdin.readline()
                if not line:
                    break
                if line[0] == "#":
                    continue
                data = line[:-1].split("\t")
                break

            if args.columns == "all":
                args.titles = data
                args.columns = list(range(len(data)))
            else:
                args.titles = [data[x] for x in args.columns]

        bins = numpy.arange(
            args.min_value, args.max_value, float(args.bin_size))
        hh = Histogram.fillHistograms(
            sys.stdin, args.columns, [bins for x in range(len(args.columns))])
        n = len(hh)

        titles = ['bin']

        if args.headers:
            titles.append(args.headers[x])
        elif args.titles:
            titles.append(args.titles[x])
        else:
            for x in args.columns:
                titles.append("col%i" % (x + 1))

        if len(titles) > 1:
            args.stdout.write("\t".join(titles) + "\n")

        for x in range(len(bins)):
            v = []
            v.append(args.bin_format % bins[x])
            for c in range(n):
                v.append(args.value_format % hh[c][x])

            args.stdout.write("\t".join(v) + "\n")

    else:
        # in-situ computation of histograms
        # retrieve data
        first = True
        vals = []

        # parse data, convert to floats
        for l in args.stdin:

            if l[0] == "#":
                continue

            data = l[:-1].split("\t")

            if first:
                first = False
                ncols = len(data)
                if args.columns == "all":
                    args.columns = list(range(ncols))

                vals = [[] for x in args.columns]

                if args.titles:
                    try:
                        args.titles = [data[x] for x in args.columns]
                    except IndexError:
                        raise IndexError("not all columns %s found in data %s" % (
                            str(args.columns), str(data)))
                    continue

            for x in range(len(args.columns)):

                try:
                    v = float(data[args.columns[x]])
                except IndexError:
                    print("# IndexError in line:", l[:-1])
                    continue
                except ValueError:
                    continue

                vals[x].append(v)

        lines = None

        hists = []
        titles = []

        if not vals:
            if args.loglevel >= 1:
                args.stdlog.write("# no data\n")
            E.stop()
            sys.exit(0)

        for x in range(len(args.columns)):

            if args.loglevel >= 1:
                args.stdlog.write(
                    "# column=%i, num_values=%i\n" % (args.columns[x], len(vals[x])))

            if len(vals[x]) < args.min_data:
                continue

            h = Histogram.Calculate(vals[x],
                                    no_empty_bins=args.no_empty_bins,
                                    increment=args.bin_size,
                                    min_value=args.min_value,
                                    max_value=args.max_value,
                                    dynamic_bins=args.dynamic_bins,
                                    ignore_out_of_range=args.ignore_out_of_range)

            if args.normalize:
                h = Histogram.Normalize(h)
            if args.cumulative:
                h = Histogram.Cumulate(h)
            if args.reverse_cumulative:
                h = Histogram.Cumulate(h, direction=0)

            hists.append(h)

            for m in args.append:
                if m == "normalize":
                    hists.append(Histogram.Normalize(h))

            if args.headers:
                titles.append(args.headers[x])
            elif args.titles:
                titles.append(args.titles[x])
            else:
                titles.append("col%i" % args.columns[x])

        if titles:
            args.stdout.write("bin\t" + "\t".join(titles) + "\n")

        if len(hists) == 1:
            Histogram.Print(hists[0], nonull=args.nonull,
                            format_bin=args.bin_format)
        else:
            combined_histogram = Histogram.Combine(
                hists, missing_value=args.missing_value)
            Histogram.Print(combined_histogram,
                            nonull=args.nonull,
                            format_bin=args.bin_format)

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
