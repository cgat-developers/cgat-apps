"""modify a fastq file.
"""

import collections
import sys
import pysam
import numpy
import CGATCore.Experiment as E
import CGATCore.IOTools as IOTools


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-i", "--input-fastq-file", dest="input_fastq_file", type="string",
        help="input fastq file. "
        "[%default]")

    parser.add_option(
        "--output-removed-tsv", dest="output_removed_tsv", type="string",
        help="if given, sequence identifiers of removed sequences will "
        "be stored in this file [%default]")

    parser.add_option(
        "--output-stats-tsv", dest="output_stats_tsv", type="string",
        help="if given, output statistics will be written to this file. "
        "[%default]")

    parser.add_option(
        "--output-removed-fastq", dest="output_removed_fastq", type="string",
        help="if given, removed fastq records will "
        "be stored in this file [%default]")

    parser.add_option(
        "-m", "--method", dest="methods", action="append", type="choice",
        choices=("filter-N", "filter-identifier", "filter-ONT",
                 "offset-quality"),
        help="methods to apply [%default]")

    parser.add_option(
        "--set-prefix", dest="set_prefix", type="string",
        help="set sequence prefix [%default]")

    parser.add_option(
        "--input-filter-tsv", dest="input_filter_tsv", type="string",
        help="list of sequence ides to filter [%default]")

    parser.add_option(
        "--min-average-quality", dest="min_average_quality", type="float",
        help="minimum average quality [%default]")

    parser.add_option(
        "--min-length", dest="min_length", type="int",
        help="minimum length [%default]")

    parser.add_option(
        "--quality-offset", dest="quality_offset", type="int",
        help="offset to modify quality values with [%default]")

    parser.set_defaults(
        methods=[],
        max_percent_N=10.0,
        input_fastq_file=None,
        set_prefix=None,
        output_removed_tsv=None,
        output_removed_fastq=None,
        output_stats_tsv=None,
        input_filter_tsv=None,
        min_average_quality=0,
        min_length=0,
        quality_offset=0,
    )

    (options, args) = E.start(parser, argv)

    if len(args) == 1:
        options.input_fastq_file = args[0]

    if options.input_fastq_file is None:
        raise ValueError("missing input fastq file")

    counter = E.Counter()

    filter_n = "filter-N" in options.methods

    filter_ont = "filter-ONT" in options.methods

    if "filter-identifier" in options.methods:
        if options.input_filter_tsv is None:
            raise ValueError("please set --input-filter-tsv for method filter-identifier")
        with IOTools.open_file(options.input_filter_tsv) as inf:
            filter_identifier = set([x.split()[0].strip() for x in inf.readlines()])
    else:
        filter_identifier = False

    if options.output_removed_tsv:
        outf_removed_tsv = IOTools.open_file(options.output_removed_tsv, "w")
    else:
        outf_removed_tsv = None

    if options.output_removed_fastq:
        outf_removed_fastq = IOTools.open_file(options.output_removed_fastq, "w")
    else:
        outf_removed_fastq = None

    if options.set_prefix:
        prefix = "{}".format(options.set_prefix)
    else:
        prefix = None

    quality_offset = options.quality_offset

    with pysam.FastxFile(options.input_fastq_file) as inf:
        for read in inf:
            counter.input += 1
            remove = False
            if filter_n:
                chars = collections.Counter(read.sequence)
                if "N" in chars and \
                   100.0 * chars["N"] / len(read.sequence) > options.max_percent_N:
                    remove = True
                    counter.filter_n += 1

            if filter_identifier:
                if read.name not in filter_identifier:
                    counter.filter_identifier += 1
                    remove = True

            if filter_ont:
                quals = read.get_quality_array()
                n = len(quals)
                if n < options.min_length or \
                        float(sum(quals)) / n < options.min_average_quality:
                    counter.remove_ont += 1
                    remove = True

            if remove:
                counter.removed += 1
                if outf_removed_tsv:
                    outf_removed_tsv.write(read.name + "\n")
                if outf_removed_fastq:
                    outf_removed_fastq.write(str(read) + "\n")
                continue

            if prefix:
                read.name = prefix + read.name[2:]

            if quality_offset:
                quals = numpy.array(read.get_quality_array())
                quals += quality_offset
                quals[quals < 0] = 0
                quals += 33
                # pysam fastq is read-only, so fudge it:
                # Note: not outputting description
                read = "@{}\n{}\n+\n{}".format(
                    read.name,
                    read.sequence,
                    "".join([chr(x) for x in quals]))

            counter.output += 1

            options.stdout.write(str(read) + "\n")

    if outf_removed_tsv:
        outf_removed_tsv.close()

    if outf_removed_fastq:
        outf_removed_fastq.close()

    if options.output_stats_tsv:
        with IOTools.open_file(options.output_stats_tsv, "w") as outf:
            outf.write(counter.asTable(as_rows=False) + "\n")

    E.info(counter)
    E.stop()
