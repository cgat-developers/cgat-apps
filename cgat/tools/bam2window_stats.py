'''compute per-windows stats from a bam-file

Purpose
-------

This script takes a bam file as input and computes a few metrics by
iterating over the file. The metrics output are:

'''
import pysam

import cgatcore.experiment as E
from cgat.BamTools.bamtools import bam2stats_window_count


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "--region", dest="region", type="string",
        help="region to restrict analysis to [%default]")

    parser.add_option(
        "--window-size", dest="window_size", type="int",
        help="window size to use [%default]")

    parser.add_option(
        "--output-all-windows", dest="output_all_windows", action="store_true",
        help="output all windows. By default, windows without reads are skipped "
        "[%default]")

    parser.add_option(
        "--reference-fasta", "--input-filename-fasta",
        dest="input_filename_fasta", type="string",
        help="filename with reference sequence. If given, used to "
        "compute G+C content in windows [%default]")

    parser.set_defaults(
        force_output=False,
        region=None,
        output_all_windows=False,
        window_size=500,
        input_filename_fasta=None,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv, add_output_options=True)

    is_stdin = True
    if len(args) > 0:
        pysam_in = pysam.AlignmentFile(args[0], "rb")
        if args[0] != "-":
            is_stdin = False
    elif options.stdin == sys.stdin:
        pysam_in = pysam.AlignmentFile("-", "rb")
    else:
        pysam_in = pysam.AlignmentFile(options.stdin, "rb")
        if options.stdin != "-":
            is_stdin = False

    if options.input_filename_fasta:
        fasta = pysam.FastaFile(options.input_filename_fasta)
    else:
        fasta = None

    counts_df = bam2stats_window_count(
        pysam_in,
        region=options.region,
        window_size=options.window_size,
        fasta=fasta)

    if not options.output_all_windows:
        counts_df = counts_df[counts_df.alignments > 0]

    # add G+C content
    if fasta:
        counts_df["percent_gc"] = 100.0 * counts_df.bases_gc / (counts_df.bases_gc + counts_df.bases_at)
        counts_df.fillna(0, inplace=True)

    counts_df.to_csv(
        options.stdout,
        sep="\t")

    E.stop()
