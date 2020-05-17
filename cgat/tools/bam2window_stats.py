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
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--region", dest="region", type=str,
        help="region to restrict analysis to ")

    parser.add_argument(
        "--window-size", dest="window_size", type=int,
        help="window size to use ")

    parser.add_argument(
        "--output-all-windows", dest="output_all_windows", action="store_true",
        help="output all windows. By default, windows without reads are skipped "
        )

    parser.add_argument(
        "--reference-fasta", "--input-filename-fasta",
        dest="input_filename_fasta", type=str,
        help="filename with reference sequence. If given, used to "
        "compute G+C content in windows ")

    parser.set_defaults(
        force_output=False,
        region=None,
        output_all_windows=False,
        window_size=500,
        input_filename_fasta=None,
    )

    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv, add_output_options=True)

    is_stdin = True
    if len(args) > 0:
        pysam_in = pysam.AlignmentFile(args[0], "rb")
        if args[0] != "-":
            is_stdin = False
    elif args.stdin == sys.stdin:
        pysam_in = pysam.AlignmentFile("-", "rb")
    else:
        pysam_in = pysam.AlignmentFile(args.stdin, "rb")
        if args.stdin != "-":
            is_stdin = False

    if args.input_filename_fasta:
        fasta = pysam.FastaFile(args.input_filename_fasta)
    else:
        fasta = None

    counts_df = bam2stats_window_count(
        pysam_in,
        region=args.region,
        window_size=args.window_size,
        fasta=fasta)

    if not args.output_all_windows:
        counts_df = counts_df[counts_df.alignments > 0]

    # add G+C content
    if fasta:
        counts_df["percent_gc"] = 100.0 * counts_df.bases_gc / (counts_df.bases_gc + counts_df.bases_at)
        counts_df.fillna(0, inplace=True)

    counts_df.to_csv(
        args.stdout,
        sep="\t")

    E.stop()
