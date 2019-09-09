'''fasta2stats - extract statistics from fasta file
===================================================

This tool outputs the number of sequences and total length of
sequences. It works on uncompressed and compressed files and will make
use of a samtools faidx index if it exists.

Command line options
--------------------

'''
import os
import sys
import numpy
import pysam

import cgatcore.experiment as E
import cgatcore.iotools as iotools


def main(argv=None):

    parser = E.ArgumentParser(descriptin=__doc__)

    parser.add_argument(
        "-f", "--fasta", dest="input_filename_fasta",
        type=str,
        help="filename with fasta sequences. ")

    parser.add_argument(
        "-o", "--output-filename-sequences", dest="output_filename_sequences",
        type=str,
        help="output per sequence information to filename")

    parser.set_defaults(
        input_filename_fasta=None,
    )

    (args, unknown) = E.start(parser,
                              argv=argv,
                              unknowns=True)

    if len(unnowns) > 0:
        args.input_filename_fasta = args[0]

    sequence_pairs = []

    if args.input_filename_fasta != "-" and os.path.exists(
            args.input_filename_fasta + ".fai"):
        has_index = 1
        fastafile = pysam.FastaFile(args.input_filename_fasta)
        sequence_pairs = list(zip(fastafile.references, fastafile.lengths))
    else:
        has_index = 0
        iterator = pysam.FastxFile(args.input_filename_fasta)
        for record in iterator:
            sequence_pairs.append(
                (record.name,
                 len(record.sequence)))

    lengths = numpy.array([x[1] for x in sequence_pairs])

    args.stdout.write("\t".join((
        "has_index", "nsequences", "total_length", "min_length",
        "max_length", "median_length", "mean_length")) + "\n")

    if len(lengths) > 0:
        args.stdout.write("\t".join(map(str, (
            has_index,
            len(sequence_pairs),
            lengths.sum(),
            lengths.min(),
            lengths.max(),
            numpy.median(lengths),
            lengths.mean()))) + "\n")
    else:
        args.stdout.write("\t".join(map(str, (
            has_index,
            len(sequence_pairs),
            0,
            "",
            "",
            "",
            ""))) + "\n")

    if args.output_filename_sequences:
        with iotools.open_file(args.output_filename_sequences, "w") as outf:
            outf.write("name\tlength\n")
            outf.write(
                "\n".join(["\t".join(map(str, x)) for x in sequence_pairs]) + "\n")

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
