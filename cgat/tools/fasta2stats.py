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

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-f", "--fasta", dest="input_filename_fasta",
        type="string",
        help="filename with fasta sequences. ")

    parser.add_option(
        "-o", "--output-filename-sequences", dest="output_filename_sequences",
        type="string",
        help="output per sequence information to filename")

    parser.set_defaults(
        input_filename_fasta=None,
    )

    (options, args) = E.start(parser, argv=argv)

    if len(args) > 0:
        options.input_filename_fasta = args[0]

    sequence_pairs = []

    if options.input_filename_fasta != "-" and os.path.exists(
            options.input_filename_fasta + ".fai"):
        has_index = 1
        fastafile = pysam.FastaFile(options.input_filename_fasta)
        sequence_pairs = list(zip(fastafile.references, fastafile.lengths))
    else:
        has_index = 0
        iterator = pysam.FastxFile(options.input_filename_fasta)
        for record in iterator:
            sequence_pairs.append(
                (record.name,
                 len(record.sequence)))

    lengths = numpy.array([x[1] for x in sequence_pairs])

    options.stdout.write("\t".join((
        "has_index", "nsequences", "total_length", "min_length",
        "max_length", "median_length", "mean_length")) + "\n")

    if len(lengths) > 0:
        options.stdout.write("\t".join(map(str, (
            has_index,
            len(sequence_pairs),
            lengths.sum(),
            lengths.min(),
            lengths.max(),
            numpy.median(lengths),
            lengths.mean()))) + "\n")
    else:
        options.stdout.write("\t".join(map(str, (
            has_index,
            len(sequence_pairs),
            0,
            "",
            "",
            "",
            ""))) + "\n")

    if options.output_filename_sequences:
        with iotools.open_file(options.output_filename_sequences, "w") as outf:
            outf.write("name\tlength\n")
            outf.write(
                "\n".join(["\t".join(map(str, x)) for x in sequence_pairs]) + "\n")

    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
