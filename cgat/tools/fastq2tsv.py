"""modify a fastq file.
"""

import collections
import sys
import pysam
import cgatcore.experiment as E
import cgatcore.iotools as iotools


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-i", "--input-fastq-file", dest="input_fastq_file", type="string",
        help="input fastq file. "
        "[%default]")

    parser.add_option(
        "-m", "--method", dest="methods", action="append", type="choice",
        choices=("length", ),
        help="methods to apply [%default]")

    parser.set_defaults(
        methods=[],
        input_fastq_file=None,
    )

    (options, args) = E.start(parser, argv)

    if len(args) == 1:
        options.input_fastq_file = args[0]

    if options.input_fastq_file is None:
        raise ValueError("missing input fastq file")

    counter = E.Counter()

    # note: complete rewrite with Counters, currently only length
    if options.methods != ["length"]:
        raise NotImplementedError()

    with pysam.FastqFile(options.input_fastq_file) as inf:

        for read in inf:
            counter.input += 1
            options.stdout.write("\t".join(
                map(str, (read.name, len(read.sequence)))) + "\n")

            counter.output += 1

    E.info(counter)
    E.stop()
