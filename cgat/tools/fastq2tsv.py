"""modify a fastq file.
"""

import collections
import sys
import pysam
import cgatcore.experiment as E
import cgatcore.iotools as iotools


def main(argv=sys.argv):

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-i", "--input-fastq-file", dest="input_fastq_file", type=str,
        help="input fastq file. "
        )

    parser.add_argument(
        "-m", "--method", dest="methods", action="append", type=str,
        choices=("length", ),
        help="methods to apply ")

    parser.set_defaults(
        methods=[],
        input_fastq_file=None,
    )

    (args, unknown) = E.start(parser,
                              argv,
                              unknowns=True)

    if len(unknown) == 1:
        args.input_fastq_file = unknown[0]

    if args.input_fastq_file is None:
        raise ValueError("missing input fastq file")

    counter = E.Counter()

    # note: complete rewrite with Counters, currently only length
    if args.methods != ["length"]:
        raise NotImplementedError()

    with pysam.FastqFile(args.input_fastq_file) as inf:

        for read in inf:
            counter.input += 1
            args.stdout.write("\t".join(
                map(str, (read.name, len(read.sequence)))) + "\n")

            counter.output += 1

    E.info(counter)
    E.stop()
