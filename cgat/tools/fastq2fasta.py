"""Convert fastq files to fasta files

Implemeted methods
==================

ont2pacbio
----------

Convert long-read fastq files, for example from ONT data,
to pacbio formatted fasta files that can be used with the
dazzler toolkit. Pacbio fasta files have the following
naming convention for the fastq header line::

    ><name>/<well>/<start>_<end> RQ=0.<qv>

where well, start, end and qv are integer values.

The quality value will be set to the average of the base qualities.

"""

import sys
import pysam
import math

import cgatcore.experiment as E


def main(argv=sys.argv):

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-i", "--input-fastq", dest="input_fastq_file", type=str,
        help="input fastq file")

    parser.add_argument(
        "-m", "--method", dest="method", type=str,
        choices=["ont2pacbio"],
        help="methods to apply ")

    parser.set_defaults(
        input_fastq_file=None,
        line_width=80,
        method=None,
    )

    (args, unknown) = E.start(parser,
                              argv,
                              add_output_options=True,
                              unknowns=True)

    if len(unknown) == 1:
        args.input_fastq_file = unknown[0]

    if args.input_fastq_file == "-":
        args.input_fastq_file = args.stdin

    outf = args.stdout
    line_width = args.line_width
    well_no = 0
    for record in pysam.FastqFile(args.input_fastq_file):
        well_no += 1
        quals = record.get_quality_array()
        seq = record.sequence
        qv = int(math.floor(sum(quals) / len(quals)))
        outf.write(">{}/{}/{}_{} RQ=0.{}\n".format(
            "test", well_no, 1, len(seq) + 1, qv))
        for x in range(0, len(seq), line_width):
            outf.write(seq[x:x + line_width] + "\n")

    E.stop()

if __name__ == "__main__":
    sys.exit(main())
