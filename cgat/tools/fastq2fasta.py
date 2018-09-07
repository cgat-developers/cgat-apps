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

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-i", "--input-fastq", dest="input_fastq_file", type="string",
        help="input fastq file")

    parser.add_option(
        "-m", "--method", dest="method", type="choice",
        choices=["ont2pacbio"],
        help="methods to apply [%default]")

    parser.set_defaults(
        input_fastq_file=None,
        line_width=80,
        method=None,
    )

    (options, args) = E.start(parser, argv, add_output_options=True)

    if len(args) == 1:
        options.input_fastq_file = args[0]

    if options.input_fastq_file == "-":
        options.input_fastq_file = options.stdin

    outf = options.stdout
    line_width = options.line_width
    well_no = 0
    for record in pysam.FastqFile(options.input_fastq_file):
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
