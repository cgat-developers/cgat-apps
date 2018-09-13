"""split mapped reads in a bam files into shorter segments.
"""


import pysam
import sys
import tempfile
import copy
import cgatcore.iotools as iotools
import cgatcore.experiment as E
from cgat.BamTools.bamtools import bam2bam_split_reads


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-i", "--input-bam", dest="input_bam_file", type="string",
        help="input bam file [%default]")

    parser.add_option(
        "-o", "--output-bam", dest="output_bam_file", type="string",
        help="input bam file [%default].")

    parser.add_option(
        "-r", "--max-read-length", dest="max_read_length", type="int",
        help="maximum read length [%default].")

    parser.add_option(
        "-m", "--output-mode", dest="output_mode", type="choice",
        choices=["buffered", "direct"],
        help="output mode for files. 'buffered' will output reads in correct "
        "sort order, 'direct' will require the output BAM file to be sorted separately."
        "[%default].")

    parser.add_option(
        "--region", dest="region", type="string",
        help="genomic region, only split in BAM file within this region [%default].")

    parser.set_defaults(
        input_bam_file="-",
        output_bam_file="-",
        max_read_length=100,
        default_quality_score=10,
        region=None,
        output_mode="buffered",
    )

    (options, args) = E.start(parser, argv)

    pysam_in = pysam.Samfile(options.input_bam_file, "rb")
    pysam_out = pysam.Samfile(options.output_bam_file, "wb", template=pysam_in)

    max_read_length = options.max_read_length

    bam2bam_split_reads(pysam_in, pysam_out,
                        default_quality_score=options.default_quality_score,
                        max_read_length=options.max_read_length,
                        output_mode=options.output_mode)

    E.stop()

if __name__ == "__main__":
    sys.exit(main())
