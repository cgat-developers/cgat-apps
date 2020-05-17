'''
fastqs2fasta.py - interleave two fastq files
=============================================

:Tags: Genomics NGS FASTQ FASTA Conversion

Purpose
-------

This script is used to interleave two :term:`fastq`-formatted files
(paired data) into a single :term:`fasta`-formatted file. Read1 is
followed by read2 in the resultant file.

:term:`fastq` files MUST be sorted by read identifier.

Usage
-----

For example::

   cgat fastqs2fasta \
         --first-fastq-file=in.fastq.1.gz \
         --second-fastq-file=in.fastq.2.gz > out.fasta

If :file:`in.fastq.1.gz` looks like this::

    @r1_from_gi|387760314|ref|NC_017594.1|_Streptococcus_saliva_#0/1
    TTCTTGTTGAATCATTTCAATTGTCTCCTTTTAGTTTTATTAGATAATAACAGCTTCTTCCACAACTTCT
    +
    ??A???ABBDDBDDEDGGFGAFHHCHHIIIDIHGIFIH=HFICIHDHIHIFIFIIIIIIHFHIFHIHHHH
    @r3_from_gi|315441696|ref|NC_014814.1|_Mycobacterium_gilvum_#0/1
    ATGAACGCGGCCGAGCAACACCGCCACCACGTGAATCGGTGGTTCTACGACTGCCCGTCGGCCTTCCACC
    +

and :file:`in.fastq.2.gz` looks like this::

    A??A?B??BDBDDDBDGGFA>CFCFIIIIIIF;HFIGHCIGHIHHEHHHIIHHFDHH-HD-IDHHHGIHG
    @r1_from_gi|387760314|ref|NC_017594.1|_Streptococcus_saliva_#0/2
    ACCTTCGTTTCCAAGGTGCAGCAGGTCAACTTGATCAAACTGCCCCTTTGAACGAAGTGAAAAAACAAAT
    +
    A????@BBDBDDADABGFGFFEHHHIEHHII@IIHIHHIDHCCIHIIIHHIEI5HIHFHIEHIH=CHHC)
    @r3_from_gi|315441696|ref|NC_014814.1|_Mycobacterium_gilvum_#0/2
    GGGAGCCTGCAGCGCCGCCGCGACTGCATCGCCGCGGCCGGCATCGTGGGATGGACGGTGCGTCAGACGC
    +
    ???A?9BBDDD5@DDDGFFGFFHIIIHHIHBFHIIHIIHHH>HEIHHFI>FFHGIIHHHDHCCFIHFIHD

then the output will be::

  >r1_from_gi|387760314|ref|NC_017594.1|_Streptococcus_saliva_#0/1
  TTCTTGTTGAATCATTTCAATTGTCTCCTTTTAGTTTTATTAGATAATAACAGCTTCTTCCACAACTTCT
  >r1_from_gi|387760314|ref|NC_017594.1|_Streptococcus_saliva_#0/2
  ACCTTCGTTTCCAAGGTGCAGCAGGTCAACTTGATCAAACTGCCCCTTTGAACGAAGTGAAAAAACAAAT
  >r3_from_gi|315441696|ref|NC_014814.1|_Mycobacterium_gilvum_#0/1
  ATGAACGCGGCCGAGCAACACCGCCACCACGTGAATCGGTGGTTCTACGACTGCCCGTCGGCCTTCCACC
  >r3_from_gi|315441696|ref|NC_014814.1|_Mycobacterium_gilvum_#0/2
  GGGAGCCTGCAGCGCCGCCGCGACTGCATCGCCGCGGCCGGCATCGTGGGATGGACGGTGCGTCAGACGC
  >r4_from_gi|53711291|ref|NC_006347.1|_Bacteroides_fragilis_#0/1
  GAGGGATCAGCCTGTTATCCCCGGAGTACCTTTTATCCTTTGAGcgatGTCCCTTCCATACGGAAACACC
  >r4_from_gi|53711291|ref|NC_006347.1|_Bacteroides_fragilis_#0/2
  CAACCGTGAGCTCAGTGAAATTGTAGTATCGGTGAAGATGCcgatTACCCGcgatGGGACGAAAAGACCC
  >r5_from_gi|325297172|ref|NC_015164.1|_Bacteroides_salanitr_#0/1
  TGCGGCGAAATACCAGCCCATGCCCCGTCCCCAGAATTCCTTGGAGCAGCCTTTGTGAGGTTCGGCTTTG
  >r5_from_gi|325297172|ref|NC_015164.1|_Bacteroides_salanitr_#0/2
  AACGGCACGCACAATGCCGACCGCTACAAAAAGGCTGCCGACTGGCTCCGCAATTACCTGGTGAACGACT


Type::

   cgat fastqs2fasta --help

for command line help.


Command line options
--------------------

'''

import sys
from itertools import zip_longest

import cgatcore.iotools as iotools
import cgat.Fastq as Fastq
import cgatcore.experiment as E


class PairedReadError(Exception):

    '''
    exception raised when reads aren't paired -
    could be not sorted or files of different lengths
    '''


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument(
        "-a", "--first-fastq-file", dest="fastq1", type=str,
        help="supply read1 fastq file")
    parser.add_argument(
        "-b", "--second-fastq-file", dest="fastq2", type=str,
        help="supply read2 fastq file")

    # add common options (-h/--help, ...) and parse command line
    (args, unknown) = E.start(parser,
                              argv=argv,
                              unknowns=True)

    if unknown and len(unknown) == 2:
        args.fastq1, args.fastq2 = unknown

    fastq1 = iotools.open_file(args.fastq1)
    fastq2 = iotools.open_file(args.fastq2)

    E.info("iterating over fastq files")
    f1_count = 0
    for f1, f2 in zip_longest(Fastq.iterate(fastq1),
                              Fastq.iterate(fastq2)):
        if not (f1 and f2) or (not f2 and f1):
            try:
                raise PairedReadError(
                    "unpaired reads detected. Are files sorted? are "
                    "files of equal length?")
            except PairedReadError as e:
                raise PairedReadError(e).with_traceback(sys.exc_info()[2])
        else:
            assert f1.identifier.endswith("/1") and \
                f2.identifier.endswith("/2"), \
                "Reads in file 1 must end with /1 and reads in file 2 with /2"
            args.stdout.write(
                ">%s\n%s\n>%s\n%s\n" %
                (f1.identifier, f1.seq, f2.identifier, f2.seq))
            f1_count += 1

    E.info("output: %i pairs" % f1_count)

    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
