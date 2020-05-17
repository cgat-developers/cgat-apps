'''compute stats from a bam-file

Purpose
-------

This script takes a bam file as input and computes a few metrics by
iterating over the file. The metrics output are:

+------------------------+------------------------------------------+
|*Category*              |*Content*                                 |
+------------------------+------------------------------------------+
|total                   |total number of alignments in bam file    |
+------------------------+------------------------------------------+
|alignments_mapped       |alignments mapped to a chromosome (bam    |
|                        |flag)                                     |
+------------------------+------------------------------------------+
|alignments_unmapped     |alignments unmapped (bam flag)            |
+------------------------+------------------------------------------+
|qc_fail                 |alignments failing QC (bam flag)          |
+------------------------+------------------------------------------+
|mate_unmapped           |alignments in which the mate is unmapped  |
|                        |(bam flag)                                |
+------------------------+------------------------------------------+
|reverse                 |alignments in which read maps to reverse  |
|                        |strand (bam flag)                         |
+------------------------+------------------------------------------+
|mate_reverse            |alignments in which mate maps to reverse  |
|                        |strand (bam flag)                         |
+------------------------+------------------------------------------+
|proper_pair             |alignments in which both pairs have been  |
|                        |mapped properly (according to the mapper) |
|                        |(bam flag)                                |
+------------------------+------------------------------------------+
|read1                   |alignments for 1st read of pair (bam flag)|
+------------------------+------------------------------------------+
|paired                  |alignments of reads that are paired (bam  |
|                        |flag)                                     |
+------------------------+------------------------------------------+
|duplicate               |read is PCR or optical duplicate (bam     |
|                        |flag)                                     |
+------------------------+------------------------------------------+
|read2                   |alignment is for 2nd read of pair (bam    |
|                        |flag)                                     |
+------------------------+------------------------------------------+
|secondary               |alignment is not primary alignment        |
+------------------------+------------------------------------------+
|alignments_duplicates   |number of alignments mapping to the same  |
|                        |location                                  |
+------------------------+------------------------------------------+
|alignments_unique       |number of alignments mapping to unique    |
|                        |locations                                 |
+------------------------+------------------------------------------+
|reads_total             |number of reads in file. Either given via |
|                        |--num-reads or deduc ed as the sum of     |
|                        |mappend and unmapped reads                |
+------------------------+------------------------------------------+
|reads_mapped            |number of reads mapping in file. Derived  |
|                        |from the total number o f alignments and  |
|                        |removing counts for multiple              |
|                        |matches. Requires the NH flag to be set   |
|                        |correctly.                                |
+------------------------+------------------------------------------+
|reads_unmapped          |number of reads unmapped in file. Assumes |
|                        |that there is only one                    |
|                        |entry per unmapped read.                  |
+------------------------+------------------------------------------+
|reads_missing           |number of reads missing, if number of     |
|                        |reads given by --input-rea ds. Otherwise  |
|                        |0.                                        |
+------------------------+------------------------------------------+
|pairs_total             |number of total pairs - this is the number|
|                        |of reads_total divided by two. If there   |
|                        |were no pairs, pairs_total will be 0.     |
+------------------------+------------------------------------------+
|pairs_mapped            |number of mapped pairs - this is the same |
|                        |as the number of proper pairs.            |
+------------------------+------------------------------------------+

Additionally, the script outputs histograms for the following tags and
scores.

* NM: number of mismatches in alignments.
* NH: number of hits of reads.
* mapq: mapping quality of alignments.

Supplying a fastq file
++++++++++++++++++++++

If a fastq file is supplied (``--fastq-file``), the script will
compute some additional summary statistics. However, as it builds a dictionary
of all sequences, it will also require a good  amount of memory. The additional
metrics output are:

+-----------------------------+----------------------------------------+
|*Category*                   |*Content*                               |
+-----------------------------+----------------------------------------+
|pairs_total                  |total number of pairs in input data     |
+-----------------------------+----------------------------------------+
|pairs_mapped                 |pairs in which both reads map           |
+-----------------------------+----------------------------------------+
|pairs_unmapped               |pairs in which neither read maps        |
+-----------------------------+----------------------------------------+
|pairs_proper_unique          |pairs which are proper and map uniquely.|
+-----------------------------+----------------------------------------+
|pairs_incomplete_unique      |pairs in which one of the reads maps    |
|                             |uniquely, but the other does not map.   |
+-----------------------------+----------------------------------------+
|pairs_incomplete_multimapping|pairs in which one of the reads maps    |
|                             |uniquely, but the other maps to multiple|
|                             |locations.                              |
+-----------------------------+----------------------------------------+
|pairs_proper_duplicate       |pairs which are proper and unique, but  |
|                             |marked as duplicates.                   |
+-----------------------------+----------------------------------------+
|pairs_proper_multimapping    |pairs which are proper, but map to      |
|                             |multiple locations.                     |
+-----------------------------+----------------------------------------+
|pairs_not_proper_unique      |pairs mapping uniquely, but not flagged |
|                             |as proper                               |
+-----------------------------+----------------------------------------+
|pairs_other                  |pairs not in any of the above categories|
+-----------------------------+----------------------------------------+

Note that for paired-end data, any ``\1`` or ``/2`` suffixes will be
removed from the read name in the assumption that these have been removed
in the bam file as well.

Usage
-----

Example::

   python bam2stats.py in.bam

This command will generate various statistics based on the supplied
BAM file, such as percentage reads mapped and percentage reads mapped
in pairs. The output looks like this:

+-----------------------------+------+-------+-----------------+
|category                     |counts|percent|of               |
+-----------------------------+------+-------+-----------------+
|alignments_total             |32018 |100.00 |alignments_total |
+-----------------------------+------+-------+-----------------+
|alignments_mapped            |32018 |100.00 |alignments_total |
+-----------------------------+------+-------+-----------------+
|alignments_unmapped          |0     | 0.00  |alignments_total |
+-----------------------------+------+-------+-----------------+
|alignments_qc_fail           |0     | 0.00  |alignments_mapped|
+-----------------------------+------+-------+-----------------+
|alignments_mate_unmapped     |241   | 0.75  |alignments_mapped|
+-----------------------------+------+-------+-----------------+
|alignments_reverse           |16016 |50.02  |alignments_mapped|
+-----------------------------+------+-------+-----------------+
|alignments_mate_reverse      |15893 |49.64  |alignments_mapped|
+-----------------------------+------+-------+-----------------+
|alignments_proper_pair       |30865 |96.40  |alignments_mapped|
+-----------------------------+------+-------+-----------------+
|alignments_read1             |16057 |50.15  |alignments_mapped|
+-----------------------------+------+-------+-----------------+
|alignments_paired            |32018 |100.00 |alignments_mapped|
+-----------------------------+------+-------+-----------------+
|alignments_duplicate         |0     | 0.00  |alignments_mapped|
+-----------------------------+------+-------+-----------------+
|alignments_read2             |15961 |49.85  |alignments_mapped|
+-----------------------------+------+-------+-----------------+
|alignments_secondary         |0     | 0.00  |alignments_mapped|
+-----------------------------+------+-------+-----------------+
|alignments_filtered          |31950 |99.79  |alignments_mapped|
+-----------------------------+------+-------+-----------------+
|reads_total                  |34250 |100.00 |reads_total      |
+-----------------------------+------+-------+-----------------+
|reads_unmapped               |0     | 0.00  |reads_total      |
+-----------------------------+------+-------+-----------------+
|reads_mapped                 |32018 |93.48  |reads_total      |
+-----------------------------+------+-------+-----------------+
|reads_missing                |2232  | 6.52  |reads_total      |
+-----------------------------+------+-------+-----------------+
|reads_mapped_unique          |32018 |100.00 |reads_mapped     |
+-----------------------------+------+-------+-----------------+
|reads_multimapping           |0     | 0.00  |reads_mapped     |
+-----------------------------+------+-------+-----------------+
|pairs_total                  |17125 |100.00 |pairs_total      |
+-----------------------------+------+-------+-----------------+
|pairs_mapped                 |17125 |100.00 |pairs_total      |
+-----------------------------+------+-------+-----------------+
|pairs_unmapped               |0     | 0.00  |pairs_total      |
+-----------------------------+------+-------+-----------------+
|pairs_proper_unique          |14880 |86.89  |pairs_total      |
+-----------------------------+------+-------+-----------------+
|pairs_incomplete_unique      |2232  |13.03  |pairs_total      |
+-----------------------------+------+-------+-----------------+
|pairs_incomplete_multimapping|0     | 0.00  |pairs_total      |
+-----------------------------+------+-------+-----------------+
|pairs_proper_duplicate       |0     | 0.00  |pairs_total      |
+-----------------------------+------+-------+-----------------+
|pairs_proper_multimapping    |0     | 0.00  |pairs_total      |
+-----------------------------+------+-------+-----------------+
|pairs_not_proper_unique      |13    | 0.08  |pairs_total      |
+-----------------------------+------+-------+-----------------+
|pairs_other                  |0     | 0.00  |pairs_total      |
+-----------------------------+------+-------+-----------------+
|read1_total                  |17125 |100.00 |read1_total      |
+-----------------------------+------+-------+-----------------+
|read1_unmapped               |0     | 0.00  |read1_total      |
+-----------------------------+------+-------+-----------------+
|read1_mapped                 |16057 |93.76  |read1_total      |
+-----------------------------+------+-------+-----------------+
|read1_mapped_unique          |16057 |100.00 |read1_mapped     |
+-----------------------------+------+-------+-----------------+
|reads_multimapping           |0     | 0.00  |read1_mapped     |
+-----------------------------+------+-------+-----------------+
|read1_missing                |1068  | 6.65  |read1_total      |
+-----------------------------+------+-------+-----------------+
|read2_total                  |17125 |100.00 |read2_total      |
+-----------------------------+------+-------+-----------------+
|read2_unmapped               |0     | 0.00  |read2_total      |
+-----------------------------+------+-------+-----------------+
|read2_mapped                 |15961 |93.20  |read2_total      |
+-----------------------------+------+-------+-----------------+
|read2_mapped_unique          |15961 |100.00 |read2_mapped     |
+-----------------------------+------+-------+-----------------+
|reads_multimapping           |0     | 0.00  |read2_mapped     |
+-----------------------------+------+-------+-----------------+
|read2_missing                |1164  | 7.29  |read2_total      |
+-----------------------------+------+-------+-----------------+

The first column contains the caterogy, the second the number of
counts and the third a percentage. The fourth column denotes the
denomiminator that was used to compute the percentage. In the table
above, wee see that 16,057 first reads in a pair map and 15,961
second reads in pair map, resulting in 14,880 proper uniquely mapped
pairs.

Type::

   cgat bam2stats --help

for command line help.

Bam2stats can read from standard input::

   cat in.bam | python bam2stats.py -


Documentation
-------------

Reads are not counted via read name, but making use of NH and HI flags
when present.  To recap, NH is the number of reported alignments that
contain the query in the current record, while HI is the hit index and
ranges from 0 to NH-1.

Unfortunately, not all aligners follow this convention. For example,
gsnap seems to set NH to the number of reportable alignments, while
the actual number of reported alignments in the file is less. Thus, if
the HI flag is present, the maximum HI is used to correct the NH
flag. The assumption is, that the same reporting threshold has been
used for all alignments.

If no NH flag is present, it is assumed that all reads have only been
reported once.

Multi-matching counts after filtering are really guesswork. Basically,
the assumption is that filtering is consistent and will tend to remove
all alignments of a query.

The error rates are computed using the following key:

substitution_rate
   Number of mismatches divided by number of aligned bases. This is the same
   as the samtools stats error rate.
insertion_rate
   Number of deletions in the read/insertions in the reference divided by the
   number of aligned bases.
deletion_rate
   Number of insertions in the read/deletions in the reference divided by the
   number of aligned bases.
error_rate
   Number of mismatches and deletions in the read divided by the number of
   aligned bases.
coverage
   Percentage of bases aligned divided by read length.

The following graphic illustrates the computation. A `.` signifies a
position that is included in the metric with `X` being an error::

   AAAAACAAAA AAAAAAAA   Reference
    AAAAAAAAAAAA AAA     Read
    ....X.... ......     substitution_rate NM / (CMATCH + CINS)
    .........X.. ...     insertion_rate CINS / (CMATCH + CINS)
    ......... ..X...     deletion_rate CDEL / (CMATCH + CDEL)
    ....X.... ..X...     error_rate NM / (CMATCH + CINS) (corresponds to samtools stats)
    .........X.. ...     match_rate CMATCH / (CMATCH + CINS)
    ....X.... .. ...     mismatch_rate NM / (CMATCH) (1 - percent_identity/100)

With CINS: Insertion into the reference (consumes read, but not
reference) and CDEL=Deletion from the reference (consumes reference,
but not read).

Command line options
--------------------

'''

import os
import sys
import cgatcore.experiment as E
import cgatcore.iotools as iotools
import numpy
import pandas
import pysam
import collections

import cgat.GTF as GTF
from cgat.BamTools.bamtools import bam2stats_count

FLAGS = {
    1: 'paired',
    2: 'proper_pair',
    4: 'unmapped',
    8: 'mate_unmapped',
    16: 'reverse',
    32: 'mate_reverse',
    64: 'read1',
    128: 'read2',
    256: 'secondary',
    512: 'qc_fail',
    1024: 'duplicate',
    2048: 'supplementary',
}


def computeMappedReadsFromAlignments(total_alignments, nh, max_hi):
    '''compute number of reads alignment from total number of alignments.
    '''
    nreads_mapped = total_alignments
    if len(nh) > 0:
        max_nh = max(nh.keys())
        if max_hi > 0:
            for x in range(2, min(max_nh + 1, max_hi)):
                nreads_mapped -= (nh[x] / x) * (x - 1)
            for x in range(max_hi, max_nh + 1):
                nreads_mapped -= (nh[x] / max_hi) * (max_hi - 1)
        else:
            for x in range(2, max(nh.keys()) + 1):
                nreads_mapped -= (nh[x] / x) * (x - 1)

    return nreads_mapped


def writeNH(outfile, nh, max_hi):
    '''output nh array, correcting for max_hi if less than nh'''

    # need to remove double counting
    # one read matching to 2 positions is only 2

    max_nh = max(nh.keys())
    if max_hi > 0:
        for x in range(1, min(max_nh + 1, max_hi)):
            if nh[x] == 0:
                continue
            outfile.write("%i\t%i\n" % (x, nh[x] / x))
        for x in range(max_hi, max_nh + 1):
            if nh[x] == 0:
                continue
            outfile.write("%i\t%i\n" % (x, nh[x] / max_hi))
    else:
        for x in range(1, max_nh + 1):
            if nh[x] == 0:
                continue
            outfile.write("%i\t%i\n" % (x, nh[x] / x))


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
        "-r", "--mask-bed-file", "--mask-gff-file", dest="filename_bed", type=str,
        metavar='GFF',
        help="gff formatted file with masking locations. The number of "
        "reads overlapping the intervals in the given file will be "
        "computed. Note that the computation currently does not take "
        "into account indels, so it is an approximate count only. "
        )

    parser.add_argument(
        "-f", "--ignore-masked-reads", dest="ignore_masked_reads", action="store_true",
        help="as well as counting reads in the file given by --mask-bed-file, "
        "also remove these reads for duplicate and match statistics. "
        )

    parser.add_argument(
        "-i", "--num-reads", dest="input_reads", type=int,
        help="the number of reads - if given, used to provide percentages "
        )

    parser.add_argument(
        "-d", "--output-details", dest="output_details", action="store_true",
        help="output per-read details into a separate file. Read names are "
        "md5/base64 encoded ")

    parser.add_argument(
        "--output-readmap", dest="output_readmap", action="store_true",
        help="output map between read name and "
        "md5/base64 encoded short name")

    parser.add_argument(
        "--add-alignment-details", dest="add_alignment_details", action="store_true",
        help="add alignment details to per-read details. Implies --output-details "
        )

    parser.add_argument(
        "-q", "--fastq-file", dest="filename_fastq",
        help="filename with sequences and quality scores. This file is only "
        "used to collect sequence identifiers. Thus, for paired end data a "
        "single file is sufficient ")

    parser.add_argument(
        "--basic-counts", dest="detailed_count", action="store_false",
        help="perform basic counting and do not compute per read stats. "
        "This is more memory efficient and faster stats computation, "
        "but only a summary counts table is output ")

    parser.set_defaults(
        filename_bed=None,
        ignore_masked_reads=False,
        input_reads=0,
        force_output=False,
        filename_fastq=None,
        detailed_count=True,
        output_details=False,
        output_readmap=False,
        add_alignment_details=False,
    )

    # add common options (-h/--help, ...) and parse command line
    (args, unknown) = E.start(parser, argv=argv, add_output_options=True, unknowns=True)

    if args.filename_bed:
        bed_mask = GTF.readAndIndex(
            GTF.iterator(iotools.open_file(args.filename_bed)))
    else:
        bed_mask = None
    
    if args.add_alignment_details:
        args.output_details = True

    is_stdin = True
    if len(unknown) > 0:
        pysam_in = pysam.AlignmentFile(unknown[0], "rb")
        if unknown[0] != "-":
            is_stdin = False
    elif args.stdin == sys.stdin:
        pysam_in = pysam.AlignmentFile("-", "rb")
    else:
        pysam_in = pysam.AlignmentFile(args.stdin, "rb")
        if args.stdin != "-":
            is_stdin = False

    if args.output_details:
        outfile_details = E.open_output_file("details", "w")
    else:
        outfile_details = None

    if args.output_readmap:
        outfile_readmap = E.open_output_file("readmap", "w")
    else:
        outfile_readmap = None

    if args.filename_fastq and not os.path.exists(args.filename_fastq):
        raise IOError("file %s does not exist" % args.filename_fastq)

    (counter, flags_counts, nh_filtered, nh_all,
     nm_filtered, nm_all, mapq, mapq_all, max_hi, details_df) = \
        bam2stats_count(pysam_in,
                        bed_mask=bed_mask,
                        ignore_masked_reads=args.ignore_masked_reads,
                        is_stdin=is_stdin,
                        filename_fastq=args.filename_fastq,
                        outfile_details=outfile_details,
                        add_alignment_details=args.add_alignment_details,
                        outfile_readmap=outfile_readmap,
                        detailed_count=args.detailed_count)

    if max_hi > 0 and max_hi != max(nh_all.keys()):
        E.warn("max_hi(%i) is inconsistent with max_nh (%i) "
               "- counts will be corrected"
               % (max_hi, max(nh_all.keys())))

    outs = args.stdout
    outs.write("category\tcounts\tpercent\tof\n")

    def _write(outs, text, numerator, denominator, base):
        percent = iotools.pretty_percent(numerator, denominator)
        outs.write('%s\t%i\t%s\t%s\n' % (text,
                                         numerator,
                                         percent,
                                         base))

    ###############################
    ###############################
    ###############################
    # Output alignment information
    ###############################
    nalignments_unmapped = flags_counts["unmapped"]
    nalignments_mapped = counter.alignments_input - nalignments_unmapped

    _write(outs,
           "alignments_total",
           counter.alignments_input,
           counter.alignments_input,
           "alignments_total")

    if counter.alignments_input == 0:
        E.warn("no alignments in BAM file - no further output")
        E.stop()
        return

    _write(outs,
           "alignments_mapped",
           nalignments_mapped,
           counter.alignments_input,
           'alignments_total')
    _write(outs,
           "alignments_unmapped",
           nalignments_unmapped,
           counter.alignments_input,
           'alignments_total')

    if nalignments_mapped == 0:
        E.warn("no mapped alignments - no further output")
        E.stop()
        return

    for flag, counts in sorted(flags_counts.items()):
        if flag == "unmapped":
            continue
        _write(outs,
               'alignments_' + flag,
               counts,
               nalignments_mapped,
               'alignments_mapped')

    if args.filename_bed:
        _write(outs,
               "alignments_masked",
               counter.alignments_masked,
               nalignments_mapped,
               'alignments_mapped')
        _write(outs,
               "alignments_notmasked",
               counter.alignments_notmasked,
               nalignments_mapped,
               'alignments_mapped')

    _write(outs,
           "alignments_filtered",
           counter.alignments_filtered,
           nalignments_mapped,
           "alignments_mapped")

    if counter.filtered == nalignments_mapped:
        normby = "alignments_mapped"
    else:
        normby = "alignments_filtered"

    if counter.filtered > 0:
        _write(outs,
               "alignments_duplicates",
               counter.alignments_duplicates,
               counter.alignments_filtered,
               normby)
        _write(outs,
               "alignments_unique",
               counter.aligmnments_filtered - counter.alignments_duplicates,
               counter.alignments_filtered,
               normby)

    ###############################
    ###############################
    ###############################
    # Output read based information
    ###############################

    # derive the number of mapped reads in file from alignment counts
    if args.filename_fastq or not is_stdin:
        nreads_total = counter.total_read
        _write(outs,
               "reads_total",
               counter.total_read,
               nreads_total,
               'reads_total')
        _write(outs,
               "reads_unmapped",
               counter.total_read_is_unmapped,
               nreads_total,
               'reads_total')
        _write(outs,
               "reads_mapped",
               counter.total_read_is_mapped,
               nreads_total,
               'reads_total')
        _write(outs,
               "reads_missing",
               counter.total_read_is_missing,
               nreads_total,
               'reads_total')
        _write(outs,
               "reads_mapped_unique",
               counter.total_read_is_mapped_uniq,
               counter.total_read_is_mapped,
               'reads_mapped')
        _write(outs,
               "reads_multimapping",
               counter.total_read_is_mmap,
               counter.total_read_is_mapped,
               'reads_mapped')
        _write(outs,
               "reads_mapped_supplementary",
               counter.total_read_has_supplementary,
               counter.total_read_is_mapped,
               'reads_mapped')
    else:
        E.warn('inferring read counts from alignments and NH tags')
        nreads_unmapped = flags_counts["unmapped"]
        nreads_mapped = computeMappedReadsFromAlignments(nalignments_mapped,
                                                         nh_all, max_hi)

        nreads_missing = 0
        if args.input_reads:
            nreads_total = args.input_reads
            # unmapped reads in bam file?
            if nreads_unmapped:
                nreads_missing = nreads_total - nreads_unmapped - nreads_mapped
            else:
                nreads_unmapped = nreads_total - nreads_mapped

        elif nreads_unmapped:
            # if unmapped reads are in bam file, take those
            nreads_total = nreads_mapped + nreads_unmapped
        else:
            # otherwise normalize by mapped reads
            nreads_unmapped = 0
            nreads_total = nreads_mapped

        outs.write("reads_total\t%i\t%5.2f\treads_total\n" %
                   (nreads_total, 100.0))
        outs.write("reads_mapped\t%i\t%5.2f\treads_total\n" %
                   (nreads_mapped, 100.0 * nreads_mapped / nreads_total))
        outs.write("reads_unmapped\t%i\t%5.2f\treads_total\n" %
                   (nreads_unmapped, 100.0 * nreads_unmapped / nreads_total))
        outs.write("reads_missing\t%i\t%5.2f\treads_total\n" %
                   (nreads_missing, 100.0 * nreads_missing / nreads_total))

        if len(nh_all) > 1:
            outs.write("reads_unique\t%i\t%5.2f\treads_mapped\n" %
                       (nh_all[1], 100.0 * nh_all[1] / nreads_mapped))

    pysam_in.close()

    ###############################
    ###############################
    ###############################
    # Output pair information
    ###############################
    if flags_counts["read2"] > 0:
        if args.filename_fastq:
            pairs_mapped = counter.total_pair_is_mapped

            # sanity check
            assert counter.total_pair_is_mapped == \
                (counter.total_pair_is_proper_uniq +
                 counter.total_pair_is_incomplete_uniq +
                 counter.total_pair_is_incomplete_mmap +
                 counter.total_pair_is_proper_duplicate +
                 counter.total_pair_is_proper_mmap +
                 counter.total_pair_not_proper_uniq +
                 counter.total_pair_is_other)

            outs.write("pairs_total\t%i\t%5.2f\tpairs_total\n" %
                       (counter.total_pairs,
                        100.0 * counter.total_pairs / counter.total_pairs))
            outs.write("pairs_mapped\t%i\t%5.2f\tpairs_total\n" %
                       (pairs_mapped,
                        100.0 * pairs_mapped / counter.total_pairs))
            outs.write(
                "pairs_unmapped\t%i\t%5.2f\tpairs_total\n" %
                (counter.total_pair_is_unmapped,
                 100.0 * counter.total_pair_is_unmapped / counter.total_pairs))
            outs.write(
                "pairs_proper_unique\t%i\t%5.2f\tpairs_total\n" %
                (counter.total_pair_is_proper_uniq,
                 100.0 * counter.total_pair_is_proper_uniq /
                 counter.total_pairs))
            outs.write(
                "pairs_incomplete_unique\t%i\t%5.2f\tpairs_total\n" %
                (counter.total_pair_is_incomplete_uniq,
                 100.0 * counter.total_pair_is_incomplete_uniq /
                 counter.total_pairs))
            outs.write(
                "pairs_incomplete_multimapping\t%i\t%5.2f\tpairs_total\n" %
                (counter.total_pair_is_incomplete_mmap,
                 100.0 * counter.total_pair_is_incomplete_mmap /
                 counter.total_pairs))
            outs.write(
                "pairs_proper_duplicate\t%i\t%5.2f\tpairs_total\n" %
                (counter.total_pair_is_proper_duplicate,
                 100.0 * counter.total_pair_is_proper_duplicate /
                 counter.total_pairs))
            outs.write(
                "pairs_proper_multimapping\t%i\t%5.2f\tpairs_total\n" %
                (counter.total_pair_is_proper_mmap,
                 100.0 * counter.total_pair_is_proper_mmap /
                 counter.total_pairs))
            outs.write(
                "pairs_not_proper_unique\t%i\t%5.2f\tpairs_total\n" %
                (counter.total_pair_not_proper_uniq,
                 100.0 * counter.total_pair_not_proper_uniq /
                 counter.total_pairs))
            outs.write(
                "pairs_other\t%i\t%5.2f\tpairs_total\n" %
                (counter.total_pair_is_other,
                 100.0 * counter.total_pair_is_other /
                 counter.total_pairs))

            nread1_total = counter.total_read1
            _write(outs,
                   "read1_total",
                   counter.total_read1,
                   nread1_total,
                   'read1_total')
            _write(outs,
                   "read1_unmapped",
                   counter.total_read1_is_unmapped,
                   nread1_total,
                   'read1_total')
            _write(outs,
                   "read1_mapped",
                   counter.total_read1_is_mapped,
                   nread1_total,
                   'read1_total')
            _write(outs,
                   "read1_mapped_unique",
                   counter.total_read1_is_mapped_uniq,
                   counter.total_read1_is_mapped,
                   'read1_mapped')
            _write(outs,
                   "reads_multimapping",
                   counter.total_read1_is_mmap,
                   counter.total_read1_is_mapped,
                   'read1_mapped')
            _write(outs,
                   "read1_missing",
                   counter.total_read1_is_missing,
                   counter.total_read1_is_mapped,
                   'read1_total')

            nread2_total = counter.total_read2
            _write(outs,
                   "read2_total",
                   counter.total_read2,
                   nread2_total,
                   'read2_total')
            _write(outs,
                   "read2_unmapped",
                   counter.total_read2_is_unmapped,
                   nread2_total,
                   'read2_total')
            _write(outs,
                   "read2_mapped",
                   counter.total_read2_is_mapped,
                   nread2_total,
                   'read2_total')
            _write(outs,
                   "read2_mapped_unique",
                   counter.total_read2_is_mapped_uniq,
                   counter.total_read2_is_mapped,
                   'read2_mapped')
            _write(outs,
                   "reads_multimapping",
                   counter.total_read2_is_mmap,
                   counter.total_read2_is_mapped,
                   'read2_mapped')
            _write(outs,
                   "read2_missing",
                   counter.total_read2_is_missing,
                   counter.total_read2_is_mapped,
                   'read2_total')

        else:
            # approximate counts
            pairs_total = nreads_total // 2
            pairs_mapped = flags_counts["proper_pair"] // 2
            _write(outs,
                   "pairs_total",
                   pairs_total,
                   pairs_total,
                   "pairs_total")
            _write(outs,
                   "pairs_mapped",
                   pairs_mapped,
                   pairs_total,
                   "pairs_total")
    else:
        # no paired end data
        pairs_total = pairs_mapped = 0
        outs.write("pairs_total\t%i\t%5.2f\tpairs_total\n" %
                   (pairs_total, 0.0))
        outs.write("pairs_mapped\t%i\t%5.2f\tpairs_total\n" %
                   (pairs_mapped, 0.0))

    outs.write("error_rate\t%i\t%5.2f\tmatches+insertions\n" %
               (counter.error_counts, counter.error_rate * 100.0))
    outs.write("insertion_rate\t%i\t%5.2f\tmatches+insertions\n" %
               (counter.insertion_counts, counter.insertion_rate * 100.0))
    outs.write("deletion_rate\t%i\t%5.2f\tmatches+deletions\n" %
               (counter.deletion_counts, counter.deletion_rate * 100.0))
    outs.write("mismatch_rate\t%i\t%5.2f\tmatches\n" %
               (counter.mismatch_counts, counter.mismatch_rate * 100.0))
    outs.write("match_rate\t%i\t%5.2f\tmatches+insertions\n" %
               (counter.match_counts, counter.match_rate * 100.0))

    if args.force_output or len(nm_filtered) > 0:
        outfile = E.open_output_file("nm", "w")
        outfile.write("NM\talignments\n")
        if len(nm_filtered) > 0:
            for x in range(0, max(nm_filtered.keys()) + 1):
                outfile.write("%i\t%i\n" % (x, nm_filtered[x]))
        else:
            outfile.write("0\t%i\n" % (counter.filtered))
        outfile.close()

    if args.force_output or len(nh_all) > 1:
        outfile = E.open_output_file("nh_all", "w")
        outfile.write("NH\treads\n")
        if len(nh_all) > 0:
            writeNH(outfile, nh_all, max_hi)
        else:
            # assume all are unique if NH flag not set
            outfile.write("1\t%i\n" % (counter.mapped_reads))
        outfile.close()

    if args.force_output or len(nh_filtered) > 1:
        outfile = E.open_output_file("nh", "w")
        outfile.write("NH\treads\n")
        if len(nh_filtered) > 0:
            writeNH(outfile, nh_filtered, max_hi)
        else:
            # assume all are unique if NH flag not set
            outfile.write("1\t%i\n" % (counter.filtered))
        outfile.close()

    if args.force_output or len(mapq_all) > 1:
        outfile = E.open_output_file("mapq", "w")
        outfile.write("mapq\tall_reads\tfiltered_reads\n")
        for x in range(0, max(mapq_all.keys()) + 1):
            outfile.write("%i\t%i\t%i\n" % (x, mapq_all[x], mapq[x]))
        outfile.close()

    if details_df is not None:
        with E.open_output_file("summaries", "w") as outf:
            details_df.describe().transpose().to_csv(
                outf, sep="\t", index_label="metric")
        bins = numpy.arange(0, 1.01, 0.01)
        histogram_df = pandas.DataFrame.from_dict(collections.OrderedDict(
            [(x, numpy.histogram(details_df[x].dropna(),
                                 bins=bins)[0]) for x in details_df.columns]))

        histogram_df.index = numpy.arange(0, 1.0, 0.01)

        row_sums = histogram_df.sum(axis=1)
        histogram_df = histogram_df[row_sums != 0]

        with E.open_output_file("histogram", "w") as outf:
            histogram_df.to_csv(outf, sep="\t", index_label="bin")

    # write footer and output benchmark information.
    E.stop()
