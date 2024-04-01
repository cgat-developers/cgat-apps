"""
BamTools - Utilities for working with BAM files
===============================================

This module brings together convenience function for working
with :term:`bam` formatted files.

"""

from pysam.libchtslib cimport *
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from pysam.libcalignedsegment cimport pysam_bam_get_cigar, \
    pysam_bam_get_qname, pysam_get_n_cigar
from pysam.libcfaidx cimport *
from libc.string cimport strchr
from libc.stdint cimport int8_t
from libc.stdio cimport puts, printf
from libc.stdlib cimport abs
from cpython cimport PyErr_SetString, PyBytes_FromStringAndSize
from cpython cimport array as c_array
from sortedcontainers import SortedList

import array
import base64
import collections
import copy
import struct
import hashlib
import itertools
import numpy
import pandas
import pysam
import sys

cimport numpy

import cgatcore.experiment as E

def parse_region_string(s):
    """parse a genomic region string.

    Returns tuple of contig, start, end. Missing values are None.
    """
    if not s:
        return None, None, None
    if ":" in s:
        contig, coords = s.split(":")
        if "-" in coords:
            start, end = list(map(int, coords.split("-")))
        else:
            start = int(coords)
            end = None
        return contig, start, end
    else:
        return s, None, None


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
    2048: 'supplementary'}

cdef struct CountsType:
    int32_t read_length
    int8_t alignments
    int8_t is_mapped
    int8_t is_unmapped
    int8_t mate_is_unmapped
    int8_t is_paired
    int8_t mapped_is_read1
    int8_t mapped_is_read2
    int8_t unmapped_is_read1
    int8_t unmapped_is_read2
    int8_t is_proper_pair
    int8_t is_secondary
    int8_t is_qcfail
    int8_t is_duplicate
    int8_t is_supplementary
    int32_t mean_quality
    int32_t median_quality

cdef int NM = 10

# pass views and fill instead of creating.
# This saves time, but also prevents possible memory leaks.
cdef int NCIGAR_CODES = 10
cdef get_cigar_stats(AlignedSegment read, uint32_t[:] base_counts, uint32_t[:] block_counts):
    cdef int nfields = NCIGAR_CODES + 1

    cdef bam1_t * src = read._delegate
    cdef int op
    cdef uint32_t l
    cdef int32_t k
    cdef uint32_t * cigar_p = pysam_bam_get_cigar(src)

    for k from 0 <= k < nfields:
        base_counts[k] = 0
        block_counts[k] = 0

    if cigar_p == NULL:
        return None

    for k from 0 <= k < pysam_get_n_cigar(src):
        op = cigar_p[k] & BAM_CIGAR_MASK
        l = cigar_p[k] >> BAM_CIGAR_SHIFT
        base_counts[op] += l
        block_counts[op] += 1

    cdef uint8_t * v = bam_aux_get(src, 'NM')
    if v != NULL:
        base_counts[nfields - 1] = <int32_t>bam_aux2i(v)

    return base_counts, block_counts


def bam2stats_count(AlignmentFile samfile,
                    bed_mask=None,
                    ignore_masked_reads=False,
                    is_stdin=True,
                    filename_fastq=None,
                    outfile_details=None,
                    add_alignment_details=False,
                    outfile_readmap=None,
                    detailed_count=None):
    '''
    '''
    cdef AlignedSegment read
    cdef bint _add_alignment_details = add_alignment_details
    # counters
    cdef int ninput = 0
    cdef int nduplicates = 0
    # number of reads present after filtering
    cdef int nfiltered = 0
    # number of reads overlapping masked regions (if bed_mask != None)
    cdef int nmasked = 0
    # number of reads not overlap RNA (if bed_mask != None)
    cdef int nnotmasked = 0

    cdef bint _ignore_masked_reads = ignore_masked_reads
    cdef int max_hi = 0
    # count nh, nm tags
    nh_filtered = collections.defaultdict(int)
    nm_filtered = collections.defaultdict(int)
    nh_all = collections.defaultdict(int)
    nm_all = collections.defaultdict(int)
    mapq_filtered = collections.defaultdict(int)
    mapq_all = collections.defaultdict(int)

    cdef int * flags_counts = <int*>calloc(len(FLAGS), sizeof(int))

    # helper variables
    cdef int last_tid = -1
    cdef int last_pos = 0
    cdef int32_t read_index

    duplicates = collections.defaultdict(int)
    counts = collections.defaultdict(int)
    cdef int count = 0

    cdef int tid = -1
    cdef int flag
    cdef uint32_t iteration = 0
    cdef uint32_t report_step = 1000000
    cdef uint32_t nalignments = 0
    cdef uint8_t * v
    cdef int32_t nm
    cdef int32_t nh
    cdef int32_t hi
    cdef int x
    cdef int lflags = len(FLAGS)
    cdef uint32_t f

    # Todo: make this filter configurable
    cdef int detail_filter_flags = 2304

    # detailed counting
    cdef FastxRecord fq
    cdef int64_t index, fastq_nreads
    cdef CountsType * fastq_counts
    cdef CountsType * fastq_count
    cdef char * read_name
    cdef char * position
    cdef bint count_fastq = False
    cdef int fastq_notfound = 0
    cdef int chop = 0
    cdef int hash_size = 5

    # number of alignment pairs that are nucleotide mismatches
    cdef long mismatch_counts = 0
    # number of alignment pairs where there is an insertion in the reference
    cdef long insertion_counts = 0
    # number of alignment pairs where there is a deletion in the reference
    cdef long deletion_counts = 0
    # number of alignment pairs that are in match state
    # (might be a nucleotide mismatch)
    cdef long match_counts = 0

    cdef uint16_t [:, :] alignment_details_view

    cdef int nfields = NCIGAR_CODES + 1

    cdef c_array.array base_counts = array.array(
        "I",
         [0] * nfields)
    cdef c_array.array block_counts = array.array(
             "I",
        [0] * nfields)

    cdef uint32_t [:] base_counts_view = base_counts
    cdef uint32_t [:] block_counts_view = block_counts

    if filename_fastq != None:
        count_fastq = True
        E.info("reading fastq file")
        # Using a python dictionary here is due
        # for a large amount of memory usage.
        # Alternatives to dictionary
        # 1. POSIX hash tables (hsearch,...) or trees: very slow
        # 2. custom hash implementation: worth the effort?
        # 3. Sorted list and binary search: too slow for many lookups
        reads = {}
        fastqfile = FastxFile(filename_fastq)
        fastq_nreads = 0
        for fq in fastqfile:
            if fastq_nreads == 0:
                # chop off /1 or /2 as mappers usually remove these
                # suffices. Test only the first.
                name = fq.name
                if name.endswith("/1") or name.endswith("/2"):
                    chop = -2

            if chop != 0:
                name = hashlib.md5(fq.name[:chop].encode("ascii")).digest()[:hash_size]
            else:
                name = hashlib.md5(fq.name.encode("ascii")).digest()[:hash_size]

            reads[name] = fastq_nreads

            fastq_nreads += 1

    elif not is_stdin and detailed_count:
        count_fastq = True
        old_pos = samfile.tell()
        E.info("two-pass processing to collect read names")
        reads = {}
        fastq_nreads = 0
        for iteration, read in enumerate(samfile):

            if iteration % report_step == 0:
                E.info("pre-pass: read {} alignments: {} reads".format(
                    iteration, len(reads)))

            read_name = pysam_bam_get_qname(read._delegate)

            # terminate string at first space to
            # truncate read names containing more than
            # just the id
            position = strchr(read_name, ' ')
            if position != NULL:
                position[0] = '\0'

            # compute hash, truncate to half the md5 digest which
            # should be 64 bits and provides 36,893,488,147,419,103,232
            # hashes.
            md5 = hashlib.md5(read_name).digest()[:hash_size]
            if md5 not in reads:
                reads[md5] = fastq_nreads
                fastq_nreads += 1
                if outfile_readmap:
                    outfile_readmap.write("{}\t{}\n".format(
                        read_name,
                        base64.encodebytes(md5)[:-1]))

        samfile.seek(old_pos)
        nalignments = iteration
    else:
        count_fastq = False
        E.info("simple counting only")

    if count_fastq:
        E.info("read names of %i reads or read pairs" % fastq_nreads)
        E.info("allocating %i bytes" % (fastq_nreads * sizeof(CountsType)))

        fastq_counts = <CountsType *>calloc(fastq_nreads, sizeof(CountsType))
        if fastq_counts == NULL:
            raise ValueError("could not allocate memory for %i bytes" %
                             (fastq_nreads * sizeof(CountsType)))

    if _add_alignment_details:
        alignment_details = numpy.zeros((fastq_nreads, NCIGAR_CODES + 1), dtype=numpy.uint16)

    E.info("starting processing of alignment file")

    for iteration, read in enumerate(samfile):

        if iteration % report_step == 0:
            if nalignments:
                E.info("read {}/{} alignments: {}%".format(
                    iteration, nalignments,
                    100.0 * iteration / nalignments))
            else:
                E.info("read {} alignments".format(
                    iteration))

        flag = read._delegate.core.flag
        ninput += 1

        f = 1
        for x from 0 <= x < lflags:
            if flag & f:
                flags_counts[x] += 1
            f = f << 1

        # get maximum NI field
        v = bam_aux_get(read._delegate, 'HI')
        if v != NULL:
            hi = <int32_t>bam_aux2i(v)
            if hi > max_hi: max_hi = hi

        v = bam_aux_get(read._delegate, 'NH')
        if v != NULL:
            nh = <int32_t>bam_aux2i(v)
            nh_all[nh] += 1
        else:
            nh = -1

        v = bam_aux_get(read._delegate, 'NM')
        if v != NULL:
            nm = <int32_t>bam_aux2i(v)
            nm_all[nm] += 1
        else:
            nm = -1

        mapq_all[read.mapq] += 1

        get_cigar_stats(read, base_counts_view, block_counts_view)

        if not read.is_secondary:
            # for consistency with samtools
            mismatch_counts += nm
            deletion_counts += base_counts_view[BAM_CDEL]
            insertion_counts += base_counts_view[BAM_CINS]
            match_counts += base_counts_view[BAM_CMATCH]

        if count_fastq:
            read_name = pysam_bam_get_qname(read._delegate)
            # terminate string at first space to
            # truncate read names containing more than
            # just the id
            position = strchr(read_name, ' ')
            if position != NULL:
                position[0] = '\0'

            md5 = hashlib.md5(read_name).digest()[:hash_size]
            read_index = reads[md5]
            try:
                fastq_count = &fastq_counts[read_index]
            except KeyError:
                fastq_notfound += 1
                continue

            # only take primary alignments for read length. read.query_length
            # includes soft-clipped sequence.
            if not read.is_secondary and not read.is_supplementary:
                fastq_count.read_length = read.query_length
                q = read.query_qualities
                # multiply by 1000 to increase significant digits
                if q:
                    fastq_count.mean_quality = round(numpy.mean(q) * 1000.0)
                    fastq_count.median_quality = round(numpy.median(q) * 1000.0)
                else:
                    fastq_count.mean_quality = 0
                    fastq_count.median_quality = 0

            # book-keeping counts
            fastq_count.alignments += 1

            if read.is_qcfail: fastq_count.is_qcfail += 1
            if read.is_duplicate: fastq_count.is_duplicate += 1
            if read.is_paired: fastq_count.is_paired += 1

            if read.is_unmapped:
                fastq_count.is_unmapped += 1
                if read.is_read1: fastq_count.unmapped_is_read1 += 1
                if read.is_read2: fastq_count.unmapped_is_read2 += 1
            else:
                # only count is_read1 and is_read2 for mapped
                fastq_count.is_mapped += 1
                if read.mate_is_unmapped: fastq_count.mate_is_unmapped += 1
                if read.is_read1: fastq_count.mapped_is_read1 += 1
                if read.is_read2: fastq_count.mapped_is_read2 += 1
                if read.is_proper_pair: fastq_count.is_proper_pair += 1
                if read.is_secondary:
                    fastq_count.is_secondary += 1
                if read.is_supplementary:
                    fastq_count.is_supplementary += 1
                if _add_alignment_details:
                    if not (flag & detail_filter_flags):
                        alignment_details[read_index, :] += base_counts_view

        # is paired and read2
        #if flag & 1 and flag & 128:
            # ignore second reads to avoid double counting pairs
            # this needs to be done properly by counting in
            # pairs
        #    continue

        # skip unmapped reads
        if read._delegate.core.flag & 4:
            continue

        if read.tid != last_tid:
            contig = samfile.getrname(read.rname)

        # note: does not take into account gaps within reads
        # or partial overlap.
        if bed_mask:
            if bed_mask.contains(contig, read.pos, read.pos + read.alen):
                nmasked += 1
                if _ignore_masked_reads:
                    continue
            else:
                nnotmasked += 1

        nfiltered += 1

        if nh >= 0: nh_filtered[nh] += 1
        if nm >= 0: nm_filtered[nm] += 1
        mapq_filtered[read.mapq] += 1

        # duplicate analysis - simply count per start position
        # ignoring sequence and strand
        if read.tid == last_tid and read.pos == last_pos:
            count += 1
            nduplicates += 1
            continue

        if count > 1:
            counts[count] += 1

        count = 1
        last_tid, last_pos = read.tid, read.pos

    E.info( "finished computing counts" )

    if fastq_notfound:
        E.warn( "could not match %i records in bam-file to fastq file" % fastq_notfound)

    counter = E.Counter()

    counter.alignments_input = ninput
    counter.alignments_filtered = nfiltered
    counter.alignments_duplicates = nduplicates
    counter.alignments_masked = nmasked
    counter.alignments_notmasked = nnotmasked

    if match_counts == 0:
        match_counts = 1

    counter.error_counts = mismatch_counts + insertion_counts
    counter.match_counts = match_counts
    counter.mismatch_counts = mismatch_counts
    counter.deletion_counts = deletion_counts
    counter.insertion_counts = insertion_counts

    counter.error_rate = float(mismatch_counts) / (match_counts + insertion_counts)
    counter.match_rate = float(match_counts) / (match_counts + insertion_counts)
    counter.mismatch_rate = float(mismatch_counts) / (match_counts)
    counter.deletion_rate = float(deletion_counts) / (match_counts + deletion_counts)
    counter.insertion_rate = float(insertion_counts) / (match_counts + insertion_counts)

    # convert flags to labels
    t = {}
    f = 1
    for x from 0 <= x < lflags:
        t[FLAGS[f]] = flags_counts[x]
        f = f << 1

    # count based on fastq data
    cdef int total_paired = 0
    cdef int total_unpaired = 0
    cdef int total_pair_is_mapped = 0
    cdef int total_pair_is_unmapped = 0
    cdef int total_pair_is_proper_uniq = 0
    cdef int total_pair_is_proper_mmap = 0
    cdef int total_pair_is_proper_duplicate = 0
    cdef int total_pair_not_proper_uniq = 0
    cdef int total_pair_is_incomplete_uniq = 0
    cdef int total_pair_is_incomplete_mmap = 0
    cdef int total_pair_is_other = 0
    # read based counting for unpaired reads
    cdef int total_reads = 0
    cdef int total_read_is_mapped_uniquely = 0
    cdef int total_read_is_mapped_multiply = 0
    cdef int total_read_is_unmapped = 0
    cdef int total_read_is_missing = 0
    cdef int total_read_has_supplementary = 0
    # read based counting for read1
    cdef int total_read1 = 0
    cdef int total_read1_is_mapped_uniquely = 0
    cdef int total_read1_is_mapped_multiply = 0
    cdef int total_read1_is_unmapped = 0
    cdef int total_read1_is_missing = 0
    # read based counting for read2
    cdef int total_read2 = 0
    cdef int total_read2_is_mapped_uniquely = 0
    cdef int total_read2_is_mapped_multiply = 0
    cdef int total_read2_is_unmapped = 0
    cdef int total_read2_is_missing = 0

    cdef float substitution_rate
    cdef float insertion_rate
    cdef float deletion_rate
    cdef float error_rate
    cdef float coverage

    if count_fastq:
        E.info("fastq counting: aggregating counts")
        for qname, index in reads.items():
            fastq_count = &fastq_counts[index]

            # paired read data
            if fastq_count.is_paired:
                total_paired += 1

                if fastq_count.is_unmapped == fastq_count.is_paired:
                    # an unmapped read pair
                    total_pair_is_unmapped += 1
                    continue

                total_pair_is_mapped += 1

                if fastq_count.is_proper_pair == 2:
                    # a unique proper pair
                    total_pair_is_proper_uniq +=1
                    # a duplicate unique proper pair
                    if fastq_count.is_duplicate == 2:
                        total_pair_is_proper_duplicate +=1
                elif fastq_count.is_proper_pair > 2:
                    # proper pairs that map to multiple locations
                    total_pair_is_proper_mmap +=1
                elif fastq_count.mapped_is_read1 == 1 \
                        and fastq_count.mapped_is_read2 == 1 \
                        and fastq_count.is_proper_pair == 0:
                    # pair which map each read once, but not is not proper
                    total_pair_not_proper_uniq += 1
                elif (fastq_count.mapped_is_read1 == 1 and
                      fastq_count.mapped_is_read2 == 0) \
                    or (fastq_count.mapped_is_read1 == 0 and
                        fastq_count.mapped_is_read2 == 1):
                    # an incomplete pair - one read of a pair matches uniquely
                    # but not the other
                    total_pair_is_incomplete_uniq += 1
                elif (fastq_count.mapped_is_read1 == 1 and
                      fastq_count.mapped_is_read2 > 1) \
                    or (fastq_count.mapped_is_read1 > 1 and
                        fastq_count.mapped_is_read2 == 1):
                    # an incomplete pair - one read of a pair matches uniquely
                    # but the other matches multiple times
                    total_pair_is_incomplete_mmap += 1
                else:
                    total_pair_is_other += 1

                # paired read counting
                if fastq_count.mapped_is_read1:
                    if fastq_count.mapped_is_read1 == 1:
                        total_read1_is_mapped_uniquely += 1
                    elif fastq_count.mapped_is_read1 > 1:
                        total_read1_is_mapped_multiply += 1
                elif fastq_count.unmapped_is_read1 > 0:
                    total_read1_is_unmapped += 1
                else:
                    total_read1_is_missing += 1

                if fastq_count.mapped_is_read2:
                    if fastq_count.mapped_is_read2 == 1:
                        total_read2_is_mapped_uniquely += 1
                    elif fastq_count.mapped_is_read2 > 1:
                        total_read2_is_mapped_multiply += 1
                elif fastq_count.unmapped_is_read2 > 0:
                    total_read2_is_unmapped += 1
                else:
                    total_read2_is_missing += 1

            else:
                # reads without data, could be unpaired or missing
                total_unpaired += 1

                if fastq_count.is_supplementary > 0:
                    total_read_has_supplementary += 1

                # counting for unpaired data
                if fastq_count.is_unmapped > 0:
                    total_read_is_unmapped += 1
                else:
                    if fastq_count.is_mapped == 1:
                        total_read_is_mapped_uniquely += 1
                    elif fastq_count.is_mapped > 1:
                        total_read_is_mapped_multiply += 1
                    else:
                        total_read_is_missing += 1

        ###################################################
        # Pair statistics
        counter.total_pairs = total_paired + total_unpaired
        counter.total_pair_is_mapped = total_pair_is_mapped
        counter.total_pair_is_unmapped = total_pair_is_unmapped
        counter.total_pair_is_proper_uniq = total_pair_is_proper_uniq
        counter.total_pair_is_incomplete_uniq = total_pair_is_incomplete_uniq
        counter.total_pair_is_incomplete_mmap = total_pair_is_incomplete_mmap
        counter.total_pair_is_proper_duplicate = total_pair_is_proper_duplicate
        counter.total_pair_is_proper_mmap = total_pair_is_proper_mmap
        counter.total_pair_not_proper_uniq = total_pair_not_proper_uniq
        counter.total_pair_is_other = total_pair_is_other

        ###################################################
        # Read Stats for SE and PE data
        # for SE data
        if total_paired == 0:
            E.info("reads: single end counting")
            counter.total_read = fastq_nreads
            counter.total_read_is_mapped_uniq = total_read_is_mapped_uniquely
            counter.total_read_is_mmap = total_read_is_mapped_multiply
            counter.total_read_is_mapped = total_read_is_mapped_multiply +\
                                           total_read_is_mapped_uniquely
            counter.total_read_is_unmapped = total_read_is_unmapped
            counter.total_read_is_missing = total_read_is_missing
            counter.total_read_has_supplementary = total_read_has_supplementary
        else:
            # for PE data
            E.info("reads: paired end counting: paired = %i" % total_paired)
            counter.total_read = fastq_nreads * 2
            counter.total_read_is_mapped_uniq = total_read1_is_mapped_uniquely +\
                                                total_read2_is_mapped_uniquely
            counter.total_read_is_mmap = total_read1_is_mapped_multiply +\
                                         total_read2_is_mapped_multiply
            counter.total_read_is_mapped = counter.total_read_is_mapped_uniq +\
                                           counter.total_read_is_mmap
            counter.total_read_is_unmapped = total_read1_is_unmapped +\
                                             total_read2_is_unmapped
            counter.total_read_is_missing = total_read1_is_missing +\
                                            total_read2_is_missing
            counter.total_read_has_supplementary = total_read_has_supplementary * 2

        ###################################################
        ## stats for 1st/2nd read separately
        counter.total_read1 = fastq_nreads
        counter.total_read1_is_mapped_uniq = total_read1_is_mapped_uniquely
        counter.total_read1_is_mmap = total_read1_is_mapped_multiply
        counter.total_read1_is_mapped = counter.total_read1_is_mapped_uniq + counter.total_read1_is_mmap
        counter.total_read1_is_unmapped = total_read1_is_unmapped
        counter.total_read1_is_missing = total_read1_is_missing
        counter.total_read2 = fastq_nreads
        counter.total_read2_is_mapped_uniq = total_read2_is_mapped_uniquely
        counter.total_read2_is_mmap = total_read2_is_mapped_multiply
        counter.total_read2_is_mapped = counter.total_read2_is_mapped_uniq + counter.total_read2_is_mmap
        counter.total_read2_is_unmapped = total_read2_is_unmapped
        counter.total_read2_is_missing = total_read2_is_missing

        substitution_rates = numpy.zeros(fastq_nreads, dtype=numpy.float64)
        insertion_rates = numpy.zeros(fastq_nreads, dtype=numpy.float64)
        deletion_rates = numpy.zeros(fastq_nreads, dtype=numpy.float64)
        error_rates = numpy.zeros(fastq_nreads, dtype=numpy.float64)
        coverages = numpy.zeros(fastq_nreads, dtype=numpy.float64)
        mask = numpy.ones(fastq_nreads, dtype=numpy.int64)

        if outfile_details:
            header = ["read_md5",
                      "read_length",
                      "alignments",
                      "mean_quality",
                      "median_quality",
                      "is_unmapped",
                      "mate_is_unmapped",
                      "is_paired",
                      "mapped_is_read1",
                      "mapped_is_read2",
                      "is_proper_pair",
                      "is_secondary",
                      "is_qcfail",
                      "is_duplicate",
                      "is_supplementary"]

            if _add_alignment_details:
                header.extend([
                    "cigar_match", "cigar_ins",
                    "cigar_del", "cigar_ref_skip",
                    "cigar_soft_clip", "cigar_hard_clip",
                    "cigar_pad", "cigar_equal",
                    "cigar_diff", "cigar_back",
                    "mismatches",
                    "coverage",
                    "substitution_rate",
                    "insertion_rate",
                    "deletion_rate",
                    "error_rate"])

            if outfile_details != sys.stdout:
                # later: get access FILE * object
                outfile_details.write("\t".join(header) + "\n")
                for qname, read_index in reads.items():
                    fastq_count = &fastq_counts[read_index]
                    # remove "\n" from base64 encoded md5
                    outfile_details.write("%s\t%s" % (
                        base64.encodebytes(qname)[:-1].decode("ascii"),
                        "\t".join( \
                                   map(str,
                                       (fastq_count.read_length,
                                        fastq_count.alignments,
                                        "{:5.2f}".format(fastq_count.mean_quality / 1000.0),
                                        "{:5.2f}".format(fastq_count.median_quality / 1000.0),
                                        fastq_count.is_unmapped,
                                        fastq_count.mate_is_unmapped,
                                        fastq_count.is_paired,
                                        fastq_count.mapped_is_read1,
                                        fastq_count.mapped_is_read2,
                                        fastq_count.is_proper_pair,
                                        fastq_count.is_secondary,
                                        fastq_count.is_qcfail,
                                        fastq_count.is_duplicate,
                                        fastq_count.is_supplementary
                                    ) )) ))
                    if _add_alignment_details:
                        if fastq_count.is_unmapped > 0:
                            substitution_rate = 0
                            insertion_rate = 0
                            deletion_rate = 0
                            error_rate = 0
                            coverage = 0
                            mask[read_index] = 0
                        else:
                            substitution_rate = float(alignment_details[read_index, NM]) /\
                                (alignment_details[read_index, BAM_CMATCH] +
                                 alignment_details[read_index, BAM_CINS])
                            insertion_rate = float(alignment_details[read_index, BAM_CINS]) /\
                                (alignment_details[read_index, BAM_CMATCH] +
                                 alignment_details[read_index, BAM_CINS])
                            deletion_rate = float(alignment_details[read_index, BAM_CDEL]) /\
                                            (alignment_details[read_index, BAM_CMATCH] +
                                             alignment_details[read_index, BAM_CDEL])
                            error_rate = float(alignment_details[read_index, NM] + \
                                          alignment_details[read_index, BAM_CINS]) /\
                                (alignment_details[read_index, BAM_CMATCH] +
                                 alignment_details[read_index, BAM_CINS])

                            # TODO: take into account supplementary aligments and soft/hard-masking?
                            if fastq_count.read_length > 0:
                                coverage = float(alignment_details[read_index, BAM_CMATCH]) /\
                                           fastq_count.read_length
                            else:
                                coverage = float(alignment_details[read_index, BAM_CMATCH]) /\
                                            (alignment_details[read_index, BAM_CMATCH] +
                                             alignment_details[read_index, BAM_CDEL])

                            # ignore reads > 32k as counters might have wrapped over
                            if fastq_count.read_length >= 32767:
                                mask[read_index] = 0

                            substitution_rates[read_index] = substitution_rate
                            insertion_rates[read_index] = insertion_rate
                            deletion_rates[read_index] = deletion_rate
                            error_rates[read_index] = error_rate
                            coverages[read_index] = coverage

                        outfile_details.write(
                            "\t" +
                            "\t".join(map(str, alignment_details[read_index, :])) +
                            "\t{:6.4f}".format(coverage) +
                            "\t{:6.4f}".format(substitution_rate) +
                            "\t{:6.4f}".format(insertion_rate) +
                            "\t{:6.4f}".format(deletion_rate) +
                            "\t{:6.4f}".format(error_rate) +
                            "\n")
                    else:
                        outfile_details.write("\n")
            else:
                # TODO: this code is incomplete.
                # output to stdout much quicker
                # use puts to avoid the following error:
                # format not a string literal and no format arguments
                puts("\t".join(header) + "\n")
                for qname, read_index in reads.items():
                    fastq_count = &fastq_counts[read_index]
                    read_name = qname
                    printf("%s\t%i\t%i\t%5.2f\t%5.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i",
                           read_name,
                           fastq_count.read_length,
                           fastq_count.alignments,
                           fastq_count.mean_quality / 1000.0,
                           fastq_count.median_quality / 1000.0,
                           fastq_count.is_unmapped,
                           fastq_count.mate_is_unmapped,
                           fastq_count.is_paired,
                           fastq_count.mapped_is_read1,
                           fastq_count.mapped_is_read2,
                           fastq_count.is_proper_pair,
                           fastq_count.is_secondary,
                           fastq_count.is_qcfail,
                           fastq_count.is_duplicate,
                           fastq_count.is_supplementary)
                    if _add_alignment_details:
                        printf("\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",
                               alignment_details_view[read_index, 0],
                               alignment_details_view[read_index, 1],
                               alignment_details_view[read_index, 2],
                               alignment_details_view[read_index, 3],
                               alignment_details_view[read_index, 4],
                               alignment_details_view[read_index, 5],
                               alignment_details_view[read_index, 6],
                               alignment_details_view[read_index, 7],
                               alignment_details_view[read_index, 8],
                               alignment_details_view[read_index, 9],
                               alignment_details_view[read_index, 10])

                    else:
                        printf("\n")

        details_df = pandas.DataFrame.from_dict(
                 collections.OrderedDict(zip(["substitution_rate",
                 "insertion_rate",
                 "deletion_rate",
                 "error_rate",
                 "coverage"],
                [substitution_rates,
                 insertion_rates,
                 deletion_rates,
                 error_rates,
                 coverages])))
        # subset to only take mapped reads
        details_df = details_df[mask == 1]
    else:
        details_df = None

    return (counter, t,
            nh_filtered,
            nh_all,
            nm_filtered,
            nm_all,
            mapq_filtered,
            mapq_all,
            max_hi,
            details_df)


class BufferedBAMOutput(object):
    """receive a list of reads and output reads in
    sorted order."""

    def __init__(self, samfile):
        self.samfile = samfile
        self.buffer = SortedList(key = lambda x: x.reference_start)

    def write(self, readlist):
        cdef AlignedSegment read
        for read in readlist:
            self.samfile.write(read)

    def add(self, read):
        self.buffer.add(read)

    def flush(self, read):

        cdef uint32_t pos
        if len(self.buffer) == 0:
            return

        # chromosome change
        if read.reference_id != self.buffer[0].reference_id:
            self.write(self.buffer)
            self.buffer = SortedList(key = lambda x: x.reference_start)
        else:
            to_output = []
            pos = read.reference_start
            while self.buffer and self.buffer[0].reference_start <= pos:
                to_output.append(self.buffer.pop(0))
            # E.debug("flushing buffer: {}, {} left".format(len(to_output),
            #                                               len(self.buffer)))
            self.write(to_output)

    def __dealloc__(self):
        self.write(self.buffer)

    def __len__(self):
        return len(self.buffer)


class DirectBAMOutput(object):
    """receive a list of reads and output reads as they arrive.
    """

    def __init__(self, samfile):
        self.samfile = samfile

    def add(self, AlignedSegment read):
        self.samfile.write(read)

    def flush(self, read):
        return

    def __len__(self):
        return 0


def bam2bam_split_reads(AlignmentFile pysam_in,
                        AlignmentFile pysam_out,
                        int max_read_length,
                        int default_quality_score,
                        output_mode="buffered"):

    if output_mode == "buffered":
        output_writer = BufferedBAMOutput(pysam_out)
    elif output_mode == "direct":
        output_writer = DirectBAMOutput(pysam_out)
    else:
        raise ValueError("unknown output mode: {}".format(output_mode))

    cdef AlignedSegment read
    cdef AlignedSegment out_read
    cdef int l, length, state
    cdef int remainder = 0
    cdef int left_bit = 0
    cdef int total_length = 0
    cdef int sequence_pos = 0
    cdef int match_length
    cdef int num_blocks
    cdef long current_reference_pos
    cdef long ninput = 0
    cdef long nskipped_no_sequence = 0
    cdef long noutput = 0

    for read in pysam_in.fetch():

        ninput += 1
        output_writer.flush(read)
        full_cigar = read.cigar
        match_length = sum([y for x, y in full_cigar if x == 0 or x == 7 or x == 8])
        num_blocks = match_length // max_read_length

        if ninput % 1000 == 0:
            E.debug("splitting read at {}:{} with {} matches into {} blocks: cache={}".format(
                    read.reference_name,
                    read.reference_start,
                    match_length,
                    num_blocks,
                    len(output_writer)))

        full_sequence = read.query_sequence
        if full_sequence is None:
            nskipped_no_sequence += 1
            continue

        full_qualities = read.query_qualities
        if full_qualities is None:
            full_qualities = [default_quality_score] * len(full_sequence)

        read.set_tag("MD", None)

        current_cigar = []
        current_sequence = []
        current_qualities = []
        current_reference_pos  = read.reference_start
        total_length = 0
        sequence_pos = 0
        for state, length in full_cigar:
            # print state, length, total_length, max_read_length
            # print current_cigar, current_sequence, current_qualities
            if state == 0 or state == 7 or state == 8:
                # test for overshoot
                if total_length + length >= max_read_length:
                    remainder = length
                    while remainder > 0:
                        # segment length to output:
                        left_bit = max_read_length - total_length

                        # add leftover segment
                        if left_bit > 0:
                            current_cigar.append((0, left_bit))
                            current_sequence.append(
                                full_sequence[sequence_pos:sequence_pos + left_bit])
                            current_qualities.extend(
                                full_qualities[sequence_pos:sequence_pos + left_bit])
                            remainder -= left_bit
                            sequence_pos += left_bit
                            current_reference_pos += left_bit

                        # output read
                        out_read = copy.copy(read)
                        out_read.cigartuples = current_cigar
                        out_read.query_sequence = "".join(current_sequence)
                        out_read.query_qualities = current_qualities
                        noutput += 1
                        output_writer.add(out_read)

                        # add remainder to read, but at most max read_length
                        l = min(remainder, max_read_length)
                        current_cigar = [(0, l)]
                        current_sequence = [
                            full_sequence[sequence_pos:sequence_pos + l]]
                        current_qualities = list(
                            full_qualities[sequence_pos:sequence_pos + l])
                        read.reference_start = current_reference_pos
                        sequence_pos += l
                        current_reference_pos += l

                        # subtract from remainder
                        total_length = l
                        remainder -= l
                    continue

                total_length += length

            if state in (0, 1, 4, 7, 8):
                current_sequence.append(
                    full_sequence[sequence_pos:sequence_pos + length])
                current_qualities.extend(
                    full_qualities[sequence_pos:sequence_pos + length])
                sequence_pos += length
            if state in (0, 2, 7, 8):
                current_reference_pos += length

            current_cigar.append((state, length))
        else:
            # emit last bit of read
            out_read = copy.copy(read)
            out_read.cigartuples = current_cigar
            out_read.query_sequence = "".join(current_sequence)
            out_read.query_qualities = current_qualities
            noutput += 1
            output_writer.add(out_read)

    E.info("splitting completed: ninput={}, noutput={}, nskipped_no_sequence={}".format(
        ninput, noutput, nskipped_no_sequence))

cdef class SetNH:

    cdef iter
    cdef stack

    def __cinit__(self, iter ):
        self.iter = itertools.groupby(iter, lambda x: x.qname)
        self.stack = []

    def __iter__(self):
        return self

    def __next__(self):
        """python version of next().
        """

        while 1:
            if self.stack:
                return self.stack.pop(0)
            else:
                key, x = next(self.iter)
                self.stack = list(x)
                nh = len(self.stack)
                for read in self.stack:
                    if not read.is_unmapped:
                        # deal with paired end reads counted
                        # as multi-mapping
                        if read.is_proper_pair and nh > 1:
                            nh -= 1
                        read.set_tag("NH", nh)


def bam2bam_filter_bam(AlignmentFile input_samfile,
                       AlignmentFile output_samfile,
                       AlignmentFile reference_samfile,
                       remove_nonunique=False,
                       remove_unique=False,
                       remove_contigs=None,
                       remove_unmapped=False,
                       remove_mismatches=False,
                       filter_mismatches=None,
                       filter_error_rate=None,
                       colour_mismatches=False,
                       int minimum_read_length=0,
                       float minimum_average_base_quality=0):

    '''To conserve memory, the tid and NM flag from *transcripts_samfile*
    are packed into memory. As a consequence, this method requires
    max(NM) < 256 (2^8) and max(tid) < 16777216 (2^24)

    If *remove_nonunique* is set, only uniquely matching reads will be
    output.

    If *remove_contigs* is given, contigs that are in remove_contigs
    will be removed. Note that this will only remove the alignment,
    but not all alignments of a particuluar read and the NH flag will
    *NOT* be updated.

    If *remove_unmapped* is given, unmapped reads will be removed.

    If *remove_mismatches* is set, only reads with number of
    mismatches better than in reference_samfile will be kept.

    If *colour_mismatches* is set, the ``CM`` tag will be used to
    count differences. By default, the ``NM`` tag is used.  The tag
    that is used needs to present in both *input_samfile* and
    *reference_samfile*.

    Detecting non-unique matches:

    This method first checks for the NH flag - if set, a unique match
    should have at most NH=1 hits.

    If not set, the method checks for BWA flags. Currently it checks
    if X0 is set (X0=Number of best hits found by BWA). Other relevant
    flags but currently not tested:

    * X1 = Number of suboptimal hits found by BWA
    * XT = Type: Unique/Repeat/N/Mate-sw
    '''

    cdef int ninput = 0
    cdef int nunmapped = 0
    cdef int nmismatches = 0
    cdef int noutput = 0
    cdef int nerror_rate = 0
    cdef int nnonunique = 0
    cdef int nunique = 0
    cdef int nremoved_contigs = 0
    cdef int nminimum_read_length = 0
    cdef int nminimum_average_base_quality = 0

    cdef bint c_remove_nonunique = remove_nonunique
    cdef bint c_remove_unique = remove_unique
    cdef bint c_remove_mismatches = remove_mismatches
    cdef bint c_remove_unmapped = remove_unmapped

    cdef int * remove_contig_tids
    cdef int nremove_contig_tids = 0
    cdef AlignedSegment read
    cdef AlignedSegment match

    cdef int32_t c_filter_error_rate = 0
    cdef int32_t error_rate = 0
    if filter_error_rate is not None:
        c_filter_error_rate = filter_error_rate * 10000

    # build index
    # this method will start indexing from the current file position
    # if you decide
    cdef int ret = 1
    cdef int x
    cdef bam1_t * b = <bam1_t*>calloc(1, sizeof( bam1_t))
    cdef uint64_t pos
    cdef uint8_t * v
    cdef int32_t nm
    cdef int32_t nh
    cdef int tid
    cdef int transript_tid
    cdef long val
    cdef bint skip = 0
    cdef float total_qualities

    cdef int32_t read_mismatches = 0

    # decide which tag to use
    cdef char * nm_tag = 'NM'
    cdef char * cm_tag = 'CM'
    cdef char * tag

    if colour_mismatches:
        tag = cm_tag
    else:
        tag = nm_tag

    if c_remove_mismatches:
        E.info( "building index" )
        if not reference_samfile:
            raise ValueError("require another bam file for mismatch filtering" )

        # L = 1 byte (unsigned char)
        def _gen(): return array.array('B')
        index = collections.defaultdict(_gen)

        while ret > 0:
            ret = bam_read1(hts_get_bgzfp(reference_samfile.htsfile),
                            b)
            if ret > 0:
                # ignore unmapped reads
                if b.core.flag & 4:
                    continue
                qname = pysam_bam_get_qname(b)
                tid = b.core.tid
                v = bam_aux_get(b, tag)
                if v != NULL:
                    nm = <int32_t>bam_aux2i(v)
                else:
                    nm = 0
                index[qname].append(nm)

        E.info( "built index for %i reads" % len(index))
        bam_destroy1(b)

    # setup list of contigs to remove:
    if remove_contigs:
        nremove_contig_tids = len(remove_contigs)
        remove_contig_tids = <int*>malloc(sizeof(int) * nremove_contig_tids)
        for x, rname in enumerate(remove_contigs):
            remove_contig_tids[x] = input_samfile.gettid(rname)

    E.info( "starting filtering" )

    cdef int nfields = NCIGAR_CODES + 1

    cdef c_array.array base_counts = array.array(
        "I",
         [0] * nfields)
    cdef c_array.array block_counts = array.array(
             "I",
        [0] * nfields)

    cdef uint32_t [:] base_counts_view = base_counts
    cdef uint32_t [:] block_counts_view = block_counts

    for read in input_samfile:

        ninput += 1
        # if ninput > 10000: break

        # remove unmapped reads
        if read._delegate.core.flag & 4:
            if c_remove_unmapped:
                nunmapped += 1
            else:
                output_samfile.write(read)
            continue

        if c_filter_error_rate > 0:
            v = bam_aux_get(read._delegate, 'NM')
            if v != NULL:
                nm = <int32_t>bam_aux2i(v)
            else:
                nm = 0

            get_cigar_stats(read, base_counts_view, block_counts_view)
            error_rate = (10000 * (nm + base_counts_view[BAM_CINS])) //\
                (base_counts_view[BAM_CMATCH] + base_counts_view[BAM_CINS])
            if error_rate > c_filter_error_rate:
                nerror_rate += 1
                continue

        # remove non-unique alignments
        if c_remove_nonunique:
            # check either NH or X0 (bwa) flag
            v = bam_aux_get(read._delegate, 'NH')
            if v == NULL:
                v = bam_aux_get(read._delegate, 'X0')

            if v != NULL:
                nh = <int32_t>bam_aux2i(v)
                if nh > 1:
                    nnonunique += 1
                    continue
        # remove unique alignments
        elif c_remove_unique:
            # check either NH or X0 (bwa) flag
            v = bam_aux_get(read._delegate, 'NH')
            if v == NULL:
                v = bam_aux_get(read._delegate, 'X0')

            if v != NULL:
                nh = <int32_t>bam_aux2i(v)
                if nh == 1:
                    nunique += 1
                    continue

        # remove reads matched to certain contigs
        if nremove_contig_tids:
            skip = 0
            for x from 0 <= x < nremove_contig_tids:
                if remove_contig_tids[x] == read.tid:
                    skip = 1
                    break
            if skip:
                nremoved_contigs += 1
                continue

        # remove reads in other bam file if their are better matching
        if c_remove_mismatches:
            if read.qname in index:
                # can compress index before, depends on
                # how many reads you expect test to filter
                nm = min(index[read.qname])
                read_mismatches = read.opt(tag)
                if nm > read_mismatches:
                    nmismatches += 1
                    continue

        if minimum_read_length and read.query_length < minimum_read_length:
            nminimum_read_length += 1
            continue

        if minimum_average_base_quality > 0:
            if read.query_qualities:
                total_qualities = sum(read.query_qualities)
                if total_qualities / read.query_length < minimum_average_base_quality:
                    nminimum_average_base_quality += 1
                    continue

        noutput += 1
        output_samfile.write(read)

    c = E.Counter()
    c.input = ninput
    c.removed_unique = nunique
    c.removed_nonunique = nnonunique
    c.removed_contigs = nremoved_contigs
    c.output = noutput
    c.removed_unmapped = nunmapped
    c.removed_mismatches = nmismatches
    c.error_rate = nerror_rate
    c.minimum_read_length = nminimum_read_length
    c.minimum_average_base_quality = nminimum_average_base_quality

    if nremove_contig_tids:
        free(remove_contig_tids)

    E.info("filtering finished")

    return c


cdef inline uint32_t get_alignment_length(bam1_t * src):
    cdef int k = 0
    cdef uint32_t l = 0
    if src == NULL:
        return 0
    cdef uint32_t * cigar_p = bam_get_cigar(src)
    if cigar_p == NULL:
        return 0
    cdef int op
    cdef int n = pysam_get_n_cigar(src)
    for k from 0 <= k < n:
        op = cigar_p[k] & BAM_CIGAR_MASK
        if op == BAM_CSOFT_CLIP or op == BAM_CHARD_CLIP:
            continue
        l += cigar_p[k] >> BAM_CIGAR_SHIFT
    return l


cdef bytes count_md_tag_mismatches(AlignedSegment read):
    """return byte string with match/mismatch information.

    matches: =
    mismatches: *
    insertions into reference: +
    deletions in reference: -


    """

    cdef bam1_t * src = read._delegate

    if src == NULL:
        return None

    cdef uint32_t * cigar_p = pysam_bam_get_cigar(src)
    if cigar_p == NULL:
        return None

    cdef uint32_t r_idx = 0
    cdef int op
    cdef uint32_t k, i, l, x
    cdef int nmatches = 0
    cdef int s_idx = 0

    cdef uint32_t max_len = get_alignment_length(src)
    if max_len == 0:
        raise ValueError("could not determine alignment length")

    cdef char * s = <char*>calloc(max_len + 1, sizeof(char))
    if s == NULL:
        raise ValueError(
            "could not allocated sequence of length %i" % max_len)

    for k from 0 <= k < pysam_get_n_cigar(src):
        op = cigar_p[k] & BAM_CIGAR_MASK
        l = cigar_p[k] >> BAM_CIGAR_SHIFT
        if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
            for i from 0 <= i < l:
                s[s_idx] = '='
                r_idx += 1
                s_idx += 1
        elif op == BAM_CDEL:
            for i from 0 <= i < l:
                s[s_idx] = '-'
                s_idx += 1
        elif op == BAM_CREF_SKIP:
            pass
        elif op == BAM_CINS:
            for i from 0 <= i < l:
                # encode insertions into reference as lowercase
                s[s_idx] = '+'
                r_idx += 1
                s_idx += 1
        elif op == BAM_CSOFT_CLIP:
            pass
        elif op == BAM_CHARD_CLIP:
            pass # advances neither
        elif op == BAM_CPAD:
            raise NotImplementedError(
                "Padding (BAM_CPAD, 6) is currently not supported. "
                "Please implement. Sorry about that.")

    cdef uint8_t * md_tag_ptr = bam_aux_get(src, "MD")
    if md_tag_ptr == NULL:
        seq = PyBytes_FromStringAndSize(s, s_idx)
        free(s)
        return seq

    cdef char * md_tag = <char*>bam_aux2Z(md_tag_ptr)
    cdef int md_idx = 0
    s_idx = 0

    while md_tag[md_idx] != 0:
        # c is numerical
        if md_tag[md_idx] >= 48 and md_tag[md_idx] <= 57:
            nmatches *= 10
            nmatches += md_tag[md_idx] - 48
            md_idx += 1
            continue
        else:
            # save matches up to this point, skipping insertions
            for x from 0 <= x < nmatches:
                while s[s_idx] == '+':
                    s_idx += 1
                s_idx += 1
            while s[s_idx] == '+':
                s_idx += 1

            r_idx += nmatches
            nmatches = 0
            if md_tag[md_idx] == '^':
                md_idx += 1
                while md_tag[md_idx] >= 65 and md_tag[md_idx] <= 90:
                    assert s[s_idx] == '-'
                    s_idx += 1
                    md_idx += 1
            else:
                # save mismatch and change to lower case
                s[s_idx] = '*'
                s_idx += 1
                r_idx += 1
                md_idx += 1

    # save matches up to this point, skipping insertions
    for x from 0 <= x < nmatches:
        while s[s_idx] == '+':
            s_idx += 1
        s_idx += 1
    while s[s_idx] == '+':
        s_idx += 1

    seq = PyBytes_FromStringAndSize(s, s_idx)
    free(s)

    return seq


def bam2stats_window_count(AlignmentFile samfile,
                           region=None,
                           chromosomes=None,
                           window_size=1000,
                           fasta=None):
    '''
    '''
    cdef AlignedSegment read
    cdef uint32_t idx = 0
    cdef uint32_t length = 0
    cdef int nflags = len(FLAGS)
    cdef int tid = -1
    cdef int flag
    cdef uint32_t f = 0
    cdef uint32_t offset = 0
    cdef uint32_t ref_window = 0
    cdef uint32_t x = 0
    cdef uint32_t _window_size = window_size
    cdef int ninput = 0
    cdef int nunmapped = 0
    cdef char code
    cdef int32_t ref_pos
    cdef int32_t next_window_start
    cdef int32_t window_start
    cdef int32_t window_end
    cdef char * sequence
    cdef int nmatches, nmismatches, ndeletions, ninsertios

    cdef int idx_alignment_starts = nflags
    cdef int idx_alignment_ends = idx_alignment_starts + 1
    cdef int idx_alignments = idx_alignment_ends + 1
    cdef int idx_matches = idx_alignments + 1
    cdef int idx_mismatches = idx_matches + 1
    cdef int idx_insertions = idx_mismatches + 1
    cdef int idx_deletions = idx_insertions + 1
    cdef int idx_hard_clipped = idx_deletions + 1
    cdef int idx_hard_clipped_bases = idx_hard_clipped + 1
    cdef int idx_soft_clipped = idx_hard_clipped_bases + 1
    cdef int idx_soft_clipped_bases = idx_soft_clipped + 1
    cdef int idx_gc = idx_soft_clipped_bases + 1
    cdef int idx_at = idx_gc + 1

    # build list of column names
    # add flags to cloumns
    columns = []
    f = 1
    for x from 0 <= x < nflags:
        columns.append(FLAGS[f])
        f = f << 1

    columns.extend([
        "alignment_starts",
        "alignment_ends",
        "alignments",
        "matched_bases",
        "mismatched_bases",
        "inserted_bases",
        "deleted_bases",
        "hard_clipped_alignments",
        "hard_clipped_bases",
        "soft_clipped_alignments",
        "soft_clipped_bases",
        "bases_gc",
        "bases_at"])

    assert idx_at == len(columns) - 1

    # build counting matrix
    contigs = samfile.references
    lengths = samfile.lengths

    if region is not None:
        itr = samfile.fetch(region)
        contig, start, end = parse_region_string(region)

        pairs = []
        for c, l in zip(contigs, lengths):
            if c == contig:
                pairs.append((c, l))
            else:
                pairs.append((c, 0))
        contigs, lengths = zip(*pairs)
    else:
        itr = samfile.fetch()

    cdef uint32_t total_length = sum(lengths) + _window_size * len(lengths)
    cdef uint32_t ncontigs = len(contigs)
    cdef int ncolumns = len(columns)
    cdef numpy.ndarray offsets = numpy.zeros(
        ncontigs + 1, dtype=numpy.int64)

    for idx, length in enumerate(lengths):
        offsets[idx] = offset
        if length > 0:
            offset += (length + _window_size) // _window_size

    cdef uint32_t total_windows = offset
    offsets[idx + 1] = offset

    cdef numpy.ndarray window_counts = numpy.zeros(
        (total_windows, ncolumns),
        dtype=numpy.int64)

    cdef int32_t [:, :] window_counts_view = window_counts

    E.debug("total_windows={}, total_length={}, nflags={}, ncolumns={}, ncontigs={}".format(
            total_windows, total_length, nflags, ncolumns, ncontigs))

    E.info("computing read counts in windows")

    for read in itr:

        flag = read._delegate.core.flag
        ninput += 1

        if ninput % 1000 == 0:
            E.debug("processing read number {} at {}:{}".format(
                    ninput, read.reference_name, read.reference_start))

        tid = read.tid
        if tid < 0:
            continue

        if flag & 4:
            nunmapped += 1
            continue

        window_start = offsets[tid] + read.reference_start // _window_size
        window_end = offsets[tid] + read.reference_end // _window_size

        window_counts[window_start][idx_alignment_starts] += 1
        window_counts[window_end][idx_alignment_ends] += 1

        # count hard/soft-clipping
        cigar = read.cigartuples
        if cigar[0][0] == 5:
            window_counts[window_start][idx_hard_clipped] += 1
            window_counts[window_start][idx_hard_clipped_bases] += cigar[0][1]
        elif cigar[0][0] == 4:
            window_counts[window_start][idx_soft_clipped] += 1
            window_counts[window_start][idx_soft_clipped_bases] += cigar[0][1]

        if cigar[-1][0] == 5:
            window_counts[window_end][idx_hard_clipped] += 1
            window_counts[window_end][idx_hard_clipped_bases] += cigar[-1][1]
        elif cigar[-1][0] == 4:
            window_counts[window_end][idx_soft_clipped] += 1
            window_counts[window_end][idx_soft_clipped_bases] += cigar[-1][1]

        for idx from window_start <= idx <= window_end:
            window_counts[idx][idx_alignments] += 1
            f = 1
            for x from 0 <= x < nflags:
                if flag & f:
                    window_counts[idx][x] += 1
                f = f << 1

        idx = window_start
        md_seq = count_md_tag_mismatches(read)
        ref_pos = read.reference_start
        nmatches = nmismatches = ndeletions = ninsertions = 0
        # border of next window
        next_window_start = offsets[tid] + ((ref_pos // _window_size + 1) * _window_size)
        for code in md_seq:

            if code == '=':
                nmatches += 1
                ref_pos += 1
            elif code == "+":
                ninsertions += 1
            elif code == "-":
                ndeletions += 1
                ref_pos += 1
            elif code == "*":
                nmismatches += 1
                ref_pos += 1

            if ref_pos >= next_window_start:
                window_counts[idx][idx_matches] += nmatches
                window_counts[idx][idx_mismatches] += nmismatches
                window_counts[idx][idx_deletions] += ndeletions
                window_counts[idx][idx_insertions] += ninsertions

                next_window_start += _window_size
                idx += 1

        window_counts[idx][idx_matches] += nmatches
        window_counts[idx][idx_mismatches] += nmismatches
        window_counts[idx][idx_deletions] += ndeletions
        window_counts[idx][idx_insertions] += ninsertions

    # compute G+C content
    if fasta:
        E.info("computing G+C content windows")
        idx = 0
        for contig, length in zip(contigs, lengths):
            if length == 0:
                idx += 1
                continue

            offset = offsets[idx]
            s = fasta.fetch(contig)
            sequence = s
            for x from 0 <= x < length:
                # do not increment at x == 0
                if x and x % _window_size == 0:
                    offset += 1

                code = sequence[x]
                if code > 97:
                    code -= 32

                if code == 'C' or code == 'G':
                    window_counts[offset][idx_gc] += 1
                elif code == 'A' or code == 'T':
                    window_counts[offset][idx_at] += 1

            idx += 1

    # set window coordinates
    E.info("setting window coordinates")
    chromosomes, starts, ends = [], [], []
    idx = 0
    for contig, length in zip(contigs, lengths):
        if length == 0:
            idx += 1
            continue

        offset = offsets[idx]
        next_offset = offsets[idx + 1]
        nwindows = next_offset - offset
        chromosomes.extend([contig] * nwindows)
        starts.extend(numpy.arange(0,
                                   window_size * nwindows,
                                   window_size))
        ends.extend(numpy.arange(window_size,
                                 window_size * (nwindows + 1),
                                 window_size))
        idx += 1

    E.info("building dataframe")
    window_df = pandas.DataFrame(window_counts,
                                 columns=columns,
                                 index=[chromosomes, starts, ends])
    window_df.index.names = ["contig", "start", "end"]

    return window_df


def is_paired(bamfile, alignments=1000):
    '''check if a `bamfile` contains paired end data

    The method reads at most the first *alignments* and returns
    True if any of the alignments are paired.
    '''

    samfile = pysam.AlignmentFile(bamfile)
    n = 0
    for read in samfile:
        if read.is_paired:
            break
        n += 1
        if n == alignments:
            break

    samfile.close()

    return n != alignments


def estimateInsertSizeDistribution(bamfile,
                                   alignments=10000,
                                   n=10,
                                   method="picard",
                                   similarity_threshold=1.0,
                                   max_chunks=1000):
    '''estimate insert size from a subset of alignments in a bam file.

    Several methods are implemented.

    picard
        The method works analogous to picard by restricting the estimates
        to a core distribution. The core distribution is defined as all
        values that lie within n-times the median absolute deviation of
        the full data set.
    convergence
        The method works similar to ``picard``, but continues reading
        `alignments` until the mean and standard deviation stabilize.
        The values returned are the median mean and median standard
        deviation encountered.

    The method `convergence` is suited to RNA-seq data, as insert sizes
    fluctuate siginificantly depending on the current region
    being looked at.

    Only mapped and proper pairs are considered in the computation.

    Returns
    -------
    mean : float
       Mean of insert sizes.
    stddev : float
       Standard deviation of insert sizes.
    npairs : int
       Number of read pairs used for the estimation
    method : string
       Estimation method
    similarity_threshold : float
       Similarity threshold to apply.
    max_chunks : int
       Maximum number of chunks of size `alignments` to be used
       in the convergence method.

    '''

    assert is_paired(bamfile), \
        'can only estimate insert size from' \
        'paired bam files'

    samfile = pysam.AlignmentFile(bamfile)

    def get_core_distribution(inserts, n):
        # compute median absolute deviation
        raw_median = numpy.median(inserts)
        raw_median_dev = numpy.median(numpy.absolute(inserts - raw_median))

        # set thresholds
        threshold_min = max(0, raw_median - n * raw_median_dev)
        threshold_max = raw_median + n * raw_median_dev

        # define core distribution
        return inserts[numpy.logical_and(inserts >= threshold_min,
                                         inserts <= threshold_max)]

    if method == "picard":

        # only get first read in pair to avoid double counting
        inserts = numpy.array(
            [read.template_length for read in samfile.head(n=alignments)
             if read.is_proper_pair
             and not read.is_unmapped
             and not read.mate_is_unmapped
             and not read.is_read1
             and not read.is_duplicate
             and read.template_length > 0])
        core = get_core_distribution(inserts, n)

        return numpy.mean(core), numpy.std(core), len(inserts)

    elif method == "convergence":

        means, stds, counts = [], [], []
        last_mean = 0
        iteration = 0
        while iteration < max_chunks:

            inserts = numpy.array(
                [read.template_length for read in samfile.head(
                    n=alignments,
                    multiple_iterators=False)
                 if read.is_proper_pair
                 and not read.is_unmapped
                 and not read.mate_is_unmapped
                 and not read.is_read1
                 and not read.is_duplicate
                 and read.template_length > 0])
            core = get_core_distribution(inserts, n)
            means.append(numpy.mean(core))
            stds.append(numpy.std(core))
            counts.append(len(inserts))
            mean_core = get_core_distribution(numpy.array(means), 2)
            mm = numpy.mean(mean_core)
            if abs(mm - last_mean) < similarity_threshold:
                break
            last_mean = mm

        return numpy.median(means), numpy.median(stds), sum(counts)
    else:
        raise ValueError("unknown method '%s'" % method)


def estimateTagSize(bamfile,
                    alignments=10,
                    multiple="error"):
    '''estimate tag/read size from first alignments in file.

    Arguments
    ---------
    bamfile : string
       Filename of :term:`bam` formatted file
    alignments : int
       Number of alignments to inspect
    multiple : string
       How to deal if there are multiple tag sizes present.
       ``error`` will raise a warning, ``mean`` will return the
       mean of the read lengths found. ``uniq`` will return a
       unique list of read sizes found. ``all`` will return all
       read sizes encountered.

    Returns
    -------
    size : int
       The read size (actual, mean or list of read sizes)

    Raises
    ------
    ValueError
       If there are multiple tag sizes present and `multiple` is set to
       `error`.

    '''
    samfile = pysam.AlignmentFile(bamfile)
    sizes = [read.rlen for read in samfile.head(alignments)]
    mi, ma = min(sizes), max(sizes)

    if mi == 0 and ma == 0:
        sizes = [read.inferred_length for read in samfile.head(alignments)]
        # remove 0 sizes (unaligned reads?)
        sizes = [x for x in sizes if x > 0]
        mi, ma = min(sizes), max(sizes)

    if mi != ma:
        if multiple == "error":
            raise ValueError('multiple tag sizes in %s: %s' % (bamfile, sizes))
        elif multiple == "mean":
            mi = int(sum(sizes) / len(sizes))
        elif multiple == "uniq":
            mi = list(sorted(set(sizes)))
        elif multiple == "all":
            return sizes

    return mi


def getNumberOfAlignments(bamfile):
    '''return number of alignments in bamfile.
    '''
    samfile = pysam.AlignmentFile(bamfile)
    return samfile.mapped


def getNumReads(bamfile):
    '''count number of reads in bam file.

    This methods works through pysam.idxstats.

    Arguments
    ---------
    bamfile : string
        Filename of :term:`bam` formatted file. The file needs
        to be indexed.
    Returns
    -------
    nreads : int
        Number of reads
    '''

    lines = pysam.idxstats(bamfile).splitlines()

    try:
        nreads = sum(
            map(int, [x.split("\t")[2]
                      for x in lines if not x.startswith("#")]))

    except IndexError as msg:
        raise IndexError(
            "can't get number of reads from bamfile, msg=%s, data=%s" %
            (msg, lines))
    return nreads


def isStripped(bamfile):
    '''
    Check if the sequence is stripped in a bam file.
    '''
    bam = pysam.AlignmentFile(bamfile)
    read = bam.next()
    if read.seq is None:
        return True
    else:
        return False


def merge_pairs(AlignmentFile input_samfile,
                outfile,
                min_insert_size = 0,
                max_insert_size = 400,
                bed_format = None ):
    '''merge paired ended data.

    For speed reasons, the aligned region is only approximated using
    the read length. Indels within reads are ignored. Thus bed coordinates
    might be a few residues off.

    The strand is always set to '+'.

    Pairs with a maximum insert size larger than *max_insert_size* are removed.

    If `bed_format` is a number, only the first x columns will be output.
    '''

    cdef int ninput = 0
    cdef int nremoved_unmapped = 0
    cdef int nremoved_insert = 0
    cdef int nremoved_contig = 0
    cdef int nremoved_unpaired = 0
    cdef int nremoved_coordinate = 0
    cdef int nremoved_take_only_second = 0
    cdef int noutput = 0
    cdef int flag
    cdef int isize

    cdef AlignedSegment read
    cdef int c_max_insert_size = max_insert_size
    cdef int c_min_insert_size = min_insert_size
    cdef int start, end, xstart, xend
    cdef int take_columns = 6

    # point to array of contig lengths
    cdef uint32_t *contig_sizes = input_samfile.header.ptr.target_len

    if bed_format != None:
        if bed_format < 3 or bed_format > 6:
            raise ValueError("a bed file must have at least 3 and at most 6 columns")
        take_columns = bed_format

    for read in input_samfile:
        ninput += 1

        flag = read._delegate.core.flag
        # remove unmapped reads
        if flag & 4:
            nremoved_unmapped += 1
            continue

        if read.pos < read.mpos:
            # lower coordinate than mate, ignore
            nremoved_coordinate += 1
            continue
        elif read.pos == read.mpos and flag & 64:
            # disambiguate, ignore first in pair
            nremoved_take_only_second += 1
            continue
        else:
            # taking the downstream pair allows to incl
            xstart = read.next_reference_start
            xend = read.reference_end
            if xstart < xend:
                start, end = xstart, xend
            else:
                start, end = xend, xstart

        # remove unpaired
        if not flag & 2:
            nremoved_unpaired += 1
            continue

        if read.tid != read.mrnm:
            nremoved_contig += 1
            continue

        # isize can be negative - depending on the pair orientation
        isize = abs(read.isize)
        if (c_max_insert_size and isize > c_max_insert_size) or \
           (c_min_insert_size and isize < c_min_insert_size):
            nremoved_insert += 1
            continue

        # truncate at contig end - overhanging reads might cause problems with chrM
        if end > contig_sizes[read.mrnm]:
            end = contig_sizes[read.mrnm]

        # count output pair as two so that it squares with ninput
        noutput += 2

        if take_columns == 3:
            outfile.write("%s\t%i\t%i\n" %
                          (input_samfile.getrname(read.tid),
                           start, end))
        elif take_columns == 4:
            outfile.write("%s\t%i\t%i\t%s\n" %
                          (input_samfile.getrname(read.tid),
                           start, end,
                           read.qname))
        elif take_columns == 5:
            outfile.write("%s\t%i\t%i\t%s\t%i\n" %
                          (input_samfile.getrname(read.tid),
                           start, end,
                           read.qname,
                           read.mapq,
                       ))
        else:
            # As we output the downstream read, reverse orientation
            if read.is_reverse:
                strand = '-'
            else:
                strand = '+'
            outfile.write("%s\t%i\t%i\t%s\t%i\t%s\n" %
                          (input_samfile.getrname(read.tid),
                           start, end,
                           read.qname,
                           read.mapq,
                           strand
                          ))

    c = E.Counter()
    c.input = ninput
    c.removed_insert = nremoved_insert
    c.removed_contig = nremoved_contig
    c.removed_unmapped = nremoved_unmapped
    c.removed_unpaired = nremoved_unpaired
    c.removed_take_only_second = nremoved_take_only_second
    c.removed_coordinate = nremoved_coordinate
    c.output = noutput

    return c


def bams2bam_filter(AlignmentFile genome_samfile,
                    AlignmentFile output_samfile,
                    AlignmentFile output_mismapped,
                    AlignmentFile transcripts_samfile,
                    AlignmentFile junctions_samfile,
                    transcripts,
                    regions = None,
                    unique = False,
                    remove_contigs = None,
                    colour_mismatches = False,
                    ignore_mismatches = False,
                    ignore_junctions = True,
                    ignore_transcripts = False ):
    '''
    To conserve memory, the tid and NM flag from *transcripts_samfile*
    are packed into memory. As a consequence, this method requires
    max(NM) < 256 (2^8) and max(tid) < 16777216 (2^24)

    If *unique* is set, only uniquely matching reads will be output.

    If *remove_contigs* is given, contigs that are in remove_contigs will
    be removed. Note that this will only remove the alignment, but not
    all alignments of a particuluar read and the NH flag will *NOT* be
    updated.

    If *colour_mismatches* is set, the ``CM`` tag will be used
    to count differences. By default, the ``NM`` tag is used.
    The tag that is used needs to present in both *transcripts_samfile*
    and *genome_samfile*.

    If *ignore_mismatches* is set, the number of mismatches is ignored.

    If *regions* is given, alignments overlapping regions will be removed.

    '''
    cdef int ninput = 0
    cdef int nunmapped_genome = 0
    cdef int nunmapped_transcript = 0
    cdef int nmismapped = 0
    cdef int noutput = 0
    cdef int nnotfound = 0
    cdef int nlocation_ok = 0
    cdef int nnonunique = 0
    cdef int ntested = 0
    cdef int nremoved_contigs = 0
    cdef int nremoved_junctions = 0
    cdef int nskipped_junctions = 0
    cdef int nadded_junctions = 0
    cdef int nremoved_regions = 0

    cdef bint c_unique = unique
    cdef bint c_test_mismatches = not ignore_mismatches
    cdef bint c_test_junctions = not ignore_junctions
    cdef bint c_test_transcripts = not ignore_transcripts
    cdef bint c_test_regions = regions
    cdef int * remove_contig_tids
    cdef int nremove_contig_tids = 0
    cdef AlignedSegment read
    cdef AlignedSegment match

    # build index
    # this method will start indexing from the current file position
    # if you decide
    cdef int ret = 1
    cdef int x
    cdef bam1_t * b = <bam1_t*>calloc(1, sizeof( bam1_t))
    cdef uint64_t pos
    cdef uint8_t * v
    cdef int32_t nm
    cdef int32_t nh
    cdef int tid
    cdef int transript_tid
    cdef long val
    cdef bint skip = 0
    cdef int * map_tid2tid

    cdef int32_t read_mismatches = 0

    # decide which tag to use
    cdef char * nm_tag = 'NM'
    cdef char * cm_tag = 'CM'
    cdef char * tag

    if colour_mismatches:
        tag = cm_tag
    else:
        tag = nm_tag

    # set with junctions that are ignored
    skip_junctions = set()

    # build list of junctions
    if c_test_junctions:
        E.info("building junction read index")

        def _gen2(): return array.array('B')
        junctions_index = collections.defaultdict(_gen2)
        ret = 1
        while ret > 0:
            ret = bam_read1(hts_get_bgzfp(junctions_samfile.htsfile),
                            b)
            if ret > 0:
                # ignore unmapped reads
                if b.core.flag & 4: continue
                qname = pysam_bam_get_qname(b)
                v = bam_aux_get(b, tag)
                nm = <int32_t>bam_aux2i(v)
                junctions_index[qname].append(nm)

        E.info( "built index for %i junction reads" % len(junctions_index))

        map_tid2tid = <int*>calloc(len(junctions_samfile.references), sizeof(int))

        for x, contig_name in enumerate(junctions_samfile.references):
            map_tid2tid[x] = genome_samfile.gettid(contig_name)

    if c_test_transcripts:
        E.info( "building transcriptome read index" )
        # L = 4 bytes
        def _gen(): return array.array('L')
        transcriptome_index = collections.defaultdict(_gen)
        ret = 1
        while ret > 0:
            ret = bam_read1(hts_get_bgzfp(transcripts_samfile.htsfile),
                            b)
            if ret > 0:
                # ignore unmapped reads
                if b.core.flag & 4: continue
                qname = pysam_bam_get_qname(b)
                tid = b.core.tid
                v = bam_aux_get(b, tag)
                nm = <int32_t>bam_aux2i(v)
                transcriptome_index[qname].append( (tid << 8) | nm )

        E.info( "built index for %i transcriptome reads" % len(transcriptome_index))

    bam_destroy1( b )

    # setup list of contigs to remove:
    if remove_contigs:
        nremove_contig_tids = len(remove_contigs)
        remove_contig_tids = <int*>malloc( sizeof(int) * nremove_contig_tids )
        for x, rname in enumerate( remove_contigs):
            remove_contig_tids[x] = genome_samfile.gettid( rname )

    E.info( "starting filtering" )

    for read in genome_samfile:

        ninput += 1
        # if ninput > 10000: break

        # is unmapped?
        if read._delegate.core.flag & 4:
            nunmapped_genome += 1
            noutput += 1
            output_samfile.write( read )
            continue

        # optionally remove non-unique reads
        if c_unique:
            v = bam_aux_get(read._delegate, 'NH')
            if v != NULL:
                nh = <int32_t>bam_aux2i(v)
                if nh > 1:
                    nnonunique += 1
                    continue

        # optionally remove reads matched to certain contigs
        if nremove_contig_tids:
            skip = 0
            for x from 0 <= x < nremove_contig_tids:
                if remove_contig_tids[x] == read.tid:
                    nremoved_contigs += 1
                    skip = 1
                    break
            if skip: continue

        g_contig = genome_samfile.getrname( read.tid )

        # optionally remove reads mapped to certain regions
        if c_test_regions:
            intervals = regions.get( g_contig, read.pos, read.aend )
            skip = 0
            for start, end, value in intervals:
                if read.overlap( start, end ):
                    nremoved_regions += 1
                    skip = 1
                    break

            if skip: continue

        if c_test_junctions:
            if read.qname in junctions_index:
                # can compress index before, depends on
                # how many reads you expect to test
                nm = min(junctions_index[read.qname])
                read_mismatches = read.opt(tag)
                if nm > read_mismatches:
                    nremoved_junctions += 1
                    continue
                else:
                    skip_junctions.add(read.qname)

        # set mapped = True, if read is mapped to transcripts
        # set location_ok = True, if read matches in expected location
        # according to transcripts
        location_ok = False
        mapped = False

        if c_test_transcripts:
            # get transcripts that read matches to
            try:
                matches = transcriptome_index[read.qname]
            except KeyError:
                nnotfound += 1
                matches = None

            if matches:

                if c_test_mismatches:
                    read_mismatches = read.opt(tag)

                for val in matches:
                    transcript_tid = val >> 8
                    # ignore reads that are mapped to transcripts with
                    # more mismatches than the genomic location
                    if c_test_mismatches:
                        nm = val & 255
                        if nm > read_mismatches: continue

                    mapped = True

                    # find transcript
                    transcript_id = transcripts_samfile._getrname( transcript_tid )
                    gtfs = transcripts[transcript_id]
                    t_contig, t_start, t_end = gtfs[0].contig, gtfs[0].start, gtfs[-1].end

                    # does read map to genomic transcript location?
                    if g_contig == t_contig and t_start <= read.pos <= t_end:
                        location_ok = True
                        break

        if location_ok:
            ntested += 1
            nlocation_ok += 1
            noutput += 1
            output_samfile.write( read )

        elif mapped:
            ntested += 1
            nmismapped += 1
            if output_mismapped:
                output_mismapped.write( read )

        else:
            nunmapped_transcript += 1
            noutput += 1
            output_samfile.write( read )

    # add junctions
    if c_test_junctions:
        E.info("adding junctions")
        junctions_samfile.reset()

        # rebuild contig map for junctions:
        if remove_contigs:
            remove_contig_tids = <int*>malloc( sizeof(int) * nremove_contig_tids )
            for x, rname in enumerate( remove_contigs):
                remove_contig_tids[x] = junctions_samfile.gettid( rname )

        for read in junctions_samfile:

            # optionally remove reads matched to certain contigs
            if nremove_contig_tids:
                skip = 0
                for x from 0 <= x < nremove_contig_tids:
                    if remove_contig_tids[x] == read.tid:
                        nremoved_contigs += 1
                        skip = 1
                    break
                if skip:
                    nskipped_junctions += 1
                    continue

            if read.qname in skip_junctions:
                nskipped_junctions += 1
            else:
                nadded_junctions += 1
                noutput += 1
                # map tid from junction database to genome database
                read.tid = map_tid2tid[read.tid]
                output_samfile.write(read)

        free( map_tid2tid )

    c = E.Counter()
    c.input = ninput
    c.removed_nonunique = nnonunique
    c.removed_mismapped = nmismapped
    c.removed_contigs = nremoved_contigs
    c.removed_junctions = nremoved_junctions
    c.removed_regions = nremoved_regions
    c.skipped_junction_reads = len(skip_junctions)
    c.skipped_junctions = nskipped_junctions
    c.added_junctions = nadded_junctions
    c.output = noutput
    c.unmapped_genome = nunmapped_genome
    c.unmapped_transcript = nunmapped_transcript
    c.notfound = nnotfound
    c.location_ok = nlocation_ok
    c.location_tested = ntested

    if nremove_contig_tids:
        free( remove_contig_tids )

    E.info( "filtering finished" )

    return c
