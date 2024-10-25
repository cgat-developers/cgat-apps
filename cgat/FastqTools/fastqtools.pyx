"""Utility functions for the bam2stats utility."""

# Specific imports for optimisation
import collections
import pysam
from pysam.libchtslib cimport cts_bamfile_t, cts_alignment_t, bam_aux_get
from pysam.libcbcf cimport VariantFile, VariantRecord, VariantRecordSample
from pysam.libcfaidx cimport FastxFile, FastqProxy
from pysam.libctabix cimport TabixFile, ctbx_index_t, ctbx_iter_t
from libc.stdint cimport uint32_t
import cgatcore.experiment as E
from cgat.Genomics import reverse_complement
import numpy
cimport numpy

# Define data type for numpy arrays
DTYPE = numpy.uint32
ctypedef numpy.uint32_t DTYPE_t

cdef count_diagonals(
    sequence,
    forward_hash,
    reverse_hash,
    uint32_t query_sequence_length,
    uint32_t kmer_size
):
    """Counts diagonal matches for k-mers in a sequence against forward and reverse hashes."""
    cdef uint32_t x, count
    cdef uint32_t offset = len(sequence)
    cdef uint32_t l = len(sequence) + query_sequence_length
    cdef numpy.ndarray[DTYPE_t, ndim=1] forward_counts = numpy.zeros(l, dtype=DTYPE)
    cdef numpy.ndarray[DTYPE_t, ndim=1] reverse_counts = numpy.zeros(l, dtype=DTYPE)
    
    for x from 0 <= x < len(sequence) - kmer_size:
        kmer = sequence[x: x + kmer_size]
        if kmer in forward_hash:
            for count in forward_hash[kmer]:
                forward_counts[count - x + offset] += 1

        if kmer in reverse_hash:
            for count in reverse_hash[kmer]:
                reverse_counts[count - x + offset] += 1

    return forward_counts.max(), reverse_counts.max()


def filter_by_sequence(
    query_sequence,
    FastxFile in_stream1,
    FastxFile in_stream2,
    outf_matched1,
    outf_matched2,
    outf_unmatched1,
    outf_unmatched2,
    uint32_t kmer_size=10,
    uint32_t min_kmer_matches=20
):
    """Filters paired-end reads based on k-mer matches with a query sequence."""
    
    # Prepare k-mer hashes
    _query_sequence = query_sequence.encode("ascii")
    reverse_sequence = reverse_complement(query_sequence)
    cdef uint32_t query_sequence_length = len(query_sequence)
    forward_hash = collections.defaultdict(list)
    reverse_hash = collections.defaultdict(list)
    
    for x from 0 <= x < query_sequence_length - kmer_size:
        forward_hash[query_sequence[x:x + kmer_size]].append(x)
        reverse_hash[reverse_sequence[x: x + kmer_size]].append(x)

    # Track counts for matched and unmatched sequences
    cdef FastqProxy read1, read2
    cdef uint32_t ninput = 0
    cdef uint32_t nmatched = 0
    cdef uint32_t nunmatched = 0
    cdef uint32_t fc1, rc1, fc2, rc2

    for read1, read2 in zip(in_stream1, in_stream2):
        # Calculate diagonal match counts
        fc1, rc1 = count_diagonals(
            read1.sequence,
            forward_hash,
            reverse_hash,
            query_sequence_length,
            kmer_size
        )
        fc2, rc2 = count_diagonals(
            read2.sequence,
            forward_hash,
            reverse_hash,
            query_sequence_length,
            kmer_size
        )

        # Write to matched or unmatched output based on k-mer match threshold
        if max(fc1, rc1, fc2, rc2) > min_kmer_matches:
            nmatched += 1
            outf_matched1.write(str(read1) + "\n")
            outf_matched2.write(str(read2) + "\n")
        else:
            nunmatched += 1
            outf_unmatched1.write(str(read1) + "\n")
            outf_unmatched2.write(str(read2) + "\n")
        
        ninput += 1
        if ninput % 1000 == 0:
            E.info(
                f"iteration: {ninput}, matched={nmatched}, unmatched={nunmatched}, "
                f"permille_matched={1000.0 * nmatched / ninput}"
            )

    # Return a summary count of inputs, matches, and unmatched reads
    c = E.Counter()
    c.input = ninput
    c.matched = nmatched
    c.unmatched = nunmatched
    return c
