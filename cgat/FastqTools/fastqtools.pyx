"""Utility functions for the bam2stats utility."""

import collections

import pysam
from pysam.libchtslib cimport *
from pysam.libcbcf cimport VariantFile, VariantRecord, VariantRecordSample
from pysam.libcfaidx cimport FastxFile, FastqProxy
from pysam.libctabix cimport TabixFile
from libc.string cimport strchr
from libc.stdint cimport int8_t
from libc.stdio cimport puts, printf
from cpython cimport array as c_array

import cgatcore.experiment as E
from cgat.Genomics import reverse_complement
import numpy
cimport numpy
DTYPE = numpy.uint32
ctypedef numpy.uint32_t DTYPE_t

cdef count_diagonals(sequence,
                     forward_hash,
                     reverse_hash,
                     uint32_t query_sequence_length,
                     uint32_t kmer_size):
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
        uint32_t min_kmer_matches=20):

    _query_sequence = query_sequence.encode("ascii")
    reverse_sequence = reverse_complement(query_sequence)
    cdef uint32_t query_sequence_length = len(query_sequence)
    cdef uint32_t x
    cdef char * s = _query_sequence

    forward_hash = collections.defaultdict(list)
    reverse_hash = collections.defaultdict(list)
    for x from 0 <= x < query_sequence_length - kmer_size:
        forward_hash[query_sequence[x:x + kmer_size]].append(x)
        reverse_hash[reverse_sequence[x: x + kmer_size]].append(x)

    cdef FastqProxy read1, read2
    cdef uint32_t ninput = 0
    cdef uint32_t nmatched = 0
    cdef uint32_t nunmatched = 0
    cdef uint32_t fc1, rc1, fc2, rc2

    for read1, read2 in zip(in_stream1, in_stream2):
    
        fc1, rc1 = count_diagonals(
            read1.sequence,
            forward_hash,
            reverse_hash,
            query_sequence_length,
            kmer_size)
        fc2, rc2 = count_diagonals(
            read2.sequence,
            forward_hash,
            reverse_hash,
            query_sequence_length,
            kmer_size)

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
            E.info("iteration: {}, matched={}, unmatched={}, permille_matched={}".format(
                ninput, nmatched, nunmatched, 1000.0 * nmatched / ninput))

    c = E.Counter()
    c.input = ninput
    c.matched = nmatched
    c.unmatched = nunmatched
    return c
