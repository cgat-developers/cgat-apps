from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
import collections
import cgatcore.experiment as E
import numpy
import cgat.Stats as Stats

# Define named tuples for results
PeakShapeResult = collections.namedtuple(
    "PeakShapeResult",
    "interval_width npeaks peak_center peak_width peak_height peak_relative_pos "
    "nreads median closest_half_height furthest_halfheight bins counts"
)
PeakShapeCounts = collections.namedtuple(
    "PeakShapeCounts",
    "nreads median counts"
)

# Define classes and methods
cdef class Counter:
    '''Base class for counters computing densities from genomic data.'''

    cdef smooth_method
    
    def __init__(self, smooth_method=None):
        self.smooth_method = smooth_method

    def countAroundPos(self, 
                       infile,
                       contig, 
                       int pos,
                       bins,
                       **kwargs):
        '''Count and bin in bins, returning a PeakShapeCounts tuple.'''
        cdef int xstart, xend, rstart, rend, i, offset
        nreads, counts = self.coverageInInterval(infile, contig, max(0, pos + bins[0]), pos + bins[-1], **kwargs)
        nbins = len(bins) - 1
        hist = numpy.zeros(nbins, dtype=numpy.int64)
        cdef int lcounts = len(counts)

        if self.smooth_method:
            smoothed_counts = self.smooth_method(counts)

        offset = -bins[0]
        xstart = offset + bins[0]
        for i from 0 <= i < nbins:
            xend = offset + bins[i+1]
            if xstart >= 0 and xend < lcounts:
                hist[i] = sum(counts[xstart:xend]) 
            xstart = xend

        result = PeakShapeCounts._make((nreads, numpy.median(counts), hist))
        return result

    # Other methods remain similar, adjusted as necessary

cdef class CounterBam(Counter):
    '''Compute densities in intervals from BAM files.'''

    cdef int shift

    def __init__(self, shift=0, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.shift = shift

    def coverageInInterval(self, AlignmentFile samfile, contig, int start, int end):
        '''Return coverage in window on *contig* bounded by *start* and *end*.'''
        cdef int shift = self.shift
        cdef int nreads = 0
        cdef int interval_width = end - start
        cdef int *ccounts = <int*>calloc(interval_width, sizeof(int))
        
        if interval_width <= 0 or contig not in samfile.references:
            return 0, numpy.zeros(0)

        if shift:
            xstart, xend = max(0, start - shift // 2), max(0, end - shift // 2)
            for read in samfile.fetch(contig, xstart, xend):
                if read.is_reverse: continue
                nreads += 1
                rstart = read.pos
                rend = rstart + 2 * shift // 2
                rstart, rend = max(0, rstart - xstart), min(interval_width, rend - xstart)
                for i from rstart <= i < rend: ccounts[i] += 1

            xstart, xend = start + shift // 2, end + shift // 2
            for read in samfile.fetch(contig, xstart, xend):
                if not read.is_reverse: continue
                nreads += 1
                rend = read.aend
                rstart = max(0, rend - 2 * shift // 2) if rend else continue
                rstart, rend = max(0, rstart - xstart), min(interval_width, rend - xstart)
                for i from rstart <= i < rend: ccounts[i] += 1
        else:
            for read in samfile.fetch(contig, start, end):
                nreads += 1
                rstart = max(0, read.pos - start)
                rend = min(interval_width, read.aend - start) if read.aend else continue
                for i from rstart <= i < rend: ccounts[i] += 1

        counts = numpy.zeros(interval_width, dtype=numpy.int64)
        for i from 0 <= i < interval_width: counts[i] = ccounts[i]
        free(ccounts)
        return nreads, counts

cdef class CounterBigwig(Counter):
    '''Compute densities in intervals from BigWig files.'''

    def __cinit__(self, shift=0, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def coverageInInterval(self, wigfile, contig, int start, int end):
        '''Return coverage in window on *contig* bounded by *start* and *end*.'''
        cdef int interval_width = end - start
        if interval_width <= 0:
            return 0, numpy.zeros(0)
        
        d = wigfile.summarize(contig, start, end, interval_width)
        nreads = sum(d.valid_count)
        return nreads, d.sum_data
