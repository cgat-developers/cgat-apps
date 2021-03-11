"""
fasta2bed.py - segment sequences
================================

:Tags: Genomics Sequences Intervals FASTA BED Conversion

Purpose
-------

This script takes a genomic sequence in :term:`fasta` format
and applies various segmentation algorithms.

The methods implemented (``--methods``) are:

cpg
   output all locations of cpg in the genome

fixed-width-windows-gc
   output fixed width windows of a certain size adding their
   G+C content as score

gaps
   ouput all locations of assembly gaps (blocks of `N`)
   in the genomic sequences

ungapped
   output ungapped locations in the genomic sequences

Usage
-----

Type::

   python fasta2bed.py --method=gap < in.fasta > out.bed


Type::

   python fasta2bed.py --help

for command line help.

Command line options
--------------------

"""

import os
import sys
import re
import tempfile
import subprocess
import glob
import collections
import pybedtools
import cgatcore.experiment as E
import cgatcore.iotools as iotools
import cgat.FastaIterator as FastaIterator
import cgatcore.iotools as iotools


def segmentWithCpG(infile, with_contig_sizes=False):
    '''segment a fasta file, output locations of CpG.'''

    ninput, nskipped, noutput = 0, 0, 0

    iterator = FastaIterator.FastaIterator(infile)

    segments, contig_sizes = [], collections.OrderedDict()

    for cur_record in iterator:
        ninput += 1
        contig = re.sub("\s.*", "", cur_record.title)
        last = None
        contig_sizes[contig] = (0, len(cur_record.sequence))
        for pos, this in enumerate(cur_record.sequence.upper()):
            if last == "C" and this == "G":
                segments.append((contig, pos - 1, pos + 1, 1.0))
            last = this

    E.info("ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped))

    if with_contig_sizes:
        return segments, contig_sizes

    return segments


def segmentWindowsCpG(infile, window_size=100, min_cpg=1):
    '''segment a fasta file based on the locations of CpG.

    Locate all CpG in sequences and centre windows of size *window_size*
    around them. Merge all windows and keep all with *min_cpg* CpG.
    '''

    cpgs, contig_sizes = segmentWithCpG(infile, with_contig_sizes=True)

    # save cpgs to temporary file
    tempf = tempfile.NamedTemporaryFile(mode="w", delete=False)
    tempf.write("\n".join(["%s\t%i\t%i\n" % (contig, start, end)
                           for contig, start, end, gc in cpgs]) + "\n")
    tempf.close()

    cpgs = pybedtools.BedTool(tempf.name)
    cpgs.set_chromsizes(contig_sizes)
    extended = cpgs.slop(b=window_size // 2)
    merged = extended.merge(o="count", c=3)
    filtered = merged.filter(lambda x: int(x.name) >= min_cpg)

    os.unlink(tempf.name)

    # return CpG content (not C+C content)
    return [(x.chrom, x.start, x.stop, float(x.name) / (x.stop - x.start) / 2)
            for x in filtered]


def segmentFixedWidthWindows(infile, window_size, window_shift):
    """return a list of fixed contig sizes."""

    ninput, nskipped, noutput = 0, 0, 0

    iterator = FastaIterator.FastaIterator(infile)
    window_shift = window_size
    # at most 50% can be gap
    gap_cutoff = int(window_size // 2)
    segments = []

    while 1:
        ninput += 1
        try:
            cur_record = next(iterator)
        except StopIteration:
            break

        if cur_record is None:
            break
        contig = re.sub("\s.*", "", cur_record.title)
        seq = cur_record.sequence
        size = len(cur_record.sequence)

        for x in range(0, size, window_shift):
            s = seq[x:x + window_size].upper()
            gc, at = 0, 0
            for c in s:
                if c in "GC":
                    gc += 1
                elif c in "AT":
                    at += 1

            # skip segments containing mostly gaps
            if window_size - (gc + at) > gap_cutoff:
                nskipped += 1
                continue

            segments.append(
                (contig, x, x + window_size, float(gc) / (gc + at)))
        noutput += 1

    E.info("ninput=%i, noutput=%i, nskipped_windows=%i" %
           (ninput, noutput, nskipped))

    return segments


def gapped_regions(seq, gap_chars):
    '''iterator yielding gapped regions in seq.'''
    # check to see if there is a gap at the beginning
    is_gap = seq[0] in gap_chars
    last = 0
    size = len(seq)
    # enumerate over the sequence
    for x, c in enumerate(seq):

        # if the base is a gapped char
        if c in gap_chars:
            # if is_gap False
            if not is_gap:
                # last base is now position of current base
                last = x
                is_gap = True
        else:
            if is_gap:
                yield(last, x)

                last = x
                is_gap = False
    if is_gap:
        yield last, size


def ungapped_regions(seq, gap_chars):
    '''iterator yielding gapped regions in seq.'''
    # check to see if there is a gap at the beginning
    is_gap = seq[0] in gap_chars
    last = 0
    first = True
    size = len(seq)
    # enumerate over the sequence

    if is_gap:
        for x, c in enumerate(seq):

            # if the base is a gapped char
            if c in gap_chars:
                # if is_gap False
                if not is_gap:
                    # last base is now position of current base
                    last = x
                    is_gap = True
            else:
                if is_gap:
                    yield(last, x)

                    last = x
                    is_gap = False
    else:
        start = True
        for x, c in enumerate(seq):

            # if the base is a gapped char
            if c in gap_chars:
                # if is_gap False
                if not is_gap:
                    # last base is now position of current base
                    last = x
                    is_gap = True
            else:
                if is_gap and start:
                    yield(last, 0)
                    start = False
                if is_gap:
                    yield(last, x)

                    last = x
                    is_gap = False

    if is_gap:
        yield last, size


def segmentGaps(infile, gap_char):

    iterator = FastaIterator.FastaIterator(infile)

    while 1:
        try:
            cur_record = next(iterator)
        except StopIteration:
            break

        if cur_record is None:
            break
        contig = re.sub("\s.*", "", cur_record.title)

        for start, end in gapped_regions(cur_record.sequence, gap_char):
            yield(contig, start, end, 0)


def segmentUngapped(infile, gap_char, min_gap_size=0):

    iterator = FastaIterator.FastaIterator(infile)

    while 1:
        try:
            cur_record = next(iterator)
        except StopIteration:
            break

        if cur_record is None:
            break
        contig = re.sub("\s.*", "", cur_record.title)
        size = len(cur_record.sequence)

        last_end = 0
        for start, end in ungapped_regions(cur_record.sequence, gap_char):

            if end == 0:
                yield(contig, end, start, 0)

            if end - start < min_gap_size:
                continue

            if last_end != 0:
                yield(contig, last_end, start, 0)
            last_end = end

        if last_end < size:
            yield(contig, last_end, size, 0)


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
        "-m", "--method", dest="method", type=str,
        choices=(
            "fixed-width-windows-gc",
            "cpg",
            "windows-cpg",
            "gaps",
            "ungapped",
            "windows"),
        help="Method to use for segmentation")

    parser.add_argument(
        "-w", "--window-size=", dest="window_size",
        type=int,
        help="window size for fixed-width windows.")

    parser.add_argument(
        "-s", "--window-shift=", dest="window_shift", type=int,
        help="shift size fixed-width windows.")

    parser.add_argument(
        "--min-cpg", dest="min_cpg", type=int,
        help="minimum number of CpG for windows-cpg")

    parser.add_argument(
        "--min-interval-length", dest="min_length", type=int,
        help="minimum length for ungapped regions")

    parser.set_defaults(
        window_size=10000,
        method="cpg",
        gap_char="NnXx",
        min_length=0,
        window_shift=10000,
        min_cpg=1,
    )

    # add common options (-h/--help, ...) and parse command line
    (args, unknown) = E.start(parser,
                              argv=argv,
                              unknowns=True)

    if args.method == "cpg":
        segments = segmentWithCpG(args.stdin)
    elif args.method == "windows-cpg":
        segments = segmentWindowsCpG(args.stdin,
                                     args.window_size,
                                     args.min_cpg)
    elif args.method == "Isoplotter":
        segments = segmentWithIsoplotter(args.stdin, args)
    elif args.method == "fixed-width-windows-gc":
        segments = segmentFixedWidthWindows(args.stdin,
                                            args.window_size,
                                            args.window_shift,
                                            )
    elif args.method == "gaps":
        segments = segmentGaps(args.stdin, args.gap_char)
    elif args.method == "ungapped":
        segments = segmentUngapped(
            args.stdin, args.gap_char, args.min_length)
    else:
        raise ValueError("unknown method %s" % (method))
    x = 0
    for contig, start, end, gc in segments:
        x += 1
        args.stdout.write("%s\n" % "\t".join(
            (contig, str(start), str(end), str(x), "%6.4f" % (100.0 * gc))))

    # write footer and output benchmark information.
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
