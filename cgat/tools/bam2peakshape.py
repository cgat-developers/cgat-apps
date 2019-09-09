'''bam2peakshape.py - compute peak shape features from a bam-file
==============================================================

:Tags: Genomics NGS Intervals BAM BED Summary

Purpose
-------

This script takes a :term:`bed` formatted file with regions of
interest, for example binding intervals from a ChIP-Seq
experiment. Using a collection of aligned reads is a :term:`bam`
formatted file or :term:`bigwig` formatted file, the script outputs a
collection of features describing the peak shape.

This script is designed with a slight emphasis on ChIP-Seq datasets.
The main reason that this script is better suited for ChIP-Seq is
that(1) it is able to center the counting window at the summit of
every individual peak; (2) it is also able to use the control ChIP-Seq
library to enable side-by-side comparison of treatment vs control;(3)
it can randomly shift the set of input regions to generate a
artificial set of regions, in the absence of real ChIP-Seq control
library, the random regions can provide a peaks profile that can be
used as the control.

For example, given the peaks regions defined by analyzing some
ChIP-Seq dataset (e.g. by using MACS), and without the need to use any
additional genomic annotations (e.g. ENSEMBL, refseq), we can
visualise the binding profiles of transcriptionfactors ChIP-Seq data
relative to the center of each peak regions.

The script outputs a tab-separated table on stdout containing features
for each interval. A peak is defined as the location of the highest
density in an interval. The width of the peak (peak_width) is defined
as the region around the peak in which the density does not drop below
a threshold of peak_heigt * 90%.

Usage
-----

Detailed usage example
++++++++++++++++++++++

The following command will generate the peak shape plot for the peak
regions defined in :file:`onepeak.bed`, using the reads stored in
:file:`small.bam`.  The command will also create a profile for the
control library.  The control library in this example is re-using the
same reads file :file:`small.bam`, however, in your actual experiment,
it should be a different library (the input library for this ChIP-Seq
experiment).::

    python ./scripts/bam2peakshape.py \
        ./tests/bam2peakshape.py/small.bam \
        ./tests/bam2peakshape.py/onepeak.bed \
        --control-bam-file=./tests/bam2peakshape.py/small.bam \
        --use-interval \
        --normalize-transcript


Output files
++++++++++++

Among the features output are:

+-------------------+---------------------------------------------------------+
|*Column*           |*Content*                                                |
+-------------------+---------------------------------------------------------+
|peak_height        |number of reads at peak                                  |
+-------------------+---------------------------------------------------------+
|peak_median        |median coverage compared to peak height                  |
+-------------------+---------------------------------------------------------+
|interval_width     |width of interval                                        |
+-------------------+---------------------------------------------------------+
|peak_width         |width of peak                                            |
+-------------------+---------------------------------------------------------+
|bins               |bins for a histogram of densities within the interval.   |
+-------------------+---------------------------------------------------------+
|npeaks             |number of density peaks in interval.                     |
+-------------------+---------------------------------------------------------+
|peak_center        |point of highest density in interval                     |
+-------------------+---------------------------------------------------------+
|peak_relative_pos  |point of highest density in interval coordinates         |
+-------------------+---------------------------------------------------------+
|counts             |counts for a histogram of densities within the interval  |
+-------------------+---------------------------------------------------------+
|furthest_half_heigh|Distance of peak center to furthest half-height position |
+-------------------+---------------------------------------------------------+
|closest_half_height|Distance of peak center to closest half-height position  |
+-------------------+---------------------------------------------------------+


Additionally, the script outputs a set of matrixes with densities over
intervals that can be used for plotting. The default filenames are
``(matrix|control)_<sortorder>.tsv.gz``, The names can be controlled
with the ``--output-filename-pattern`` option.


Type::

   python bam2peakshape.py --help

for command line help.


Options
-------

Option: Shift
+++++++++++++

shift the each read by a certain distance, because in a ChIP-Seq
experment, the read is always at the edge of an sonicated fragment,
the actual binding site is usually L/2 distance away from the read,
where L is the length of sonicated fragment (determined either
experimentally or computationally).

This option is used only if the input reads are in :term:`bam` formatted file.
If input reads are :term:`bigwig` formatted file, this option is ignored.

Option: Random shift
++++++++++++++++++++

randomly shift the set of input regions to generate a artificial set
of regions. In the absence of real ChIP-Seq control library, the
random regions can provide a peaks profile that can be used as the
control.

Option: Centring method
+++++++++++++++++++++++

"reads" will output in the way that the summit of the peaks are
aligned. "middle" will output in the way that the middle of the input
bed intervals are aligned.

Option: Only interval
+++++++++++++++++++++

Only count reads that are in the interval as defined by the input bed file.

Option: normalization=sum
+++++++++++++++++++++++++

normalize counts such that the sum of all counts in all features are
exactly 1000000.

The detail normalization algorithm as follows: norm = sum(all counts
in all features)/1000000.0 normalized count = normalized count / norm

.. todo::

   paired-endedness is not fully implemented.

Command line options
--------------------

'''

import sys
import os
import re
import cgatcore.experiment as E
import cgatcore.iotools as iotools
import pysam
import cgat.Bed as Bed
import numpy
import collections
import pyBigWig

import cgat.BamTools.peakshape as bam2peakshape


def buildArgumentParser(argv):

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("-f", "--format", dest="format", type=str,
                        choices=("bam", "bigwig"),
                        help="format of genomic input files for densities "
                        )

    parser.add_argument(
        "-o", "--use-interval", dest="use_interval", action="store_true",
        help="only count tags that are in interval given "
        "in bed file. Otherwise, use a fixed width window (see --window-size) "
        "around peak ")

    parser.add_argument(
        "-w", "--window-size", dest="window_size", type=int,
        help="window size in bp on either side of a peak used for getting "
        "read densities. If ``--window-size`` is 1000, the actual window size"
        "will be 2kb, 1kb on either side of the peak in an interval"
        )

    parser.add_argument(
        "-b", "--bin-size", dest="bin_size", type=int,
        help="bin-size in bp for computing read densities. "
        "If ``--window-size`` is set to 1000 and ``--bin-size`` to 10, "
        "there will be 100 bins on either side of a peak. "
        )

    parser.add_argument(
        "--smooth-method", dest="smooth_method", type=str,
        choices=("none", "sum", "sg"),
        help="smooting method to apply to density data before sampling "
        "according to ``bin-size``. sg=SavitzkyGolay, sum=sum density in bin, "
        "none=no smoothing "
        )

    parser.add_argument("-s", "--sort-order", dest="sort_orders",
                        type=str,
                        action="append",
                        choices=("peak-height", "peak-width", "unsorted",
                                 "interval-width", "interval-score"),
                        help="output sort order for matrices. "
                        )

    parser.add_argument(
        "-c", "--control-bam-file", "--control-bigwig-file",
        action="append",
        dest="control_files",
        type=str,
        help="control file. If given, two peakshapes are computed, "
        "one for the primary data and one for the control data. "
        "The control file is centered around the same "
        "base as the primary file and output in the same "
        "sort order as the primary profile to all side-by-side. "
        "comparisons. Multiple control files can be given. The "
        "control files should have the same format as the "
        "principal input file "
        )

    parser.add_argument(
        "-r", "--random-shift", dest="random_shift", action="store_true",
        help="shift intervals in random direction up/downstream of interval "
        )

    parser.add_argument(
        "-e", "--centring-method", dest="centring_method", type=str,
        choices=("reads", "middle"),
        help="centring method. Available are: "
        "reads=use density to determine peak, "
        "middle=use middle of interval "
        )

    parser.add_argument(
        "-n", "--normalize-matrix", dest="normalization", type=str,
        choices=("none", "sum"),
        help="matrix normalisation to perform. "
        )

    parser.add_argument(
        "--use-strand", dest="strand_specific", action="store_true",
        help="use strand information in intervals. Intervals on the "
        "negative strand are flipped "
        )

    parser.add_argument(
        "-i", "--shift-size", dest="shift", type=int,
        help="shift for reads. When processing bam files, "
        "reads will be shifted upstream/downstream by this amount. "
        )

    parser.set_defaults(
        bin_size=10,
        shift=0,
        window_size=1000,
        sort_orders=[],
        centring_method="reads",
        control_files=[],
        random_shift=False,
        strand_specific=False,
        format="bam",
        report_step=100,
        use_interval=False,
        smooth_method=None,
    )

    return parser


IntervalData = collections.namedtuple(
    "IntervalData",
    "foreground interval controls shifted")


def outputFeatureTable(outfile, features_per_interval, bins):
    '''ouput results from density profiles.'''

    outfile.write("\t".join(
        ("contig",
         "start",
         "end",
         "name",
         "\t".join(bam2peakshape.PeakShapeResult._fields))) + "\n")

    # output principal table
    n = 0
    for foreground, bed, controls, shifted in features_per_interval:
        n += 1
        if "name" in bed:
            name = bed.name
        else:
            name = str(n)
        outfile.write("%s\t%i\t%i\t%s\t" %
                      (bed.contig, bed.start, bed.end, name))

        outfile.write("\t".join(map(str, foreground[:-2])))
        bins, counts = foreground[-2], foreground[-1]
        outfile.write("\t%s" % ",".join(map(str, bins)))
        outfile.write("\t%s" % ",".join(map(str, counts)))
        outfile.write("\n")


def writeMatricesForSortOrder(features_per_interval,
                              bins,
                              foreground_track,
                              control_tracks,
                              shifted,
                              sort_order):
    '''output one or more matrices for each sort sorder.

    For each sort order output the forerground. If there
    are additional controls and shifted section, output
    these as well

    The files will named:
    matrix_<track>_<sortorder>

    '''
    if "name" in features_per_interval[0].interval:
        names = [x.interval.name for x in features_per_interval]
    else:
        names = list(map(str, list(range(1, len(features_per_interval) + 1))))

    bins = ["%i" % x for x in bins]
    sort_order = re.sub("-", "_", sort_order)

    # write foreground
    iotools.write_matrix(
        E.open_output_file("matrix_%s_%s.gz" % (foreground_track, sort_order)),
        [x.foreground.counts for x in features_per_interval],
        row_headers=names,
        col_headers=bins,
        row_header="name")

    # write controls
    for idx, track in enumerate(control_tracks):
        iotools.write_matrix(
            E.open_output_file("matrix_%s_%s.gz" % (track, sort_order)),
            [x.controls[idx].counts for x in features_per_interval],
            row_headers=names,
            col_headers=bins,
            row_header="name")

    # write shifted matrix
    if shifted:
        iotools.write_matrix(
            E.open_output_file("matrix_shift_%s.gz" % (sort_order)),
            [x.shifted.counts for x in features_per_interval],
            row_headers=names,
            col_headers=bins,
            row_header="name")

    # output a combined matrix
    if len(control_tracks) > 0 or shifted:
        rows = []
        for row in features_per_interval:
            l = [row.foreground.counts]
            l.extend([row.controls[x].counts for x in
                      range(len(control_tracks))])
            if shifted:
                l.append(row.shifted.counts)
            rows.append(numpy.concatenate(l))

        n = 1 + len(control_tracks)
        if shifted:
            n += 1

        # make column names unique and make sure they can be sorted
        # lexicographically
        all_bins = []
        for x in range(n):
            all_bins.extend(["%i:%s" % (x, b) for b in bins])

        iotools.write_matrix(
            E.open_output_file("matrix_sidebyside_%s.gz" % (sort_order)),
            rows,
            row_headers=names,
            col_headers=all_bins,
            row_header="name")


def outputMatrices(features_per_interval,
                   bins,
                   foreground_track,
                   control_tracks=None,
                   shifted=False,
                   sort_orders=None):
    '''ouput matrices from density profiles
    in one or more sort_orders.
    '''

    # output sorted matrices
    if not sort_orders:
        writeMatricesForSortOrder(features_per_interval,
                                  bins,
                                  foreground_track,
                                  control_tracks,
                                  shifted,
                                  "unsorted")

    for sort_order in sort_orders:

        if sort_order == "peak-height":
            features_per_interval.sort(
                key=lambda x: x.foreground.peak_height)

        elif sort_order == "peak-width":
            features_per_interval.sort(
                key=lambda x: x.foreground.peak_width)

        elif sort_order == "interval-width":
            features_per_interval.sort(
                key=lambda x: x.interval.end - x.interval.start)

        elif sort_order == "interval-score":
            try:
                features_per_interval.sort(
                    key=lambda x: float(x.interval.score))
            except IndexError:
                E.warn("score field not present - no output")
                continue
            except TypeError:
                E.warn("score field not a valid number - no output")
                continue

        writeMatricesForSortOrder(features_per_interval,
                                  bins,
                                  foreground_track,
                                  control_tracks,
                                  shifted,
                                  sort_order)


def buildDensityMatrices(bedfile,
                         fg_file,
                         control_files,
                         counter,
                         window_size=1000,
                         bin_size=10,
                         strand_specific=False,
                         centring_method="reads",
                         use_interval=False,
                         random_shift=False,
                         smooth_method="none",
                         report_step=1000):
    '''compute densities and peakshape parameters
    in intervals given by *bedfile* using reads in *fg_file*.

    If *control_files* are given, densities are produced for
    these as well.

    Returns a list of results for each interval in *bedfile* of
    type IntervalData and an array of bin-values.
    '''

    if window_size:
        # bins are centered at peak-center and then stretching outwards.
        bins = numpy.arange(-window_size + bin_size // 2,
                            +window_size,
                            bin_size)

    result = []
    c = E.Counter()
    c.input = 0

    for bed in bedfile:
        c.input += 1

        # if bed.contig not in contigs:
        #    c.skipped += 1
        #    continue

        if c.input % report_step == 0:
            E.info("iteration: %i" % c.input)

        features = counter.countInInterval(
            fg_file,
            bed.contig, bed.start, bed.end,
            window_size=window_size,
            bins=bins,
            use_interval=use_interval,
            centring_method=centring_method)

        if features is None:
            c.skipped += 1
            continue

        if control_files:
            control = []
            for control_file in control_files:
                control.append(counter.countAroundPos(
                    control_file,
                    bed.contig,
                    features.peak_center,
                    bins=features.bins))

        else:
            control = None

        if random_shift:
            direction = numpy.random.randint(0, 2)
            if direction:
                pos = features.peak_center + 2 * bins[0]
            else:
                pos = features.peak_center + 2 * bins[-1]
            shifted = counter.countAroundPos(fg_file,
                                             bed.contig,
                                             pos,
                                             bins=features.bins)
        else:
            shifted = None

        if strand_specific and bed.strand == "-":
            features._replace(hist=features.hist[::-1])
            if control:
                for c in control:
                    c._replace(hist=c.hist[::-1])
            if shifted:
                shifted._replace(hist=shifted.hist[::-1])

        result.append(IntervalData._make((features, bed, control, shifted)))
        c.added += 1

    E.info("interval processing: %s" % c)

    return result, bins


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    parser = buildArgumentParser(argv)

    # add common options (-h/--help, ...) and parse command line
    (args, unknown) = E.start(parser, argv=argv, add_output_options=True, unknowns=True)

    if len(unknown) != 2:
        raise ValueError(
            "please specify one bam- or wig-file and one bed file")

    if args.control_files:
        E.info("using control files: %s" % ",".join(args.control_files))

    infile, bedfile = unknown
    control_files = []

    if args.format == "bigwig":
        fg_file = pyBigWig.open(infile)
        for control_file in args.control_files:
            control_files.append(pyBigWig.open(control_file))
        counter = bam2peakshape.CounterBigwig(
            smooth_method=args.smooth_method)

    elif args.format == "bam":
        fg_file = pysam.AlignmentFile(infile, "rb")
        for control_file in args.control_files:
            control_files.append(pysam.AlignmentFile(control_file, "rb"))
        counter = bam2peakshape.CounterBam(
            shift=args.shift,
            smooth_method=args.smooth_method)

    features_per_interval, bins = buildDensityMatrices(
        Bed.iterator(iotools.open_file(bedfile)),
        fg_file,
        control_files,
        counter,
        window_size=args.window_size,
        bin_size=args.bin_size,
        strand_specific=args.strand_specific,
        centring_method=args.centring_method,
        use_interval=args.use_interval,
        random_shift=args.random_shift,
        smooth_method=args.smooth_method,
        report_step=args.report_step)

    if len(features_per_interval) == 0:
        E.warn("no data - no output")
        E.stop()
        return

    outputFeatureTable(args.stdout, features_per_interval, bins)

    # apply normalization
    # Note: does not normalize control?
    # Needs reworking, currently it does not normalize across
    # all samples nor does the work "sum" reflect the per million
    # normalization.
    if args.normalization == "sum":
        E.info("starting sum normalization")
        # get total counts across all intervals
        norm = 0.0
        for foreground, bed, controls, shifted in features_per_interval:
            norm += sum(foreground.counts)
        # per million
        norm /= float(1000000)
        E.info("sum/million normalization with %f" % norm)

        # normalise
        new_data = []
        for foreground, bed, controls, shifted in features_per_interval:
            foreground = foreground._replace(
                counts=numpy.array(foreground.counts,
                                   dtype=numpy.float) / norm)
            new_controls = []
            for control in controls:
                new_controls.append(
                    control._replace(
                        counts=numpy.array(control.counts,
                                           dtype=numpy.float) / norm))
            if shifted:
                shifted = shifted._replace(
                    counts=numpy.array(shifted.counts,
                                       dtype=numpy.float) / norm)
            new_data.append(IntervalData._make((
                foreground, bed, new_controls, shifted)))
        features_per_interval = new_data
    else:
        E.info("no normalization performed")

    # center bins
    out_bins = bins[:-1] + args.bin_size

    # build tracks
    def _toTrack(filename):
        return os.path.splitext(os.path.basename(filename))[0]

    outputMatrices(features_per_interval,
                   out_bins,
                   foreground_track=_toTrack(infile),
                   control_tracks=[_toTrack(x) for x in args.control_files],
                   shifted=args.random_shift,
                   sort_orders=args.sort_orders)

    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
