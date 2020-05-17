'''bam2geneprofile.py - build meta-gene profile for a set of transcripts/genes
===========================================================================

:Tags: Genomics NGS Genesets Intervals GTF BAM Summary

Purpose
-------

This script takes a :term:`gtf` formatted file, a short reads
:term:`bam` formatted file and computes meta-gene profiles over
various annotations derived from the :term:`gtf` file.

A meta-gene profile is an abstract genomic entity over which reads
stored in a :term:`bam` formatted file have been counted. A meta-gene
might be an idealized eukaryotic gene (upstream, exonic sequence,
downstream) or any other genomic landmark of interest such as
transcription start sites.

The script can be used to visualize binding profiles of a chromatin
mark in gene bodies, binding of transcription factors in promotors or
sequencing bias (e.g. 3' bias) in RNA-Seq data.

This script is designed with a slight emphasis on RNA-Seq
datasets. For example, it takes care of spliced reads, by using the
CIGAR string in the BAM file to accurately define aligned bases (if
the --base-accurate is specified, currently this feature is turned off
by default for speed reasons).

Alternatively, for the purpose of visualizing binding profiles of
transcription factor ChIP-Seq without the need to use any genomic
annotations (ENSEMBL, or refseq), you may also consider using
:doc:`bam2peakshape`, which is designed with a slight emphasis on
Chip-Seq datasets. For example, :doc:`bam2peakshape` is able to center
the counting window to the summit of every individual peak.
:doc:`bam2peakshape` is also able to: (1) plot the control ChIP-Seq
library to enable side-by-side comparison; (2) randomize the given
regions to provide a semi-control.

Usage
-----

Quick start examples
++++++++++++++++++++

The following command will generate the gene profile plot similar to
Fig 1(a) in the published cgat paper, but using a test dataset that is
much smaller and simpler than the dataset used for publishing the cgat
paper. ::

    python ./scripts/bam2geneprofile.py
        --bam-file=./tests/bam2geneprofile.py/multipleReadsSplicedOutAllIntronsAndSecondExon.bam
        --gtf-file=./tests/bam2geneprofile.py/onegeneWithoutAnyCDS.gtf.gz
        --method=geneprofile
        --reporter=gene

In the following, a slightly more involved example will use more
features of this script. The following command generate the gene
profile showing base accuracy of upstream (500bp), exons, introns and
downstream(500bp) of a gene model from some user supplied RNA-Seq data
and geneset. ::

    python ./scripts/bam2geneprofile.py
        --bam-file=./rnaseq.bam
        --gtf-file=./geneset.gtf.gz
        --method=geneprofilewithintrons
        --reporter=gene
        --extension-upstream=500
        --resolution-upstream=500
        --extension-downstream=500
        --resolution-downstream=500

The output will contain read coverage over genes. The profile will
contain four separate segments:

1. the upstream region of a gene ( set to be 500bp ),
   (``--extension-upstream=500``).

2. the transcribed region of a gene. The transcribed region of every gene will
   be scaled to 1000 bp (default), shrinking longer transcripts and
   expanding shorter transcripts.

3. the intronic regions of a gene. These will be scaled to 1000b (default).

4. the downstream region of a gene (set to be 500bp),
   (``--extension-downstream=500``).



Detailed explaination
+++++++++++++++++++++

The :file:`bam2geneprofile.py` script reads in a set of transcripts
from a :term:`gtf` formatted file. For each transcript, overlapping
reads from the provided :term:`bam` file are collected. The counts
within the transcript are then mapped onto the meta-gene structure and
counts are aggregated over all transcripts in the :term:`gtf` file.

:term:`Bam` files need to be sorted by coordinate and indexed.

A meta-gene structure has two components - regions of variable size,
such as exons, introns, etc, which nevertheless have a fixed start and
end coordinate in a transcript. The other component are regions of
fixed width, such a regions of a certain size upstream or downstream
of a landmark such as a transcription start site.

The size of the former class, regions of variable size, can be varied
with ``--resolution`` options. For example, the option
``--resolution-upstream-utr=1000`` will create a meta-gene with a
1000bp upstream UTR region. UTRs that are larger will be compressed,
and UTRs that are smaller, will be stretched to fit the 1000bp
meta-gene UTR region.

The size of fixed-width regions can be set with ``--extension``
options. For example, the options ``--extension-upstream`` will set
the size of the uptsream extension region to 1000bp. Note that no
scaling is required when counting reads towards the fixed-width
meta-gene profile.

Type::

   python bam2geneprofile.py --help

for command line help.

Options
-------

The script provides a variety of different meta-gene structures i.e.
geneprofiles, selectable via using the option: (``--method``).

Profiles
++++++++

Different profiles are accessible through the ``--method`` option. Multiple
methods can be applied at the same time. While ``upstream`` and ``downstream``
typically have a fixed size, the other regions such as ``CDS``, ``UTR`` will be
scaled to a common size.

utrprofile
    UPSTREAM - UTR5 - CDS - UTR3 - DOWNSTREAM
    gene models with UTR. Separate the coding section from the non-coding part.

geneprofile
    UPSTREAM - EXON - DOWNSTREAM
    simple exonic gene models

geneprofilewithintrons
    UPSTREAM - EXON - INTRON - DOWNSTREAM

    gene models containing also intronic sequence, only correct if
    used with ``--use-base-accuracy`` option.

separateexonprofile
    UPSTREAM - FIRST EXON - EXON - LAST EXON - DOWNSTREAM

    gene models with the first and last exons separated out from all
    other exons.  Only applicable to gene models with > 1 exons.

separateexonprofilewithintrons
    UPSTREAM - FIRST EXON - EXON - INTRON - LAST EXON - DOWNSTREAM

    gene models with first and last exons separated out, and includes
    all introns together.  Excludes genes with < 2 exons and no introns.

geneprofileabsolutedistancefromthreeprimeend

    UPSTREAM - EXON (absolute distance, see below) - INTRON (absolute
    distance, see below) - DOWNSTREAM (the downstream of the exons)
    region, the script counts over the mRNA transcript only, skipping
    introns. Designed to visualize the 3 prime bias in RNASeq data,
    only correct if used together with ``--use-base-accuracy`` option.

    absolute distance: In order to to visualize the 3 prime bias,
    genes are not supposed to be streched to equal length as it did in
    all other counting methods. In this counting method, we first set
    a fix length using
    ``--extension-exons-absolute-distance-topolya``, the script will
    discard genes shorter than this fixed length. For genes (when all
    the exons stitched together) longer than this fixed length, the
    script will only count over this fixed length ( a absolute
    distance ) from three prime end, instead of compress the longer
    genes. Same goes for absolute distance intron counting.

tssprofile
    UPSTREAM - DOWNSTREAM
    transcription start/stop sites

intervalprofile
    UPSTREAM - INTERVAL - DOWNSTREAM
    Similar to geneprofile, but count over the complete span of the gene
    (including introns).

midpointprofile
    UPSTREAM  - DOWNSTREAM
    aggregate over midpoint of gene model


Normalization
+++++++++++++

Normalization can be applied in two stages of the computation.

Count vector normalization
^^^^^^^^^^^^^^^^^^^^^^^^^^

Before adding counts to the meta-gene profile, the profile for the
individual transcript can be normalized. Without normalization, highly
expressed genes will contribute more to the meta-gene profile than
lowly expressed genes.  Normalization can assure that each gene
contributes an equal amount.

Normalization is applied to the vector of read counts that is computed
for each transcript. Normalization can be applied for the whole
transcript (``total``) or on a per segment basis depending on the
counter. For example, in the gene counter, exons, upstream and
downstream segments can be normalized independently.

Counts can be normalized either by the maximum or the sum of all
counts in a segment or across the whole transcript. Normalization is
controlled with the command line option ``--normalize-trancript``. Its
arguments are:

* ``none``: no normalization
* ``sum``: sum of counts within a region (exons, upstream, ...).
  The area under the curve will sum to 1 for each region.
* ``max``: maximum count within a region (exons,upstream, ...).
* ``total-sum``: sum of counts across all regions. The area
  under the curve will sum to 1 for
  the complete transcript.
* ``total-max``: maximum count across all regions.

The options above control the contribution of individual transcripts
to a meta-gene profile and are thus suited for example for RNA-Seq data.

The options above do not control for different read-depths or any
local biases. To compare meta-gene profiles between samples,
additional normalization is required.

Meta-gene profile normalization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To enable comparison between experiments, the meta-gene profile itself
can be normalized.  Normalization a profile can help comparing the
shapes of profiles between different experiments independent of the
number of reads or transcripts used in the construction of the
meta-gene profile.

Meta-gene profile normalization is controlled via the
``--normalize-profile`` option. Possible normalization are:

* none: no normalization
* area: normalize such that the area under the meta-gene profile is 1.
* counts: normalize by number of features (genes,tss) that have been counted.
* background: normalize with background (see below).

A special normalization is activated with the ``background`` option.
Here, the counts at the left and right most regions are used to
estimate a background level for each transcript. The counts are then
divided by this background-level. The assumption is that the meta-gene
model is computed over a large enough area to include genomic
background.

Genes versus transcripts
++++++++++++++++++++++++

The default is to collect reads on a per-transcript
level. Alternatively, the script can merge all transcripts of a gene
into a single virtual transcript. Note that this virtual transcript
might not be a biologically plausible transcript. It is usually better
to provide :file:`bam2geneprofile.py` with a set of representative
transcripts per gene in order to avoid up-weighting genes with
multiple transcripts.

Control
+++++++

If control files (chip-seq input tracks) are supplied, counts in the
control file can be used to compute a fold-change.

Bed and wiggle files
++++++++++++++++++++

The densities can be computed from :term:`bed` or :term:`wiggle`
formatted files. If a :term:`bed` formatted file is supplied, it must
be compressed with and indexed with :file:`tabix`.

.. note::

   Paired-endedness is ignored. Both ends of a paired-ended read are
   treated individually.


Command line options
--------------------

'''

import os
import sys
import cgatcore.experiment as E
import cgatcore.iotools as iotools
import pysam
import cgat.GTF as GTF
import numpy
import pandas
import pyBigWig

from cgat.BamTools import geneprofile


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
        """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("-m", "--method", dest="methods", type=str,
                        action="append",
                        choices=("geneprofile", "tssprofile", "utrprofile",
                                 "intervalprofile", "midpointprofile",
                                 "geneprofilewithintrons",
                                 "geneprofileabsolutedistancefromthreeprimeend",
                                 "separateexonprofile",
                                 "separateexonprofilewithintrons",
                                 ),
                        help='counters to use. Counters describe the '
                        'meta-gene structure to use. '
                        'Note using geneprofilewithintrons, or '
                        'geneprofileabsolutedistancefromthreeprimeend will '
                        'automatically turn on the --use-base-accuracy option')

    parser.add_argument("-b", "--bam-file", "--bedfile", "--bigwigfile",
                        dest="infiles",
                        metavar="BAM",
                        type=str, action="append",
                        help="BAM/bed/bigwig files to use. Do not mix "
                        "different types ")

    parser.add_argument("-c", "--control-bam-file", dest="controlfiles",
                        metavar="BAM",
                        type=str, action="append",
                        help="control/input to use. Should be of the same "
                        "type as the bam/bed/bigwig file")

    parser.add_argument("-g", "--gtf-file", dest="gtffile", type=str,
                        metavar="GTF",
                        help="GTF file to use. ")

    parser.add_argument(
        "--normalize-transcript",
        dest="transcript_normalization",
        type=str,
        choices=("none", "max", "sum", "total-max", "total-sum"),
        help="normalization to apply on each transcript "
        "profile before adding to meta-gene profile. "
        )

    parser.add_argument(
        "--normalize-profile",
        dest="profile_normalizations",
        type=str, action="append",
        choices=("all", "none", "area", "counts", "background"),
        help="normalization to apply on meta-gene "
        "profile normalization. "
        )

    parser.add_argument(
        "-r", "--reporter", dest="reporter", type=str,
        choices=("gene", "transcript"),
        help="report results for genes or transcripts."
        " When 'genes` is chosen, exons across all transcripts for"
        " a gene are merged. When 'transcript' is chosen, counts are"
        " computed for each transcript separately with each transcript"
        " contributing equally to the meta-gene profile."
        " ")

    parser.add_argument("-i", "--shift-size", dest="shifts", type=int,
                        action="append",
                        help="shift reads in :term:`bam` formatted file "
                        "before computing densities (ChIP-Seq). "
                        )

    parser.add_argument("-a", "--merge-pairs", dest="merge_pairs",
                        action="store_true",
                        help="merge pairs in :term:`bam` formatted "
                        "file before computing "
                        "densities (ChIP-Seq). "
                        )

    parser.add_argument("-u", "--use-base-accuracy", dest="base_accuracy",
                        action="store_true",
                        help="compute densities with base accuracy. The default "
                        "is to only use the start and end of the aligned region "
                        "(RNA-Seq) "
                        )

    parser.add_argument("-e", "--extend", dest="extends", type=int,
                        action="append",
                        help="extend reads in :term:`bam` formatted file "
                        "(ChIP-Seq). "
                        )

    parser.add_argument("--resolution-upstream", dest="resolution_upstream",
                        type=int,
                        help="resolution of upstream region in bp "
                        )

    parser.add_argument("--resolution-downstream", dest="resolution_downstream",
                        type=int,
                        help="resolution of downstream region in bp "
                        )

    parser.add_argument("--resolution-upstream-utr",
                        dest="resolution_upstream_utr",
                        type=int,
                        help="resolution of upstream UTR region in bp "
                        )

    parser.add_argument("--resolution-downstream-utr",
                        dest="resolution_downstream_utr",
                        type=int,
                        help="resolution of downstream UTR region in bp "
                        )

    parser.add_argument("--resolution-cds", dest="resolution_cds", type=int,
                        help="resolution of cds region in bp "
                        )

    parser.add_argument("--resolution-first-exon", dest="resolution_first",
                        type=int,
                        help="resolution of first exon in gene, in bp"
                        )

    parser.add_argument("--resolution-last-exon", dest="resolution_last",
                        type=int,
                        help="resolution of last exon in gene, in bp"
                        )

    parser.add_argument("--resolution-introns",
                        dest="resolution_introns", type=int,
                        help="resolution of introns region in bp "
                        )

    parser.add_argument("--resolution-exons-absolute-distance-topolya",
                        dest="resolution_exons_absolute_distance_topolya",
                        type=int,
                        help="resolution of exons absolute distance "
                        "topolya in bp "
                        )

    parser.add_argument("--resolution-introns-absolute-distance-topolya",
                        dest="resolution_introns_absolute_distance_topolya",
                        type=int,
                        help="resolution of introns absolute distance "
                        "topolya in bp "
                        )

    parser.add_argument("--extension-exons-absolute-distance-topolya",
                        dest="extension_exons_absolute_distance_topolya",
                        type=int,
                        help="extension for exons from the absolute "
                        "distance from the topolya in bp "
                        )

    parser.add_argument(
        "--extension-introns-absolute-distance-topolya",
        dest="extension_introns_absolute_distance_topolya", type=int,
        help="extension for introns from the absolute distance from "
        "the topolya in bp ")

    parser.add_argument(
        "--extension-upstream", dest="extension_upstream", type=int,
        help="extension upstream from the first exon in bp"
        )

    parser.add_argument(
        "--extension-downstream", dest="extension_downstream", type=int,
        help="extension downstream from the last exon in bp"
        )

    parser.add_argument(
        "--extension-inward", dest="extension_inward", type=int,
        help="extension inward from a TSS start site in bp"
        )

    parser.add_argument(
        "--extension-outward", dest="extension_outward", type=int,
        help="extension outward from a TSS start site in bp"
        )

    parser.add_argument("--scale-flank-length", dest="scale_flanks", type=int,
                        help="scale flanks to (integer multiples of) gene length"
                        )

    parser.add_argument(
        "--control-factor", dest="control_factor", type=float,
        help="factor for normalizing control and foreground data. "
        "Computed from data if not set. "
        )

    parser.add_argument("--output-all-profiles", dest="output_all_profiles",
                        action="store_true",
                        help="keep individual profiles for each "
                        "transcript and output. "
                        )

    parser.add_argument("--counts-tsv-file", dest="input_filename_counts",
                        type=str,
                        help="filename with count data for each transcript. "
                        "Use this instead "
                        "of recomputing the profile. Useful for plotting the "
                        "meta-gene profile "
                        "from previously computed counts "
                        )

    parser.add_argument(
        "--background-region-bins",
        dest="background_region_bins",
        type=int,
        help="number of bins on either end of the profile "
        "to be considered for background meta-gene normalization "
        )

    parser.add_argument("--output-res",
                        dest="resolution_images", type=int,
                        help="the output dpi for the figure plot - will default to "
                        )

    parser.add_argument("--image-format", dest="image_format", type=str,
                        help="The output format for the figure plot - defaults to "
                        )

    parser.set_defaults(
        remove_rna=False,
        ignore_pairs=False,
        force_output=False,
        bin_size=10,
        extends=[],
        shifts=[],
        sort=[],
        reporter="transcript",
        resolution_cds=1000,
        resolution_introns=1000,
        # 3kb is a good balance of seeing long enough 3 prime bias and not omit
        # too many genes. Tim 31th Aug 2013
        resolution_exons_absolute_distance_topolya=3000,
        # introns is only for assess the noise level, thus do ont need a long
        # region, a long region has the side effect of omit more genes. Tim
        # 31th Aug 2013
        resolution_introns_absolute_distance_topolya=500,
        # extension can simply just be the same as resolution
        extension_exons_absolute_distance_topolya=3000,
        extension_introns_absolute_distance_topolya=500,
        resolution_upstream_utr=1000,
        resolution_downstream_utr=1000,
        resolution_upstream=1000,
        resolution_downstream=1000,
        resolution_first=1000,
        resolution_last=1000,
        # mean length of transcripts: about 2.5 kb
        extension_upstream=2500,
        extension_downstream=2500,
        extension_inward=3000,
        extension_outward=3000,
        plot=True,
        methods=[],
        infiles=[],
        controlfiles=[],
        gtffile=None,
        profile_normalizations=[],
        transcript_normalization=None,
        scale_flanks=0,
        merge_pairs=False,
        min_insert_size=0,
        max_insert_size=1000,
        base_accuracy=False,
        matrix_format="single",
        control_factor=None,
        output_all_profiles=False,
        background_region_bins=10,
        input_filename_counts=None,
        resolution_images=None,
        image_format="png",
    )

    # add common options (-h/--help, ...) and parse command line
    (args, unknown) = E.start(parser, argv=argv, add_output_options=True, unknowns=True)

    # Keep for backwards compatability
    if len(unknown) == 2:
        infile, gtf = unknown
        args.infiles.append(infile)
        args.gtffile = gtf

    if not args.gtffile:
        raise ValueError("no GTF file specified")

    if args.gtffile == "-":
        args.gtffile = args.stdin
    else:
        args.gtffile = iotools.open_file(args.gtffile)

    if len(args.infiles) == 0:
        raise ValueError("no bam/wig/bed files specified")

    for methodsRequiresBaseAccuracy in [
            "geneprofilewithintrons",
            "geneprofileabsolutedistancefromthreeprimeend",
    ]:
        # If you implemented any methods that you do not want the
        # spliced out introns or exons appear to be covered by
        # non-existent reads, it is better you let those methods imply
        # --base-accurarcy by add them here.
        if methodsRequiresBaseAccuracy in args.methods:
            args.base_accuracy = True

    if args.reporter == "gene":
        gtf_iterator = GTF.flat_gene_iterator(GTF.iterator(args.gtffile))
    elif args.reporter == "transcript":
        gtf_iterator = GTF.transcript_iterator(GTF.iterator(args.gtffile))

    # Select rangecounter based on file type
    if len(args.infiles) > 0:
        if args.infiles[0].endswith(".bam"):
            bamfiles = [pysam.AlignmentFile(x, "rb") for x in args.infiles]

            if args.controlfiles:
                controlfiles = [pysam.AlignmentFile(x, "rb")
                                for x in args.controlfiles]
            else:
                controlfiles = None

            format = "bam"
            if args.merge_pairs:
                range_counter = geneprofile.RangeCounterBAM(
                    bamfiles,
                    shifts=args.shifts,
                    extends=args.extends,
                    merge_pairs=args.merge_pairs,
                    min_insert_size=args.min_insert_size,
                    max_insert_size=args.max_insert_size,
                    controfiles=controlfiles,
                    control_factor=args.control_factor)

            elif args.shifts or args.extends:
                range_counter = geneprofile.RangeCounterBAM(
                    bamfiles,
                    shifts=args.shifts,
                    extends=args.extends,
                    controlfiles=controlfiles,
                    control_factor=args.control_factor)

            elif args.base_accuracy:
                range_counter = geneprofile.RangeCounterBAMBaseAccuracy(
                    bamfiles,
                    controlfiles=controlfiles,
                    control_factor=args.control_factor)
            else:
                range_counter = geneprofile.RangeCounterBAM(
                    bamfiles,
                    controlfiles=controlfiles,
                    control_factor=args.control_factor)

        elif args.infiles[0].endswith(".bed.gz"):
            bedfiles = [pysam.Tabixfile(x) for x in args.infiles]

            if args.controlfiles:
                controlfiles = [pysam.Tabixfile(x)
                                for x in args.controlfiles]
            else:
                controlfiles = None

            range_counter = geneprofile.RangeCounterBed(
                bedfiles,
                controlfiles=controlfiles,
                control_factor=args.control_factor)

        elif args.infiles[0].endswith(".bw"):
            wigfiles = [pyBigWig.open(x) for x in args.infiles]
            range_counter = geneprofile.RangeCounterBigWig(wigfiles)

        else:
            raise NotImplementedError(
                "can't determine file type for %s" % str(args.infiles))

    counters = []
    for method in args.methods:
        if method == "utrprofile":
            counters.append(geneprofile.UTRCounter(
                range_counter,
                args.resolution_upstream,
                args.resolution_upstream_utr,
                args.resolution_cds,
                args.resolution_downstream_utr,
                args.resolution_downstream,
                args.extension_upstream,
                args.extension_downstream,
            ))

        elif method == "geneprofile":
            counters.append(geneprofile.GeneCounter(
                range_counter,
                args.resolution_upstream,
                args.resolution_cds,
                args.resolution_downstream,
                args.extension_upstream,
                args.extension_downstream,
                args.scale_flanks))

        elif method == "geneprofilewithintrons":
            counters.append(geneprofile.GeneCounterWithIntrons(
                range_counter,
                args.resolution_upstream,
                args.resolution_cds,
                args.resolution_introns,
                args.resolution_downstream,
                args.extension_upstream,
                args.extension_downstream,
                args.scale_flanks))

        elif method == "geneprofileabsolutedistancefromthreeprimeend":
            # args.extension_exons_absolute_distance_tostartsite,
            # args.extension_introns_absolute_distance_tostartsite,
            # Tim 31th Aug 2013: a possible feature for future,  if five prime
            # bias is of your interest.
            # (you need to create another class). It is not very difficult to
            # derive from this class, but is not implemented yet
            # This future feature is slightly different the TSS profile
            # already implemented, because in this future feature introns are
            # skipped,
            counters.append(
                geneprofile.GeneCounterAbsoluteDistanceFromThreePrimeEnd(
                    range_counter, args.resolution_upstream,
                    args.resolution_downstream,
                    args.resolution_exons_absolute_distance_topolya,
                    args.resolution_introns_absolute_distance_topolya,
                    args.extension_upstream,
                    args.extension_downstream,
                    args.extension_exons_absolute_distance_topolya,
                    args.extension_introns_absolute_distance_topolya,
                    args.scale_flanks))

        elif method == "tssprofile":
            counters.append(geneprofile.TSSCounter(
                range_counter,
                args.extension_outward,
                args.extension_inward))

        elif method == "intervalprofile":
            counters.append(geneprofile.RegionCounter(
                range_counter,
                args.resolution_upstream,
                args.resolution_cds,
                args.resolution_downstream,
                args.extension_upstream,
                args.extension_downstream))

        elif method == "midpointprofile":
            counters.append(geneprofile.MidpointCounter(
                range_counter,
                args.resolution_upstream,
                args.resolution_downstream,
                args.extension_upstream,
                args.extension_downstream))

        # add new method to split 1st and last exons out
        # requires a representative transcript for reach gene
        # gtf should be sorted gene-position
        elif method == "separateexonprofile":
            counters.append(geneprofile.SeparateExonCounter(
                range_counter,
                args.resolution_upstream,
                args.resolution_first,
                args.resolution_last,
                args.resolution_cds,
                args.resolution_downstream,
                args.extension_upstream,
                args.extension_downstream))

        elif method == "separateexonprofilewithintrons":
            counters.append(geneprofile.SeparateExonWithIntronCounter(
                range_counter,
                args.resolution_upstream,
                args.resolution_first,
                args.resolution_last,
                args.resolution_cds,
                args.resolution_introns,
                args.resolution_downstream,
                args.extension_upstream,
                args.extension_downstream))

    # set normalization
    for c in counters:
        c.setNormalization(args.transcript_normalization)
        if args.output_all_profiles:
            c.setOutputProfiles(iotools.open_file(E.get_output_file(c.name) +
                                                  ".profiles.tsv.gz", "w"))

    if args.input_filename_counts:
        # read counts from file
        E.info("reading counts from %s" % args.input_filename_counts)
        all_counts = pandas.read_csv(
            iotools.open_file(args.input_filename_counts),
            sep='\t', header=0, index_col=0)

        if len(counters) != 1:
            raise NotImplementedError(
                'counting from matrix only implemented for 1 counter.')
        # build counter based on reference counter
        counter = geneprofile.UnsegmentedCounter(counters[0])
        counters = [counter]
        geneprofile.countFromCounts(counters, all_counts)

    else:
        E.info("starting counting with %i counters" % len(counters))
        feature_names = geneprofile.countFromGTF(counters,
                                                 gtf_iterator)

    # output matrices
    if not args.profile_normalizations:
        args.profile_normalizations.append("none")
    elif "all" in args.profile_normalizations:
        args.profile_normalizations = ["none",
                                       "area",
                                       "counts",
                                       "background"]

    for method, counter in zip(args.methods, counters):
        profiles = []
        for norm in args.profile_normalizations:
            # build matrix, apply normalization
            profile = counter.getProfile(
                normalize=norm,
                background_region_bins=args.background_region_bins)
            profiles.append(profile)

        for x in range(1, len(profiles)):
            assert profiles[0].shape == profiles[x].shape

        # build a single matrix of all profiles for output
        matrix = numpy.concatenate(profiles)
        matrix.shape = len(profiles), len(profiles[0])
        matrix = matrix.transpose()

        with iotools.open_file(E.get_output_file(counter.name) +
                               ".matrix.tsv.gz", "w") as outfile:
            outfile.write("bin\tregion\tregion_bin\t%s\n" % "\t".join(
                args.profile_normalizations))
            fields = []
            bins = []
            for field, nbins in zip(counter.fields, counter.nbins):
                fields.extend([field] * nbins)
                bins.extend(list(range(nbins)))

            for row, cols in enumerate(zip(fields, bins, matrix)):
                outfile.write("%i\t%s\t" %
                              (row, "\t".join([str(x) for x in cols[:-1]])))
                outfile.write("%s\n" %
                              ("\t".join([str(x) for x in cols[-1]])))

        with iotools.open_file(E.get_output_file(counter.name) +
                               ".lengths.tsv.gz", "w") as outfile:
            counter.writeLengthStats(outfile)

        if args.output_all_profiles:
            counter.closeOutputProfiles()

    if args.plot:

        import matplotlib
        # avoid Tk or any X
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for method, counter in zip(args.methods, counters):

            if method in ("geneprofile",
                          "geneprofilewithintrons",
                          "geneprofileabsolutedistancefromthreeprimeend",
                          "utrprofile",
                          "intervalprofile",
                          "separateexonprofile",
                          "separateexonprofilewithintrons"):

                plt.figure()
                plt.subplots_adjust(wspace=0.05)
                max_scale = max([max(x) for x in counter.aggregate_counts])

                for x, counts in enumerate(counter.aggregate_counts):
                    plt.subplot(6, 1, x + 1)
                    plt.plot(list(range(len(counts))), counts)
                    plt.title(counter.fields[x])
                    plt.ylim(0, max_scale)

                figname = counter.name + ".full"

                fn = E.get_output_file(figname) + "." + args.image_format
                plt.savefig(os.path.expanduser(fn), format=args.image_format, dpi=args.resolution_images)

                plt.figure()

                points = []
                cuts = []
                for x, counts in enumerate(counter.aggregate_counts):
                    points.extend(counts)
                    cuts.append(len(counts))

                plt.plot(list(range(len(points))), points)

                xx, xxx = 0, []
                for x in cuts:
                    xxx.append(xx + x // 2)
                    xx += x
                    plt.axvline(xx,
                                color="r",
                                ls="--")

                plt.xticks(xxx, counter.fields)

                figname = counter.name + ".detail"

                fn = E.get_output_file(figname) + "." + args.image_format
                plt.savefig(os.path.expanduser(fn), format=args.image_format, dpi=args.resolution_images)

            elif method == "tssprofile":

                plt.figure()
                plt.subplot(1, 3, 1)
                plt.plot(list(range(-args.extension_outward,
                                    args.extension_inward)),
                         counter.aggregate_counts[0])
                plt.title(counter.fields[0])
                plt.subplot(1, 3, 2)
                plt.plot(list(range(-args.extension_inward,
                                    args.extension_outward)),
                         counter.aggregate_counts[1])
                plt.title(counter.fields[1])
                plt.subplot(1, 3, 3)
                plt.title("combined")
                plt.plot(list(range(-args.extension_outward,
                                    args.extension_inward)),
                         counter.aggregate_counts[0])
                plt.plot(list(range(-args.extension_inward,
                                    args.extension_outward)),
                         counter.aggregate_counts[1])
                plt.legend(counter.fields[:2])

                fn = E.get_output_file(counter.name) + "." + args.image_format
                plt.savefig(os.path.expanduser(fn), format=args.image_format, dpi=args.resolution_images)

            elif method == "midpointprofile":

                plt.figure()
                plt.plot(numpy.arange(-args.resolution_upstream, 0),
                         counter.aggregate_counts[0])
                plt.plot(numpy.arange(0, args.resolution_downstream),
                         counter.aggregate_counts[1])

                fn = E.get_output_file(counter.name) + "." + args.image_format
                plt.savefig(os.path.expanduser(fn), format=args.image_format, dpi=args.resolution_images)

    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
