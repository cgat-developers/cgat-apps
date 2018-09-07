"""Compute statistics from a VCF file.

Implemeted methods
==================

mutation-profile
----------------

Compute a mutation profile for each sample in the VCF file. For each
central base in a tri-nucleotide, the following equivalence classes
are used::

    R    A   R    A   Code     Forward strand substitutions
    ---------------
    C -> A   G    T   C>A    = C->A, G->T
    G    T   C -> A

    C -> G   C    G   C>G    = C->G, G->C
    G    C   G -> C

    C -> T   G    A   C>T    = C->T, G->A
    G    A   C -> T

    T -> A   T    A   T>A    = T->A, A->T
    A    T   A -> T

    T -> C   A    G   T>C    = T->C, A->G
    A    G   T -> C

    T -> G   A    C   T>G    = T->G, A->C
    A    C   T -> G

kinship
-------

Compute kinship coefficient for all pairs of samples in VCF.
The kinship coefficient is estimated using the robust estimator
by Manichaikul et al (2010).

The robust estimator for between family relationship is (eq. 11):

phi_ij = N_{Aa, Aa} - 2 N_{AA,aa} / (2 N_{Aa}(i))
         + 1/2
         - 1/4 * N_{Aa}(i) + N_{Aa}(j) / N_{Aa}(i)

The robust estimator for within family relationship is (eq. 9):

phi_ij = N_{Aa,Aa} - 2 N_{AA,aa} / (N_{Aa}(i) + N_{Aa}(j))

with:

N_{Aa,Aa}: # of variants where both individuals are heterozygous
N_{AA,aa}: # of variants where both individuals are homozygous different
N_{Aa}(i): # of heterozygous variants in individual i

format-distribution
-------------------

Compute distribution of one or more metrics in the format field. This method
outputs several tables:

format_per_sample

   Histogram over the FORMAT field. The first two columns in the table
   are `FORMAT`, the FORMAT field specifier, and `bin`. These columns
   are followed by each sample as a column.

format_unset_samples

   Table showing the number of unset FORMAT fields per sample. The
   first column is FORMAT followed by columns for each sample.

format_unset_sites

   Table showing the distribution of sites for each FORMAT field that
   have no annotation for a particular field. This table has the
   column `bin` followed by one column for each FORMAT field.

gc-depth-profile
----------------

gc-context
----------



"""

import os
import sys
import pysam

import cgatcore.experiment as E

from cgat.VCFTools import vcf2stats_count


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-i", "--input-vcf", dest="input_vcf_file", type="string",
        help="input vcf file")

    parser.add_option(
        "-f", "--input-fasta", dest="input_fasta_file", type="string",
        help="input fasta file. faidx indexed reference sequence file to "
        "determine INDEL context [%default]")

    parser.add_option(
        "-e", "--input-bed", dest="input_bed_file", type="string",
        help="input file with intervals. Tab-delimited file of intervals "
        "in bed format to restrict analysis to. [%default]")

    parser.add_option(
        "-r", "--region", dest="region", type="string",
        help="Region string to restrict analysis to. Takes precedence "
        "over --input-bed. [%default]")

    parser.add_option(
        "-m", "--method", dest="methods", action="append", type="choice",
        choices=("mutational-signature",
                 "mutational-signature-profile",
                 "kinship",
                 "format-distribution",
                 "gc-context",
                 "gc-depth-profile"),
        help="methods to apply [%default]")

    parser.add_option(
        "--format-distribution", dest="format_distributions", action="append",
        type="string",
        help="format to compute histograms on. Option can specified multiple times. "
        "At the moment, only integer metrics are supported [%default]")

    parser.add_option(
        "--format-distribution-nbins", dest="format_distributions_nbins", type="int",
        help="number of bins to use for histograms [%default]")

    parser.add_option(
        "--only-variant-positions", dest="only_variant_positions",
        action="store_true",
        help="only use variant positions [%default]")

    parser.add_option(
        "--gc-window-size", dest="gc_window_size", type="int",
        help="(half) window size to use for G+C computation. A size "
        "of 50 means that 50 bases on either side of the variant are "
        "used to compute the G+C content [%default]")

    parser.set_defaults(
        methods=[],
        input_vcf_file=None,
        input_bed_file=None,
        region=None,
        input_fasta_file=None,
        format_distributions=[],
        format_distribution_nbins=1000,
        gc_window_size=50,
        report_step=1000000,
    )

    (options, args) = E.start(parser, argv, add_output_options=True)

    if len(args) == 1:
        options.input_vcf_file = args[0]

    if options.input_vcf_file is None:
        raise ValueError("please supply a VCF file")

    if options.input_fasta_file is None:
        raise ValueError("please supply a FASTA file")

    if "format-distribution" in options.methods and not options.format_distributions:
        raise ValueError("please supply at least one FORMAT field (DP, GQ) "
                         "when --method=format-distribution has been selected")

    if not os.path.exists(options.input_vcf_file):
        raise OSError("input vcf file {} does not exist".format(
            options.input_vcf_file))

    if not os.path.exists(options.input_vcf_file + ".tbi") and not \
       os.path.exists(options.input_vcf_file + ".csi"):
        raise OSError("input vcf file {} needs to be indexed".format(
            options.input_vcf_file))

    if not os.path.exists(options.input_fasta_file):
        raise OSError("input fasta file {} does not exist".format(
            options.input_fasta_file))

    if not os.path.exists(options.input_fasta_file + ".fai"):
        raise OSError("input fasta file {} needs to be indexed".format(
            options.input_fasta_file))

    # update paths to absolute
    options.input_fasta_file = os.path.abspath(options.input_fasta_file)
    options.input_vcf_file = os.path.abspath(options.input_vcf_file)

    # catch issue with empty variant files
    try:
        vcf_in = pysam.VariantFile(options.input_vcf_file)
    except (OSError, ValueError):
        E.warn("could not open variant file - likely to be empty")
        E.stop()
        return 0

    fasta_in = pysam.FastaFile(options.input_fasta_file)

    if options.input_bed_file:
        if not os.path.exists(options.input_bed_file):
            raise OSError("input bed file {} does not exist".format(
                options.input_bed_file))
        bed_in = pysam.TabixFile(options.input_bed_file)
    else:
        bed_in = None

    vcf2stats_count(
        vcf_in, fasta_in, bed_in, options)

    E.stop()


if __name__ == "__main__":
    sys.exit(main())
