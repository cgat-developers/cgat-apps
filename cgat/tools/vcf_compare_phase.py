"""Compare phasing information in two VCF files.

"""

import sys
import os
import pysam
import collections
import string
import cgatcore.experiment as E


def main(argv=sys.argv):

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument(
        "-i", "--input-vcf", dest="input_vcf_file", type=str,
        help="input vcf file")

    parser.add_argument(
        "-t", "--truth-vcf", dest="truth_vcf_file", type=str,
        help="truth vcf file")

    parser.add_argument(
        "-f", "--input-fasta", dest="input_fasta_file", type=str,
        help="input fasta file. faidx indexed reference sequence file to "
        "determine INDEL context ")

    parser.add_argument(
        "-e", "--input-bed", dest="input_bed_file", type=str,
        help="input file with intervals. Tab-delimited file of intervals "
        "in bed format to restrict analysis to. ")

    parser.add_argument(
        "-m", "--method", dest="methods", action="append", type=str,
        choices=("mutational-signature", "kinship"),
        help="methods to apply ")

    parser.set_defaults(
        methods=[],
        input_vcf_file=None,
        input_bed_file=None,
        input_fasta_file=None,
        truth_vcf_file=None,
    )

    (args, unknown) = E.start(parser,
                              argv,
                              add_output_options=True,
                              unknowns=True)

    if len(unknown) == 1:
        args.input_vcf_file = unknown[0]

    if args.input_vcf_file is None:
        raise ValueError("please supply a VCF file")

    if args.truth_vcf_file is None:
        raise ValueError("please supply a VCF file with truth data")

    if args.input_fasta_file is None:
        raise ValueError("please supply a fasta file with the reference genome")

    if not os.path.exists(args.input_vcf_file):
        raise OSError("input vcf file {} does not exist".format(
            args.input_vcf_file))

    if not os.path.exists(args.input_vcf_file + ".tbi"):
        raise OSError("input vcf file {} needs to be indexed".format(
            args.input_vcf_file))

    if not os.path.exists(args.truth_vcf_file):
        raise OSError("truth vcf file {} does not exist".format(
            args.truth_vcf_file))

    if not os.path.exists(args.truth_vcf_file + ".tbi"):
        raise OSError("truth vcf file {} needs to be indexed".format(
            args.truth_vcf_file))

    if not os.path.exists(args.input_fasta_file):
        raise OSError("input fasta file {} does not exist".format(
            args.input_fasta_file))

    if not os.path.exists(args.input_fasta_file + ".fai"):
        raise OSError("input fasta file {} needs to be indexed".format(
            args.input_fasta_file))

    # update paths to absolute
    args.input_fasta_file = os.path.abspath(args.input_fasta_file)
    args.input_vcf_file = os.path.abspath(args.input_vcf_file)
    args.truth_vcf_file = os.path.abspath(args.truth_vcf_file)

    test_vcf = pysam.VariantFile(args.input_vcf_file)
    truth_vcf = pysam.VariantFile(args.truth_vcf_file)
    contigs = test_vcf.header.contigs
    truth_contigs = set(truth_vcf.header.contigs)

    test_vcf_samples = set(test_vcf.header.samples)
    truth_vcf_samples = set(truth_vcf.header.samples)

    common_samples = test_vcf_samples.intersection(truth_vcf_samples)
    if len(common_samples) == 0:
        raise ValueError("no common samples in test/truth VCFs")

    def pair_iterator(test_vcf, truth_vcf, contig):
        counter = E.Counter()
        test_iter = test_vcf.fetch(contig)
        truth_iter = truth_vcf.fetch(contig)

        test_record = next(test_iter)
        truth_record = next(truth_iter)
        try:
            while 1:
                if test_record.pos < truth_record.pos:
                    test_record = next(test_iter)
                    continue

                elif test_record.pos > truth_record.pos:
                    truth_record = next(truth_iter)
                    continue

                elif len(test_record.alts) > 1:
                    counter.skip_test_truth += 1
                    test_record = next(test_iter)
                    continue

                elif len(truth_record.alts) > 1:
                    counter.skip_multiallelic_truth += 1
                    truth_record = next(truth_iter)
                    continue

                elif test_record.alts != truth_record.alts:
                    counter.skip_genotype_difference += 1
                    test_record = next(test_iter)
                    truth_record = next(truth_iter)
                    continue

                if test_record.ref != truth_record.ref:
                    # todo: deal with indels
                    raise ValueError(
                        "mismatching reference bases at position "
                        "{}:{}".format(test_record.chrom, test_record.pos))

                yield test_record, truth_record
                test_record = next(test_iter)
                truth_record = next(truth_iter)

        except StopIteration:
            pass

        E.debug(str(counter))

    counters_per_contig = {}

    for contig in contigs:
        counter_contig = collections.defaultdict(E.Counter)
        counters_per_contig[contig] = counter_contig

        E.info("processing contig {}".format(contig))

        if contig not in truth_contigs:
            E.warn(
                "skipping contig {} as it is not in truth data".format(contig))
            continue

        switch = False
        last_is_unphased = True

        for test_record, truth_record in pair_iterator(test_vcf, truth_vcf, contig):

            for sample in common_samples:
                counter = counter_contig[sample]

                truth_phased = truth_record.samples[sample].phased
                test_phased = test_record.samples[sample].phased
                truth_genotype = truth_record.samples[sample]["GT"]
                test_genotype = test_record.samples[sample]["GT"]
                truth_alleles = set(truth_genotype)
                test_alleles = set(test_genotype)

                ignore = False
                if not truth_phased:
                    counter.truth_unphased += 1
                    ignore = True
                if not test_phased:
                    counter.test_unphased += 1
                    ignore = True
                    last_is_unphased = True
                else:
                    last_is_unphased = False

                if len(test_alleles) == 1:
                    counter.test_homozygous += 1
                    ignore = True
                else:
                    if not test_phased:
                        counter.test_unphased_hets += 1

                if len(truth_alleles) == 1:
                    counter.truth_homozygous += 1
                    ignore = True

                if ignore:
                    counter.ignore += 1
                    continue

                E.debug("comparing: {}:{} {} -> {}: {} {}".format(
                    test_record.chrom, test_record.pos,
                    test_record.ref, test_record.alts,
                    test_genotype,
                    truth_genotype))

                if switch:
                    truth_genotype = truth_genotype[::-1]

                counter.test_phased_hets += 1

                if truth_genotype != test_genotype:
                    if not last_is_unphased:
                        E.debug("SWITCH: {}".format(switch))
                        counter.switch += 1
                    switch = not switch

    outf = args.stdout
    outf.write("\t".join(("contig",
                          "sample",
                          "switch_error_percent",
                          "false_negative_rate",
                          "switches",
                          "test_phased_hets",
                          "test_unphased_hets",
                          "test_unphased",
                          "truth_unphased",
                          "test_homozygous",
                          "truth_homozygous")) + "\n")

    for contig, contig_dict in list(counters_per_contig.items()):
        for sample, c in list(contig_dict.items()):
            outf.write("\t".join(
                map(str, (
                    contig,
                    sample,
                    "{:6.4f}".format(100.0 * c.switch / (c.test_phased_hets + 1)),
                    "{:6.4f}".format(100.0 * c.test_unphased_hets /
                                     (c.test_phased_hets + c.test_unphased_hets)),
                    c.switch,
                    c.test_phased_hets,
                    c.test_unphased_hets,
                    c.test_unphased,
                    c.truth_unphased,
                    c.test_homozygous,
                    c.truth_homozygous))) + "\n")

    E.stop()
    # TODO: aggregate stats per samples
