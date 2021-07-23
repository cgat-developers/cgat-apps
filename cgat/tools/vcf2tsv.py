"""convert fasta to VCF
==========================

Output a file in VCF format with variants according
to a fasta file.

Note: currently only takes first sample in VCF and assumes
single-sample VCF

"""

import sys
import random
import collections
import cgatcore.experiment as E
import numpy
import pysam


class Counter(object):
    pass


class CounterIndelType(Counter):

    header = ["indel_type",
              "indel_length",
              "indel_delta",
              "indel_class"]

    def count(self, record):
        self.ignore = False
        self.indel_length = ""
        self.indel_class = ""
        self.indel_type = ""
        self.indel_delta = ""
        ref = record.ref
        alt = record.alts[0]
        if len(record.alts) > 1:
            return
        lref = len(ref)
        lalt = len(alt)
        if lref == lalt:
            return
        elif lref > lalt:
            self.indel_type = "deletion"
            self.indel_delta = ref[1:]
        else:
            self.indel_type = "insertion"
            self.indel_delta = alt[1:]

        self.indel_length = len(self.indel_delta)
        counts = collections.Counter(self.indel_delta)
        if len(counts) == 1 and self.indel_length > 2:
            self.indel_class = "mononucleotide-instability"
        elif len(counts) == 2 and self.indel_length > 4:
            self.indel_class = "dinucleotide-instability"
        else:
            self.indel_class = "other"

    def __str__(self):
        return "\t".join(map(
            str,
            (self.indel_type, self.indel_length,
             self.indel_delta, self.indel_class)))


class CounterContext(Counter):

    header = ["class", "left_class", "right_class", "context"]

    def __init__(self, fasta):

        if fasta is None:
            raise ValueError("CounterContext requires an indexed fasta file")
        self.fasta = fasta

        self.region = 20

    def count(self, record):

        pos = record.pos - 1
        s = self.fasta.fetch(record.chrom,
                             pos - self.region,
                             pos + self.region + 1)

        def classify(region, threshold_low_complexity_region=5, threshold_homopolymer=10):
            c = collections.Counter(region)

            # monomers define homopolymers
            if len(list(c.values())) == 0:
                return "empty"

            if max(c.values()) > threshold_homopolymer:
                return "hompolymer"

            # dimers define low-complexity regions
            c = collections.Counter([region[x:x+2] for x in range(len(region)-1)])
            if len(list(c.values())) == 0:
                return "empty"

            if max(c.values()) > threshold_low_complexity_region:
                return "low-complexity"

            return "-"

        left_region = s[:self.region]
        right_region = s[self.region + 1:]
        self.cls_left = classify(left_region)
        self.cls_right = classify(right_region)

        self.context = left_region.lower() + s[self.region] + right_region.lower()

        full_cls = (self.cls_left, self.cls_right)

        if self.cls_left == self.cls_right:
            self.cls = self.cls_left
            return

        if "-" in full_cls:
            self.cls = "border"
        else:
            self.cls = "repetetive"

    def __str__(self):
        return "\t".join((self.cls, self.cls_left, self.cls_right, self.context))


class CounterBAM(Counter):

    def __init__(self, bam):
        if bam is None:
            raise ValueError("CounterBAM requires an indexed BAM file")

        self.bam = bam


class CounterBAMIndels(CounterBAM):

    header = ["mean_depth",
              "mean_freq_deletions",
              "mean_freq_insertions",
              "mean_freq_indels",
              "positions_with_deletions",
              "positions_with_insertions",
              "positions with_indels",
              "indel_status_string"]

    def __init__(self, *args, **kwargs):
        CounterBAM.__init__(self, *args, **kwargs)

        self.region = 5

    def count(self, record):

        pos = record.pos - 1
        varset = None
        depths = []
        deletions = []
        insertions = []
        for column in self.bam.pileup(record.chrom,
                                      pos - self.region,
                                      pos + self.region + 1,
                                      truncate=True,
                                      stepper="all"):
            depths.append(len(column.pileups))
            deletions.append(sum([x.is_del == 1 for x in column.pileups]))
            insertions.append(sum([x.indel > 0 for x in column.pileups]))

        self.depths = numpy.array(depths)
        self.deletions = numpy.array(deletions, dtype=numpy.float)
        self.insertions = numpy.array(insertions, dtype=numpy.float)

    def __str__(self):

        indels = self.deletions + self.insertions
        indel_status = numpy.floor(numpy.true_divide(indels, self.depths) * 10.0)
        indel_status_string = "".join(map(str, list(map(int, indel_status))))
        return "\t".join(map(str, (
                "{:.2f}".format(numpy.mean(self.depths)),
                "{:.4f}".format(numpy.mean(self.deletions / self.depths)),
                "{:.4f}".format(numpy.mean(self.insertions / self.depths)),
                "{:.4f}".format(numpy.mean(indels / self.depths)),
                sum(self.deletions > 0),
                sum(self.insertions > 0),
                sum(indels > 0),
                indel_status_string,
                )))


class CounterBAMAllelicDepth(CounterBAM):

    header = ["nbases", "ngaps",
              "nref", "nalt", "nother",
              "freq_gaps",
              "freq_ref", "freq_alt", "freq_other",
              "genotype"]

    def __init__(self, *args, **kwargs):
        CounterBAM.__init__(self, *args, **kwargs)

        self.region = 5

    def count(self, record):

        pos = record.pos - 1
        varset = None
        depths = []
        deletions = []
        insertions = []
        for column in self.bam.pileup(record.chrom,
                                      pos,
                                      pos+1,
                                      truncate=True,
                                      stepper="all"):

            bases = [x.alignment.query_sequence[x.query_position]
                     for x in column.pileups if x.query_position is not None]
            self.ngaps = len([x for x in column.pileups if x.query_position is None])
            self.nbases = len(bases)
            self.nreference = len([x for x in bases if x == record.ref])
            if len(record.alts) > 1:
                self.is_multiallelic = True
            alt = record.alts[0]
            self.nalt = len([x for x in bases if x == alt])
            self.nother = self.nbases - self.nreference - self.nalt
            gt = sorted(set(record.samples[0]["GT"]))
            if len(gt) == 1:
                if gt[0] == 1:
                    self.genotype = "hom-alt"
                elif gt[0] == 0:
                    self.genotype = "hom-ref"
                else:
                    self.genotype = "hom-other"
            else:
                self.genotype = "het"

    def __str__(self):
        if self.nbases > 0:
            f_ref = "{:.4f}".format(float(self.nreference) / self.nbases)
            f_alt = "{:.4f}".format(float(self.nalt) / self.nbases)
            f_other = "{:.4f}".format(float(self.nother) / self.nbases)
        else:
            f_ref = f_alt = f_other = "na"

        if self.ngaps + self.nbases > 0:
            f_gaps = "{:.4f}".format(float(self.ngaps) / (self.nbases + self.ngaps))
        else:
            f_gaps = "na"

        return "\t".join(map(str, (
                    self.nbases,
                    self.ngaps,
                    self.nreference,
                    self.nalt,
                    self.nother,
                    f_gaps,
                    f_ref, f_alt, f_other,
                    self.genotype
                    )))


def main(argv=None):

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument(
        "-s", "--sample-size", dest="sample_size", type=float,
        help="sample size. If less than 0, take a proportion of the chromosome size. "
        "If greater than 0, take a fixed number of variants ")

    parser.add_argument(
        "--input-filename-fasta", dest="input_filename_fasta", type=str,
        help="filename with reference sequence in fasta format ")

    parser.add_argument(
        "--input-filename-bam", dest="input_filename_bam", type=str,
        help="filename with aligned reads ")

    parser.add_argument(
        "--no-vcf-columns", dest="no_vcf_columns", action="store_true",
        help="do not output vcf columns")

    parser.add_argument(
        "--counter", dest="counters", type=str, action="append",
        choices=["context", "bam-indels", "bam-allelic-depth", "indel-type"],
        help="counters to apply ")

    parser.set_defaults(
        input_filename_fasta=None,
        input_filename_bam=None,
        input_filename_vcf=None,
        sample_size=0.001,
        sample_name="NA12878",
        region_size=20,
        threshold_homopolymer=12,
        threshold_repeat=5,
        no_vcf_columns=False,
        counters=[],
    )

    (args, unknown) = E.start(parser,
                              argv=argv,
                              add_output_options=True,
                              unknowns=True)

    if len(unknown) > 0:
        args.input_filename_vcf = unknown[0]

    vcf_in = pysam.VariantFile(args.input_filename_vcf)

    counters = []

    if args.input_filename_fasta:
        fasta = pysam.FastaFile(args.input_filename_fasta)
    else:
        fasta = None

    if args.input_filename_bam:
        bam = pysam.AlignmentFile(args.input_filename_bam)
    else:
        bam = None

    for counter in args.counters:
        if counter == "context":
            counters.append(CounterContext(fasta))
        elif counter == "bam-indels":
            counters.append(CounterBAMIndels(bam))
        elif counter == "bam-allelic-depth":
            counters.append(CounterBAMAllelicDepth(bam))
        elif counter == "indel-type":
            counters.append(CounterIndelType())

    outf = args.stdout
    if not args.no_vcf_columns:
        header = str(vcf_in.header).strip().split("\n")[-1].strip()[1:].split("\t")

    else:
        header = ["chrom", "pos"]

    outf.write("\t".join(header))

    for counter in counters:
        outf.write("\t" + "\t".join(counter.header))

    outf.write("\n")
    for record in vcf_in:

        for counter in counters:
            counter.count(record)

        if not args.no_vcf_columns:
            outf.write("{}\t".format(
                str(record).strip()))
        else:
            outf.write("{}\t{}\t".format(
                record.chrom, record.pos))

        outf.write("\t".join(map(str, counters)) + "\n")

    E.stop()


if __name__ == "__main__":
    sys.exit(main())
