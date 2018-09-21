"""Utility functions for the bam2stats utility."""

import pysam
from pysam.libchtslib cimport *
from pysam.libcbcf cimport VariantFile, VariantRecord, VariantRecordSample
from pysam.libcfaidx cimport FastaFile
from pysam.libctabix cimport TabixFile
from libc.string cimport strchr
from libc.stdint cimport int8_t
from libc.stdio cimport puts, printf
from cpython cimport array as c_array
import array
import base64
import collections
import collections
import hashlib
import itertools
import math
import numpy
import pandas
import re
import string
import sys
import cgatcore.experiment as E


# deprecated: avoid R dependencies in cgat-apps
# import rpy2.robjects
# from rpy2.robjects import r as R
# from rpy2.robjects.packages import importr
# from rpy2.robjects import pandas2ri

cimport numpy


def generate_from_bed(vcf_file, bed_file):
    for bed in bed_file.fetch(parser=pysam.asBed()):
        for v in vcf_file.fetch(bed.contig, bed.start, bed.end):
            yield v

def generate_from_vcf(vcf_file):
    for v in vcf_file:
        yield v

def generate_from_region(vcf_file, region):
    for v in vcf_file.fetch(region=region):
        yield v


ACGT = str.maketrans("ACGT", "TGCA")


cdef class Counter(object):

    cdef int nsamples
    cdef object samples
    cdef bint only_variant_positions

    def __init__(self, samples, only_variant_positions=False):
        self.nsamples = len(samples)
        self.samples = samples
        self.only_variant_positions = only_variant_positions

    cdef process_record(self, VariantRecord record, bint is_snp):
        raise NotImplementedError("base class must implement process_record")

    
cdef class CounterKinship(Counter):

    cdef int32_t * data_genotype_ptr

    cdef numpy.ndarray n_Aa
    cdef numpy.ndarray n_AAaa
    cdef numpy.ndarray n_AaAa
    cdef numpy.ndarray genotype_values
    
    def __init__(self, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs)

        self.data_genotype_ptr = NULL
        self.n_Aa = numpy.zeros(self.nsamples, dtype=numpy.int64)
        # if memory becomes an issue, use triangular matrices
        self.n_AAaa = numpy.zeros((self.nsamples, self.nsamples), dtype=numpy.int64)
        self.n_AaAa = numpy.zeros((self.nsamples, self.nsamples), dtype=numpy.int64)
        self.genotype_values = numpy.zeros(self.nsamples, dtype=numpy.int8)

    cdef process_record(self, VariantRecord record, bint is_snp):

        if not is_snp:
            return

        cdef int mdat = 0
        cdef int32_t * ptr = NULL
        cdef int allele
        cdef bint is_het_i, is_het_j
        cdef int nret
        cdef int _i, _j

        cdef int8_t [:] genotype_values_view = self.genotype_values

        nret = bcf_get_genotypes(
            record.header.ptr,
            record.ptr,
            &self.data_genotype_ptr,
            &mdat)

        if nret == 0:
            raise ValueError("no genotypes")

        nret /= self.nsamples
        ptr = <int32_t*>self.data_genotype_ptr
        for _i from 0 <= _i < self.nsamples:
            # sum up alleles. Missing alleles or more than 2 alleles 
            # will set the value to a large negative value
            genotype_values_view[_i] = 0
            for _j from 0 <= _j < nret:
                if bcf_gt_is_missing(ptr[_j]):
                    genotype_values_view[_i] = -10
                    continue
                if ptr[_j]==bcf_int32_vector_end:
                    break
                allele = bcf_gt_allele(ptr[_j])
                if allele > 1:
                    genotype_values_view[_i] = -10
                genotype_values_view[_i] += allele
            ptr += nret

        for _i from 0 <= _i < self.nsamples:
            if genotype_values_view[_i] < 0:
                continue
            is_het_i = genotype_values_view[_i] == 1
            if is_het_i:
                self.n_Aa[_i] += 1
            for _j from 0 <= _j < _i:
                if genotype_values_view[_j] < 0:
                    continue
                is_het_j = genotype_values_view[_j] == 1
                if is_het_i and is_het_j:
                    self.n_AaAa[_i][_j] += 1
                elif (not is_het_i and
                      not is_het_j and
                      genotype_values_view[_i] != genotype_values_view[_j]):
                    self.n_AAaa[_i][_j] += 1

    def output(self):

        cdef int _i, _j
        with E.open_output_file("kinship") as outf:
            outf.write("sample_i\tsample_j\twithin_kinship\t"
                       "between_kinship\tn_het_i\tn_het_j\tn_homhom\tn_hethet\n")
            for _i from 0 <= _i < self.nsamples:
                for _j from 0 <= _j < _i:
                    try:
                        between_kinship = (self.n_AaAa[_i][_j] - 2.0 * self.n_AAaa[_i][_j]) \
                            / (2.0 * self.n_Aa[_i]) \
                            + 0.5 - 0.25 * (self.n_Aa[_i] + self.n_Aa[_j]) / self.n_Aa[_i]
                        between_kinship_f = "{:6.4f}".format(between_kinship)
                    except ZeroDivisionError:
                        between_kinship_f = ""

                    try:
                        within_kinship = (self.n_AaAa[_i][_j] - 2.0 * self.n_AAaa[_i][_j]) \
                            / (self.n_Aa[_i] + self.n_Aa[_j])
                        within_kinship_f = "{:6.4f}".format(within_kinship)
                    except ZeroDivisionError:
                        within_kinship_f = ""

                    outf.write("\t".join(map(
                        str,
                        (self.samples[_i],
                         self.samples[_j],
                         within_kinship_f,
                         between_kinship_f,
                         self.n_Aa[_i],
                         self.n_Aa[_j],
                         self.n_AAaa[_i][_j],
                         self.n_AaAa[_i][_j]))) + "\n")

    def __dealloc__(self):
        if self.data_genotype_ptr is not NULL:
            free(self.data_genotype_ptr)


cdef class CounterFormatDistributions(Counter):
    
    cdef numpy.ndarray counts
    cdef numpy.ndarray unset_samples
    cdef numpy.ndarray unset_sites
    cdef numpy.ndarray values

    cdef object codes
    cdef int ncodes
    cdef int nbins

    def __init__(self, codes, nbins, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs)

        self.nbins = nbins
        # use byte values to facilitate internal processing
        self.codes = [x.encode("ascii") for x in codes]
        self.ncodes = len(codes)

        self.counts = numpy.zeros(
            (self.ncodes,
             self.nsamples,
             int(nbins) + 1), dtype=numpy.int64)

        self.unset_samples = numpy.zeros(
            (self.ncodes,
             self.nsamples),
            dtype=numpy.int64)

        self.unset_sites = numpy.zeros(
            (self.ncodes,
             self.nsamples + 1),
            dtype=numpy.int64)
        
        self.values = numpy.zeros(
            self.nsamples,
            dtype=numpy.int32)

    cdef process_record(self, VariantRecord record, bint is_snp):

        cdef int64_t [:, :, :] counts_view = self.counts
        cdef int64_t [:, :] unset_samples_view = self.unset_samples
        cdef int64_t [:, :] unset_sites_view = self.unset_sites
        cdef int32_t [:] values_view = self.values
        cdef int32_t * data_ptr = &values_view[0]
        cdef int code_idx, sample_idx, ival
        cdef char * code
        cdef int ndest

        unset = numpy.zeros(self.ncodes, dtype=numpy.int64)
        cdef int64_t [:] unset_view = unset

        # quickly iterate over all values for a particular format
        for code_idx, code in enumerate(self.codes):
            bcf_get_format_int32(
                record.header.ptr,
                record.ptr,
                code,
                &data_ptr,
                &ndest)

            for sample_idx from 0 <= sample_idx < self.nsamples:
                ival = data_ptr[sample_idx]
                if ival == bcf_int32_missing:
                    unset_view[code_idx] += 1
                    unset_samples_view[code_idx][sample_idx] += 1
                elif ival < 0:
                    raise ValueError("received negative value {} for format field {}".format(
                        ival, code))
                else:
                    if ival > self.nbins:
                        ival = self.nbins
                    counts_view[code_idx][sample_idx][ival] += 1

        for code_idx, ival in enumerate(unset_view):
            unset_sites_view[code_idx][ival] += 1

    def output(self):

        str_codes = [x.decode("ascii") for x in self.codes]

        with E.open_output_file("format_per_sample") as outf:
            outf.write("FORMAT\tbin\t{}\n".format(
                "\t".join(self.samples)))

            for code_idx, code in enumerate(str_codes):
                bins = numpy.arange(0, self.nbins + 1, 1)
                bins.shape = (self.nbins + 1, 1)
                mm = self.counts[code_idx].transpose()
                mm = numpy.hstack([bins, mm])
                assert mm.shape == (self.nbins + 1, self.nsamples + 1)

                # some metrics are discretized, so remove any all-0 rows            
                row_sums = mm.sum(axis=1)
                comp = mm[:, 0]
                mm = mm[comp != row_sums]
                for row in mm:
                    outf.write("{}\t{}\n".format(code, "\t".join(map(str, row))))

        with E.open_output_file("format_unset_samples") as outf:
            outf.write("FORMAT\t{}\n".format(
                "\t".join(self.samples)))
            for code_idx, code in enumerate(str_codes):
                outf.write("{}\t{}\n".format(
                    code,
                    "\t".join(map(str, self.unset_samples[code_idx]))))
        bins = numpy.arange(0, self.nsamples + 1)
        bins.shape = (self.nsamples + 1, 1)

        # output as binary for numpy savetxt function
        with E.open_output_file("format_unset_sites", mode="wb") as outf:
            header = "bin\t{}\n".format(
                "\t".join(str_codes))
            outf.write(header.encode("ascii"))
            mm = numpy.hstack([bins, self.unset_sites.transpose()])
            numpy.savetxt(outf,
                          mm,
                          delimiter="\t",
                          fmt="%i")


cdef class CounterMutationalSignature(Counter):

    cdef object signatures
    cdef object profile

    cdef FastaFile fasta_in

    def __init__(self, fasta_in, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs)
        self.fasta_in = fasta_in

        # Mutation profile
        #
        # C -> A   G    T   C>A = C->A, G->T
        # G    T   C -> A
        #
        # C -> G   C    G   C>G = C->G, G->C
        # G    C   G -> C
        #
        # C -> T   G    A   C>T = C->T, G->A
        # G    A   C -> T
        #
        # T -> A   T    A   T>A = T->A, A->T
        # A    T   A -> T
        #
        # T -> C   A    G   T>C = T->C, A->G
        # A    G   T -> C

        # T -> G   A    C   T>G = T->G, A->C
        # A    C   T -> G

        self.profile = {
            "CA": ("C>A", 0),
            "GT": ("C>A", 1),
            "CG": ("C>G", 0),
            "GC": ("C>G", 1),
            "CT": ("C>T", 0),
            "GA": ("C>T", 1),
            "TA": ("T>A", 0),
            "AT": ("T>A", 1),
            "TC": ("T>C", 0),
            "AG": ("T>C", 1),
            "TG": ("T>G", 0),
            "AC": ("T>G", 1)
        }
        assert len(self.profile.keys()) == 12

        self.signatures = collections.defaultdict(
            lambda: collections.defaultdict(
                lambda: collections.defaultdict(int)))

    cdef process_record(self, VariantRecord record, bint is_snp):

        # note: pos in VariantFile is 1-based
        if not is_snp:
            return

        alt = record.alts[0]
        context = self.fasta_in.fetch(record.chrom,
                                      record.pos - 2,
                                      record.pos + 1).upper()

        assert context[1] == record.ref, \
            "reference sequence mismatch? expected {} at {}:{}, got {}".format(
                context[1], record.chrom, record.pos, record.ref)

        # skip any non ACGT bases
        if not(context[0] in "ACGT" and
               context[1] in "ACGT" and
               context[2] in "ACGT"):
            return

        profile_class, reverse = self.profile[record.ref + alt]
        if reverse:
            context = context.translate(ACGT)[::-1]

        c = "{}.{}".format(context[0], context[2])

        all_alleles = collections.defaultdict(int)
        first_variants = None
        is_unique = False
        for s in self.samples:
            ai = list(record.samples[s].allele_indices)
            variants = [x for x in ai if x and x > 0]

            if len(variants) > 0:
                self.signatures[s][profile_class][c] += 1
                for v in variants:
                    all_alleles[v] += 1

            if first_variants is None:
                first_variants = tuple(sorted(variants))
            else:
                if tuple(sorted(variants)) != first_variants:
                    is_unique = True

        if is_unique:
            self.signatures["unique"][profile_class][c] += 1

    def output(self):

        with E.open_output_file("mutation_profile") as outf:
            outf.write("sample\tsignature\tcontext\tcount\t"
                       "percent_sample\tpercent_context\n")
            for s, dd in sorted(self.signatures.items()):
                total_class = collections.defaultdict(int)
                for cls, d in dd.items():
                    t = 0
                    for alt, count in d.items():
                        t += count
                    total_class[cls] = float(t)
                total_sample = sum(total_class.values())

                for cls, d in sorted(dd.items()):
                    for context, count in sorted(d.items()):
                        outf.write("\t".join(map(
                            str,
                            (s, cls, context, count,
                             100.0 * count / total_sample,
                             100.0 * count / total_class[cls])
                        )) + "\n")



cdef class CounterMutationalSignatureProfile(CounterMutationalSignature):
    """compute signutares from mutational profile.

    This module requires the deconstructSigs R library. To install:
    source("https://bioconductor.org/biocLite.R")
    biocLite()
    biocLite(c("BSgenome", "BSgenome.Hsapiens.UCSC.hg19", "GenomeInfoDb"))
    install.packages("deconstructSigs")
    """

    cdef object signatures_database
    cdef object counts_method

    def __init__(self, *args, **kwargs):
        raise NotImplementedError("R dependencies are not resolved")
    #     CounterMutationalSignature.__init__(self, *args, **kwargs)
    #     deconstructSigs = importr('deconstructSigs')
    #     R('''library(deconstructSigs)''')

    #     self.signatures_database = "signatures.cosmic"
    #     self.counts_method = "genome"
        
    # def output(self):

    #     colnames = list(R('''colnames(randomly.generated.tumors)'''))

    #     # colnames are of format:
    #     # [1] "A[C>A]A" "A[C>A]C" "A[C>A]G" "A[C>A]T" "C[C>A]A" "C[C>A]C" "C[C>A]G"
    #     # [8] "C[C>A]T" "G[C>A]A" "G[C>A]C" "G[C>A]G" "G[C>A]T" "T[C>A]A" "T[C>A]C"

    #     rows = []
    #     samples = []
    #     for sample, dd in self.signatures.items():
    #         values = {}
    #         for cls, d in dd.items():
    #             for context, count in d.items():
    #                 # context = "G.G", cls = C>T
    #                 values["{}[{}]{}".format(
    #                     context[0], cls, context[2])] = count

    #         extra = set(values.keys()).difference(colnames)
    #         missing = set(colnames).difference(values.keys())
    #         if len(extra) > 0:
    #             raise AssertionError(
    #                 "unexpected additional columns in sample {}: {}".format(
    #                     sample, sorted(extra)))
    #         if len(missing) > 0:
    #             E.warn("missing counts for sample {}: {} will be set to 0".format(
    #                 sample, sorted(missing)))

    #         samples.append(sample)
    #         rows.append([values.get(x, 0) for x in colnames])

    #     df = pandas.DataFrame.from_records(rows,
    #                                        columns=colnames,
    #                                        index=samples)

    #     # convert to R dataframe. The conversion translates special charactes
    #     # such as []< to ., so set column names explicitely.
    #     rdf = pandas2ri.py2ri(df)
    #     rdf.colnames = rpy2.robjects.StrVector(colnames)
    #     R.assign("rdf", rdf)
    #     results = []
    #     for sample in df.index:
    #         R('''result = whichSignatures(rdf,
    #         sample.id='{sample}',
    #         signatures.ref={signatures},
    #         contexts.needed=TRUE,
    #         signature.cutoff=0,
    #         tri.counts.method="{counts_method}")'''.format(
    #             sample=sample,
    #             signatures=self.signatures_database,
    #             counts_method=self.counts_method))
    #         results.append(R('''result$weights'''))

    #     results = pandas.concat(results, axis=0)
    #     results.columns = [re.sub("[.]", "", x).lower() for x in results.columns]
    #     results.index.name = "sample"
    #     results = results.reset_index()

    #     with E.open_output_file("mutation_profile_signatures") as outf:
    #         results.to_csv(outf, sep="\t", index=False)


cdef class CounterGCContext(Counter):

    cdef FastaFile fasta_in
    cdef int nbins
    cdef int window_size
    cdef numpy.ndarray counts

    def __init__(self, fasta_in, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs)
        self.fasta_in = fasta_in
        self.nbins = 100
        self.counts = numpy.zeros((self.nsamples, self.nbins + 1))
        self.window_size = 50

    cdef process_record(self, VariantRecord record, bint is_snp):

        if not is_snp:
            return
        context = self.fasta_in.fetch(record.chrom,
                                      max(0, record.pos - self.window_size),
                                      record.pos + self.window_size).upper()
        counts = collections.Counter(context)
        gc = counts["G"] + counts["C"]
        at = counts["A"] + counts["T"]
        gc_content = int(math.floor(100.0 * gc / (gc + at)))
        c = "{}.{}".format(context[0], context[2])
        for idx, s in enumerate(self.samples):
            try:
                ai = list(record.samples[s].allele_indices)
            except TypeError:
                # no GT field
                continue

            if self.only_variant_positions:
                if len([x for x in ai if x > 0]):
                    continue

            self.counts[idx][gc_content] += 1

    def output(self):
        with E.open_output_file("gc_context") as outf:
            outf.write("percent_gc\t{}\n".format("\t".join(self.samples) ))
            bins = numpy.arange(0, self.nbins + 1, 1)
            bins.shape = (self.nbins + 1, 1)
            mm = self.counts.transpose()
            mm = numpy.hstack([bins, mm])
            assert mm.shape == (self.nbins + 1, self.nsamples + 1)
            for row in mm:
                outf.write("{}\n".format("\t".join(map(str, row))))


cdef class CounterGCDepthProfile(Counter):

    cdef FastaFile fasta_in
    cdef int nbins_gc
    cdef int nbins_dp
    cdef int window_size
    cdef numpy.ndarray counts

    def __init__(self, fasta_in, gc_window_size=50, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs)
        self.fasta_in = fasta_in
        self.nbins_gc = 100
        self.nbins_dp = 1000
        self.counts = numpy.zeros(
            (self.nsamples, self.nbins_dp + 1, self.nbins_gc + 1))
        self.window_size = gc_window_size

    cdef process_record(self, VariantRecord record, bint is_snp):

        if not is_snp:
            return

        context = self.fasta_in.fetch(record.chrom,
                                      max(0, record.pos - self.window_size),
                                      record.pos + self.window_size).upper()

        counts = collections.Counter(context)
        gc = counts["G"] + counts["C"]
        at = counts["A"] + counts["T"]
        gc_content = int(math.floor(100.0 * gc / (gc + at)))

        c = "{}.{}".format(context[0], context[2])

        for idx, s in enumerate(self.samples):
            try:
                ai = list(record.samples[s].allele_indices)
            except TypeError:
                # no GT field
                continue

            if self.only_variant_positions:
                if len([x for x in ai if x > 0]):
                    continue
                    
            depth = record.samples[s]["DP"]
            if depth != None:
                depth = min(depth, self.nbins_dp)
                self.counts[idx][depth][gc_content] += 1

    def output(self):
        with E.open_output_file("gc_dp_prof") as outf:
            gc_bins = numpy.arange(0, self.nbins_gc + 1, 1)

            dp_bins = numpy.arange(0, self.nbins_dp + 1, 1)
            dp_bins.shape = (self.nbins_dp + 1, 1)

            outf.write("sample\tdp\tbin\t{}\n".format("\t".join(map(str, gc_bins))))

            for sample_idx, sample in enumerate(self.samples):
                mm = self.counts[sample_idx]
                mm = numpy.hstack([dp_bins, mm])
                assert mm.shape == (self.nbins_dp + 1, self.nbins_gc + 2)
                # some metrics are discretized, so remove any all-0 rows
                row_sums = mm.sum(axis=1)
                comp = mm[:, 0]
                mm = mm[comp != row_sums]
                for row in mm:
                    outf.write("{}\t{}\n".format(sample, "\t".join(map(str, row))))

        with E.open_output_file("gc_dp_profile_stats") as outf:
            gc_bins = numpy.arange(0, self.nbins_gc + 1, 1)
            dp_bins = numpy.arange(0, self.nbins_dp + 1, 1)
            dp_bins.shape = (self.nbins_dp + 1, 1)
            outf.write ("sample\tgc_bin\tmean\n")
            for sample_idx, sample in enumerate(self.samples):
                mm = self.counts[sample_idx]
                if mm.sum() == 0:
                    continue

                sums = (dp_bins * mm).sum(axis=0)
                counts = mm.sum(axis=0)
                means = numpy.nan_to_num((sums / counts))

                for gc_bin, mean in zip(gc_bins, means):
                    outf.write("{}\t{}\t{:.2f}\n".format(sample, gc_bin, mean))


def vcf2stats_count(VariantFile vcf_in,
                    FastaFile fasta_in,
                    TabixFile bed_in,
                    options):

    if options.region is not None:
        vcf_records = generate_from_region(vcf_in, options.region)
    elif bed_in is not None:
        vcf_records = generate_from_bed(vcf_in, bed_in)
    else:
        vcf_records = generate_from_vcf(vcf_in)

    samples = list(vcf_in.header.samples)

    counters = []

    for method in options.methods:
        if method == "mutational-signature":
            counters.append(CounterMutationalSignature(
                fasta_in=fasta_in,
                samples=samples))
        elif method == "mutational-signature-profile":
            counters.append(CounterMutationalSignatureProfile(
                fasta_in=fasta_in,
                samples=samples))
        elif method == "kinship":
            counters.append(CounterKinship(
                samples=samples))
        elif method == "format-distribution":
            counters.append(CounterFormatDistributions(
                nbins=options.format_distribution_nbins,
                codes=options.format_distributions,
                samples=samples))
        elif method == "gc-context":
            counters.append(CounterGCContext(
                fasta_in=fasta_in,
                samples=samples))
        elif method == "gc-depth-profile":
            counters.append(CounterGCDepthProfile(
                fasta_in=fasta_in,
                samples=samples,
                gc_window_size=options.gc_window_size,
                only_variant_positions=options.only_variant_positions))

    cdef bint is_snp

    cdef int report_step = options.report_step
    cdef int record_idx
    cdef VariantRecord record

    cdef Counter counter

    for record_idx, record in enumerate(vcf_records):
        if (record_idx % options.report_step) == 0:
            E.debug("iteration {}: {}:{}".format(record_idx, record.chrom, record.pos))

        if record.filter.values()[0].name != "PASS":
            continue

        # skip <NON_REF>
        alts = record.alts
        if alts[0][0] == "<":
            continue

        ref = record.ref
        is_snp = (len(ref) == 1 and
                  len(alts) == 1 and
                  len(alts[0]) == 1)

        for counter in counters:
            counter.process_record(record, is_snp)

    for counter in counters:
        counter.output()
