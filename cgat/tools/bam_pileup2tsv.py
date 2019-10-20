"""Per position statistics from a BAM file
==========================================

Compute per-position statistics.

Methods
-------

read-variant:
    Output for each position the bases of the aligned reads.

read-list:
    Output for each position the reads that are aligned to this
    position.

depth-vcf:
    Base coverage at position (excluding deletions). Note that
    the output of this method is in :term:`VCF` format.

coverage-vcf:
    Read coverage at position (including reads that have a gap
    at site). Note that the output of this method is in :term:`VCF`
    format.

barcode:

"""

import sys
import os
import json
import re
import pandas
import pysam
import cgatcore.experiment as E


def generate_from_bed(bam_file, bed_file, **kwargs):
    for bed in bed_file.fetch(parser=pysam.asBed()):
        for v in bam_file.pileup(bed.contig, bed.start, bed.end,
                                 **kwargs,
                                 truncate=True):
            yield v


def generate_from_bam(bam_file, **kwargs):
    for v in bam_file.pileup(**kwargs):
        yield v


def generate_from_region(bam_file, region, **kwargs):
    for v in bam_file.pileup(region=region, **kwargs):
        yield v


def main(argv=None):

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser()

    parser.add_argument(
        "-i", "--input-fastq-file", dest="input_fastq_file", type=str,
        help="input fastq file. "
        )

    parser.add_argument(
        "-m", "--method", dest="method", type=str,
        choices=("read-variant", "depth-vcf", "read-list", "coverage-vcf", "barcode"),
        help="method to apply ")

    parser.add_argument(
        "-e", "--input-bed", dest="input_bed_file", type=str,
        help="input file with intervals. Tab-delimited file of intervals "
        "in bed format to restrict analysis to. ")

    parser.add_argument(
        "-r", "--region-string", dest="region_string", type=str,
        help="region string. Only apply method in specified region. "
        )

    parser.add_argument(
        "-f", "--reference-fasta-file", dest="reference_fasta_file",
        help="reference genomic sequence in fasta format. "
        )

    parser.add_argument(
        "--min-base-quality", dest="min_base_quality", type=int,
        help="minimum base quality for barcode analysis. "
        )

    parser.add_argument(
        "-s", "--stepper", dest="stepper", type=str,
        choices=("nofilter", "samtools", "all"))

    parser.set_defaults(
        method="read-variant",
        reference_fasta_file=None,
        input_bed_file=None,
        regex_sample_name="([^/]+).bam",
        stepper="nofilter",
        min_base_quality=13,
        region_string=None)

    # add common options (-h/--help, ...) and parse command line
    (args, unknowns) = E.start(parser, argv=argv,
                               add_output_options=True,
                               unknowns=True)

    pysam_in = pysam.AlignmentFile(unknowns[0], "rb")

    if args.input_bed_file:
        if not os.path.exists(args.input_bed_file):
            raise OSError("input bed file {} does not exist".format(
                args.input_bed_file))
        bed_in = pysam.TabixFile(args.input_bed_file)
    else:
        bed_in = None

    if args.region_string is not None:
        itr = generate_from_region(pysam_in, args.region,
                                   stepper=args.stepper,
                                   min_base_quality=args.min_base_quality)
    elif bed_in is not None:
        itr = generate_from_bed(pysam_in, bed_in,
                                stepper=args.stepper,
                                min_base_quality=args.min_base_quality)
    else:
        itr = generate_from_bam(pysam_in,
                                stepper=args.stepper,
                                min_base_quality=args.min_base_quality)

    reference_fasta = pysam.FastaFile(args.reference_fasta_file)

    outf = args.stdout
    counter = E.Counter()

    if args.method == "read-variant":
        outf.write("chromosome\tposition\tref\ttypes\n")

        for pileupcolumn in itr:
            counter.positions_pileup += 1
            reference_base = reference_fasta.fetch(
                pileupcolumn.reference_name,
                pileupcolumn.reference_pos,
                pileupcolumn.reference_pos + 1)
            matches = []
            bases = set()
            for read in pileupcolumn.pileups:
                qpos = read.query_position
                if qpos is not None:
                    base = read.alignment.query_sequence[qpos]
                else:
                    base = "-"

                matches.append((base,
                                read.alignment.query_name))
                bases.add(base)

            bases = list(bases)
            if len(bases) == 1:
                counter.position_noninformative += 1
                if bases[0] == reference_base:
                    counter.position_reference += 1
                continue

            counter.position_informative += 1

            d = {}
            for base in bases:
                d[base] = ",".join([x[1] for x in matches if x[0] == base])

            outf.write("{}\t{}\t{}\t{}\n".format(
                pileupcolumn.reference_name,
                pileupcolumn.reference_pos,
                reference_base,
                json.dumps(d)))

    elif args.method in ("depth-vcf", "coverage-vcf"):
        if args.regex_sample_name:
            sample_name = re.search(args.regex_sample_name, unknowns[0]).groups()[0]
        else:
            sample_name = "unknown"

        outf.write("##fileformat=VCFv4.1\n")
        outf.write("##FORMAT=<ID=GT,Number=1,Type=String,"
                   "Description=\"Genotype\">\n")
        outf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,"
                   "Description=\"Genotype\">\n")
        outf.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\t"
            "FILTER\tINFO\tFORMAT\t{}\n".format(sample_name))

        is_depth = args.method == "depth-vcf"

        for idx, pileupcolumn in enumerate(itr):

            if idx % 1000 == 0:
                E.info("processed {} positions".format(idx))

            reference_base = reference_fasta.fetch(
                pileupcolumn.reference_name,
                pileupcolumn.reference_pos,
                pileupcolumn.reference_pos + 1).upper()

            if reference_base == 'A':
                alt_base = 'C'
            else:
                alt_base = 'A'

            if is_depth:
                n = sum([1 for x in pileupcolumn.pileups if not (x.is_del or x.is_refskip)])
            else:
                n = pileupcolumn.n

            outf.write("{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT:DP\t0/1:{}\n".format(
                pileupcolumn.reference_name,
                pileupcolumn.reference_pos,
                reference_base,
                alt_base,
                n))

    elif args.method == "read-list":
        outf.write("chromosome\tposition\treference_base\tbase\tquality\tquery_name\n")

        for pileupcolumn in itr:
            reference_base = reference_fasta.fetch(pileupcolumn.reference_name,
                                                   pileupcolumn.reference_pos,
                                                   pileupcolumn.reference_pos + 1)
            matches = []
            for read in pileupcolumn.pileups:
                qpos = read.query_position
                if qpos is not None:
                    base = read.alignment.query_sequence[qpos]
                    quality = read.alignment.query_qualities[qpos]
                else:
                    base = "-"
                    quality = ""

                outf.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    pileupcolumn.reference_name,
                    pileupcolumn.reference_pos,
                    reference_base,
                    base,
                    quality,
                    read.alignment.query_name))

    elif args.method == "barcode":

        rows = []
        for c in itr:
            rows.append((c.reference_pos,
                         c.n,
                         "".join(c.get_query_sequences()),
                         pysam.qualities_to_qualitystring(c.get_query_qualities())))
        df = pandas.DataFrame.from_records(
            rows,
            columns=["pos", "gapped_depth", "bases", "qualities"])

        df["depth"] = df.bases.str.len()
        bases = ["A", "C", "G", "T"]
        for b in bases:
            df[b] = df.bases.str.upper().str.count(b)
        df["consensus"] = df[bases].idxmax(axis=1)
        df["consensus_counts"] = df.lookup(df.index, df.consensus)
        df["consensus_support"] = df.consensus_counts / df.depth
        df["offconsensus_counts"] = df.depth - df.consensus_counts
        df.loc[df.consensus_counts == 0, "consensus"] = "N"

        df.to_csv(outf, sep="\t", index=False)

    E.info(counter)
    # write footer and output benchmark information.
    E.stop()
