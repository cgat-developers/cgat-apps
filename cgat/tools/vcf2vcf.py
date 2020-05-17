'''vcf2vcf.py - manipulate vcf files
=================================

Purpose
-------

Manipulate vcf-formatted files.


Usage
-----

Type::

   python vcf2vcf.py --help

for command line usage.

Methods
-------

This script provides the following methods:

re-order
   reorder sample columns in vcf formatted file according to a given sort order

Documentation
-------------

This is a tool for manipulating vcf-formatted files.  The following
options are available:

+-----------+-------------------------+
+-----------+-------------------------+

lift-over
^^^^^^^^^

Command line options
--------------------

'''

import os
import sys
import random
import collections
import numpy
import re
import pysam
import quicksect
import cgatcore.experiment as E
import cgatcore.iotools as iotools


def read_liftover_chain(infile):

    E.debug("started reading mapping information")

    map_id2chromosome = ["", ]
    map_chromosome2id = {}
    n = 0

    Chain = collections.namedtuple(
        "Chain",
        ["score",
         "target_name", "target_size", "target_strand",
         "target_start", "target_end",
         "query_name", "query_size", "query_strand",
         "query_start", "query_end", "chainid"])

    def blocks(infile):

        keep = False
        for line in infile:
            if line.startswith("chain"):
                chain_data = Chain._make(line[:-1].split(" ")[1:])

                if chain_data.target_strand == "-":
                    raise NotImplementedError("target strand is negative")
                alignment_data = []
            elif line.strip() == "":
                yield chain_data, alignment_data
            else:
                alignment_data.append(list(map(int, line.split(("\t")))))

    map_chromosomes = collections.defaultdict(
        quicksect.IntervalTree)
    map_contig2length = collections.defaultdict(int)

    for chain_data, alignment_data in blocks(infile):

        map_contig2length[chain_data.query_name] = int(chain_data.query_size)

        # target maps to query
        # coordinates are zero-based, half-open When
        # the strand value is "-", position coordinates are listed in
        # terms of the reverse-complemented sequence
        x = int(chain_data.target_start)
        y = int(chain_data.query_start)
        # revert coordinates for negative strands (it seems that
        # the mapping file uses reverse coordinates, while liftover
        # output doesn't)
        invert = chain_data.query_strand == "-"
        mm = map_chromosomes[chain_data.target_name]

        for d in alignment_data:
            if len(d) == 3:
                size, increment_x, increment_y = d
            else:
                size, increment_x, increment_y = d[0], 0, 0

            mm.add(x, x + size,
                   (chain_data.query_name,
                    y,
                    y + size,
                    invert))

            x += increment_x + size
            y += increment_y + size

            if y < 0:
                raise ValueError(
                    "illegal mapping in chain {}".format(chain_data))

    return map_chromosomes, map_contig2length


def main(argv=None):

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument(
        "--input-filename-fasta", dest="input_filename_fasta", type=str,
        help="filename with reference sequence in fasta format ")

    parser.add_argument(
        "--input-filename-bam", dest="input_filename_bam", type=str,
        help="filename with aligned reads ")

    parser.add_argument(
        "--method", dest="methods", type=str, action="append",
        choices=["add-strelka-genotype",
                 "lift-over"],
        help="methods to apply ")

    parser.add_argument(
        "--input-filename-chain", dest="input_filename_chain", type=str,
        help="filename with alignment chain for lift-over ")

    parser.add_argument(
        "--normal-sample-regex", dest="normal_sample_regex", type=str,
        help="regular expression to apply to header to identify normal "
        "sample id ")

    parser.add_argument(
        "--output-filename-unmapped", dest="output_filename_unmapped", type=str,
        help="filename with variants that could not be lifted over ")

    parser.set_defaults(
        input_filename_fasta=None,
        input_filename_bam=None,
        input_filename_vcf="-",
        sample_size=0.001,
        region_size=20,
        methods=[],
        normal_sample_regex=None,
        input_filename_chain=None,
        output_filename_unmapped=None,
    )

    (args, unknown) = E.start(parser,
                              argv=argv,
                              add_output_options=True,
                              unknowns=True)

    if len(unknown) > 0:
        args.input_filename_vcf = unknown[0]

    vcf_in = pysam.VariantFile(args.input_filename_vcf)

    if "lift-over" in args.methods:
        if args.input_filename_chain is None:
            raise ValueError("--method=lift-over requires --input-filename-chain")
        if not os.path.exists(args.input_filename_chain):
            raise OSError("file {} with chain data does not exist".format(
                args.input_filename_chain))
        E.info("reading chain from {}".format(args.input_filename_chain))
        with iotools.open_file(args.input_filename_chain) as inf:
            map_chain, map_contig2length = read_liftover_chain(inf)

    if args.input_filename_fasta:
        fasta = pysam.FastaFile(args.input_filename_fasta)
    else:
        fasta = None

    if args.input_filename_bam:
        bam = pysam.AlignmentFile(args.input_filename_bam)
    else:
        bam = None

    outf = args.stdout

    c = E.Counter()

    if "add-strelka-genotype" in args.methods:
        map_nt2gt = {"ref": "0/0",
                     "het": "0/1",
                     "hom": "1/1",
                     "conflict": "."}

        map_tumour2gt = {"ref": "0/0",
                         "het": "0/1",
                         "hom": "1/1"}

        header = str(vcf_in.header).splitlines()

        header.insert(
            len(header) - 1,
            '##FORMAT=<ID=GT,Number=1,Type=String,Description='
            '"Genotypes of reference and alternative alleles, '
            'added by cgatcore vcf2vcf.">')

        header = "\n".join(header)
        if args.normal_sample_regex:
            normal_sample = re.search(" -bam-file \S+/([^/]+)_S\d+.bam", header).groups()[0]
        else:
            normal_sample = "NORMAL"

        is_first = True

        for record in vcf_in:
            c.input += 1

            if "GT" in record.format:
                if is_first:
                    outf.write(header + "\n")
                    is_first = False
                outf.write(str(record))
                c.has_gt += 1
                continue

            gt_normal = map_nt2gt[record.info["NT"]]
            gt_tumour = record.info["SGT"]
            norm, tumour = gt_tumour.split("->")
            if gt_tumour[0] in "ACGT":
                alts = record.alts
                if alts is None:
                    c.no_alt += 1
                    continue

                if len(record.alts) > 1:
                    c.multi_allelic += 1
                    continue

                _map_tumour2gt = {
                    record.alts[0]: "1",
                    record.ref: "0"}
                try:
                    gt_tumour = "/".join(
                        sorted([_map_tumour2gt[x] for x in tumour]))
                except KeyError:
                    gt_tumour = "."
                    c.ambigous_genotype += 1
            else:
                gt_tumour = map_tumour2gt[tumour]

            fields = str(record)[:-1].split("\t")
            # FORMAT
            fields[8] = ":".join(("GT", fields[8]))
            # SAMPLES
            # makes a few assumptions, fix!
            header_insert_normal = False
            if len(fields) == 11:
                fields[9] = ":".join((gt_normal, fields[9]))
                fields[10] = ":".join((gt_tumour, fields[10]))
            elif len(fields) == 10:
                header_insert_normal = True
                values = fields[9].split(":")
                fields.append(":".join((gt_tumour, fields[9])))
                fields[9] = ":".join([gt_normal] + ["."] * len(values))
            else:
                raise NotImplementedError()

            if is_first:
                if not header_insert_normal:
                    outf.write(header + "\n")
                else:
                    header = re.sub(r"\tFORMAT\t",
                                    "\tFORMAT\t%s\t" % normal_sample, header)
                    outf.write(header + "\n")
                is_first = False
            outf.write("\t".join(fields) + "\n")
            c.output += 1

    elif "lift-over" in args.methods:
        header = str(vcf_in.header).splitlines()

        if fasta:
            # validate contig size
            expected_lengths = dict(list(zip(fasta.references, fasta.lengths)))
        else:
            expected_lengths = map_contig2length

        # update contig names and sizes in VCF header
        header = [x for x in header if not x.startswith("##contig")]
        header[-1:-1] = ["##contig=<ID={},length={}>".format(
            contig, length) for contig, length in sorted(expected_lengths.items())]

        header.insert(
            len(header) - 1,
            '##liftover=<CHAIN={},REFERENCE={}>'.format(
                args.input_filename_chain,
                args.input_filename_fasta))
        outf.write("\n".join(header) + "\n")

        unmapped_contigs = set()
        unknown_contigs = set()

        trans_genotypes = str.maketrans("01", "10")

        if fasta:
            # validate contig size
            expected_lengths = dict(list(zip(fasta.references, fasta.lengths)))
            for contig, length in list(map_contig2length.items()):
                if contig in expected_lengths:
                    if length != expected_lengths[contig]:
                        raise ValueError(
                            "contig lengths mismatch. For contig {} chain files "
                            "says {}, but fasta files says {}".format(
                                contig, length, expected_lengths[contig]))
            E.info("contig sizes in chain file and fasta files correspond.")

        if args.output_filename_unmapped:
            outfile_unmapped = iotools.open_file(args.output_filename_unmapped, "w")
            outfile_unmapped.write("\n".join(header) + "\n")
        else:
            outfile_unmapped = None

        for record in vcf_in:
            c.input += 1

            try:
                mm = map_chain[record.contig]
            except KeyError:
                c.skipped_unmapped_contig += 1
                unmapped_contigs.add(record.contig)
                if outfile_unmapped:
                    outfile_unmapped.write("skipped_unmapped_contig\t{}".format(str(record)))
                continue

            try:
                m = mm.search(record.start, record.stop)
            except AttributeError:
                c.skipped_mapping_error += 1
                if outfile_unmapped:
                    outfile_unmapped.write("skipped_mapping_error\t{}".format(str(record)))
                continue

            if len(m) == 0:
                c.skipped_unmapped_position += 1
                if outfile_unmapped:
                    outfile_unmapped.write("skipped_unmapped_position\t{}".format(str(record)))
                continue
            elif len(m) > 1:
                c.skipped_multimapping_position += 1
                if outfile_unmapped:
                    outfile_unmapped.write("skipped_multimapping_position\t{}".format(str(record)))
                continue

            m = m[0]
            y_contig, y_start, y_end, y_invert = m.data

            if y_invert:
                y_pos = y_end - (record.start - m.start)
            else:
                y_pos = (record.start - m.start) + y_start

            if fasta:
                try:
                    ref_base = fasta.fetch(y_contig, y_pos, y_pos + len(record.ref)).upper()
                except KeyError:
                    c.skipped_unknown_contig += 1
                    unknown_contigs.add(y_contig)
                    ref_base = None
                    continue

            swap_alleles = False
            if ref_base:
                error = False
                if ref_base == record.ref:
                    c.matches += 1
                else:
                    if len(record.alts) == 1:
                        alt_base = record.alts[0]
                        if ref_base == alt_base:
                            swap_alleles = True
                            c.allele_swap_variant += 1
                        else:
                            c.error_mismatch_variant += 1
                            error = "mismatch"
                    else:
                        error = "multi-mismatch"
                        c.error_multi_mismatch_variant += 1

                if error:
                    if outfile_unmapped:
                        outfile_unmapped.write("{}\t{}".format(error, str(record)))
                    c.skipped_error_variant += 1
                    continue

            fields = str(record)[:-1].split("\t")
            fields[0] = y_contig
            fields[1] = str(y_pos)

            if swap_alleles:
                fields[4] = alt_base
                fields[5] = ref_base
                # update genotype fields
                keep = False
                for idx in range(9, len(fields)):
                    gt, rest = fields[idx].split(":", 1)
                    keep = keep or "0" in gt
                    fields[idx] = ":".join((gt.translate(trans_genotypes), rest))

                # remove reference only calls
                if not keep:
                    if outfile_unmapped:
                        outfile_unmapped.write("reference_call\t{}".format(str(record)))
                    c.skipped_allele_swap_reference += 1
                continue

            c.output += 1
            outf.write("\t".join(fields) + "\n")

        c.unmapped_contigs = len(unmapped_contigs)
        c.unknown_contigs = len(unknown_contigs)

        E.info(c.asTable())
        if unknown_contigs:
            E.info("unknown contigs: {}".format(",".join(sorted(unknown_contigs))))
        if unmapped_contigs:
            E.info("unmapped contigs: {}".format(",".join(sorted(unmapped_contigs))))

    E.stop()


if __name__ == "__main__":
    sys.exit(main())
