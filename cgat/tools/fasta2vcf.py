"""convert fasta to VCF
==========================

Output a file in VCF format with variants according
to a fasta file.


"""

import sys
import random
import cgatcore.experiment as E
import pysam


def main(argv=None):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-s", "--sample-size", dest="sample_size", type="float",
        help="sample size. If less than 0, take a proportion of the chromosome size. "
        "If greater than 0, take a fixed number of variants [%default]")

    parser.set_defaults(
        input_filename_fasta=None,
        sample_size=0.001,
        sample_name="NA12878"
    )

    (options, args) = E.start(parser,
                              argv=argv,
                              add_output_options=True)

    if len(args) > 0:
        options.input_filename_fasta = args[0]

    if options.input_filename_fasta == "-":
        options.input_filename_fasta = options.stdin

    outf = options.stdout
    outf.write("##fileformat=VCFv4.1\n")
    outf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    outf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format(options.sample_name))

    with pysam.FastxFile(options.input_filename_fasta) as inf:
        for record in inf:
            contig = record.name
            sequence = record.sequence
            if options.sample_size < 1.0:
                nsamples = int(float(len(sequence)) * options.sample_size)
            else:
                nsamples = int(options.sample_size)
            E.info("generating {} sampled variants for contig {}".format(nsamples, contig))
            sampled_positions = set()
            missing_nsamples = nsamples
            while len(sampled_positions) < nsamples:
                raw_positions = random.sample(list(range(len(sequence))), nsamples - len(sampled_positions))
                filtered_positions = [x for x in raw_positions if sequence[x] != "N"]
                sampled_positions.update(filtered_positions)
                E.debug("sample update: total={}, raw={}, filtered={}".format(
                        len(sampled_positions),
                        len(raw_positions),
                        len(filtered_positions)))

            sampled_positions = sorted(sampled_positions)

            for position in sampled_positions:
                base = sequence[position]
                outf.write("{}\t{}\t.\t{}\t{}\t.\t.\t.\tGT\t0/0\n".format(
                        contig, position + 1, base, base))

    E.stop()


if __name__ == "__main__":
    sys.exit(main())
