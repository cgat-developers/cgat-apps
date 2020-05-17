"""compare two bed files.
"""

import collections
import pysam
import sys
import cgatcore.experiment as E
import cgatcore.iotools as iotools
import quicksect


def get_size_bin(size, size_bins):
    bin = 0
    nsize_bins = len(size_bins)
    while bin < nsize_bins and size > size_bins[bin]:
        bin += 1
    return bin


def main(argv=sys.argv):

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-b", "--reference-bed-file", dest="reference_bed_file", type=str,
        help="reference bed file "
        )

    parser.add_argument(
        "-m", "--method", dest="method", type=str,
        choices=("lvc-comparison", ),
        help="methods to apply ")

    parser.set_defaults(
        method="lvc-comparison",
        reference_fasta_file=None,
        input_bed_file=None,
        size_bins=(1000, 10000, 100000),
        output_sets=True,
        region_string=None)

    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv, add_output_options=True)

    reference_set = collections.defaultdict(
        quicksect.IntervalTree)

    E.info("reading reference bed file from {}".format(args.reference_bed_file))
    with iotools.open_file(args.reference_bed_file) as inf:
        for record in pysam.tabix_iterator(inf, pysam.asBed()):
            mm = reference_set[record.contig]
            mm.add(record.start,
                   record.end)
    E.info("read reference intervals on {} contigs: {}".format(
            len(list(reference_set.keys())), ",".join(list(reference_set.keys()))))

    if args.output_sets:
        output_tp = E.open_output_file("tp")
        output_fp = E.open_output_file("fp")
        output_fn = E.open_output_file("fn")
    else:
        output_tp = None
        output_fp = None
        output_fn = None

    if args.method == "lvc-comparison":
        c = E.Counter()

        found = set()
        counts = {}
        names = set()
        nsize_bins = len(args.size_bins)
        for bin in range(len(args.size_bins) + 1):
            counts[bin] = dict([(x, collections.defaultdict(int)) for x in
                                ("tp", "fn", "fp", "test", "truth")])

        for record in pysam.tabix_iterator(args.stdin, pysam.asBed()):
            if record.contig not in reference_set:
                c.ignored_no_contig += 1
                continue

            c.test += 1
            matches = reference_set[record.contig].search(record.start, record.end)
            size = record.end - record.start
            bin = get_size_bin(size, args.size_bins)

            if len(matches) == 0:
                c.fp += 1
                status = "fp"
                if output_fp:
                    output_fp.write(str(record) + "\n")
            elif len(matches) >= 1:
                c.tp += 1
                status = "tp"
                if output_tp:
                    output_tp.write(str(record) + "\n")
                # todo: overlap criteria

                # record found
                for match in matches:
                    found.add((record.contig, match.start, match.end))

            name = record.name.split(",")[0]
            names.add(name)
            counts[bin]["test"][name] += 1
            counts[bin][status][name] += 1

        outf = args.stdout

        with iotools.open_file(args.reference_bed_file) as inf:
            for record in pysam.tabix_iterator(inf, pysam.asBed()):
                c.truth += 1
                bin = get_size_bin(record.end - record.start, args.size_bins)
                counts[bin]["truth"]["all"] += 1

                key = (record.contig, record.start, record.end)
                if key not in found:
                    c.fn += 1
                    counts[bin]["fn"]["all"] += 1

        outf.write("\t".join(("category",
                              "size",
                              "test",
                              "tp",
                              "fp",
                              "truth",
                              "fn")) + "\n")

        for name in sorted(names):
            for bin in range(len(args.size_bins) + 1):
                if bin == len(args.size_bins):
                    size_bin = ">={}".format(args.size_bins[-1])
                else:
                    size_bin = "<{}".format(args.size_bins[bin])
                outf.write("\t".join(map(str, (
                                name,
                                size_bin,
                                counts[bin]["test"][name],
                                counts[bin]["tp"][name],
                                counts[bin]["fp"][name],
                                counts[bin]["truth"]["all"],
                                counts[bin]["fn"]["all"],
                                ))) + "\n")

    E.info(str(c))
    E.stop()
