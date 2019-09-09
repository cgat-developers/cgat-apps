'''bam2fasta.py - sequences aligned to a particular region
==========================================================

:Tags: Genomics NGS Sequences BAM FASTA Conversion

Purpose
-------

This script takes a :term:`bam` formatted file and outputs all reads
overlapping a certain region.

Example
-------

For example::

   cat in.bam | cgat bam2fasta

This command converts the :term:`bam` formatted file in.bam into

Type::

   cgat bam2fasta --help

for command line help.

Command line options
--------------------

'''

import os
import re
import sys
import collections
import pysam
import numpy
import pandas
import cgatcore.experiment as E
import cgatcore.iotools as iotools


# adapted from here: https://bitbucket.org/brentp/biostuff/src/282b504ac9020fe1449e23f800b20b5bd7d12061/nwalign/pairwise.py?at=default&fileviewer=file-view-default
# replace with alignlib version
# gop/gep would be good to group gaps
def global_align(seqj, seqi, gap=-1, match=1, mismatch=-1, nmatch=0):
    """
    """
    assert "-" not in seqi
    assert "-" not in seqj
    UP, LEFT, DIAG, NONE = range(4)
    max_j = len(seqj)
    max_i = len(seqi)

    score = numpy.zeros((max_i + 1, max_j + 1), dtype='f')
    pointer = numpy.zeros((max_i + 1, max_j + 1), dtype='i')
    pointer[0, 0] = NONE
    score[0, 0] = 0.0
    pointer[0, 1:] = LEFT
    pointer[1:, 0] = UP

    score[0, 1:] = gap * numpy.arange(max_j)
    score[1:, 0] = gap * numpy.arange(max_i)

    for i in range(1, max_i + 1):
        ci = seqi[i - 1]
        for j in range(1, max_j + 1):
            cj = seqj[j - 1]

            if cj == ci:
                diag_score = score[i - 1, j - 1] + match
            elif cj == 'N' or ci == 'N':
                diag_score = score[i - 1, j - 1] + nmatch
            else:
                diag_score = score[i - 1, j - 1] + mismatch
            up_score = score[i - 1, j] + gap
            left_score = score[i, j - 1] + gap
            if diag_score >= up_score:
                if diag_score >= left_score:
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                else:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
            else:
                if up_score > left_score:
                    score[i, j] = up_score
                    pointer[i, j] = UP
                else:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT

    align_j = ""
    align_i = ""
    while True:
        p = pointer[i, j]
        if p == NONE:
            break
        s = score[i, j]
        if p == DIAG:
            align_j += seqj[j - 1]
            align_i += seqi[i - 1]
            i -= 1
            j -= 1
        elif p == LEFT:
            align_j += seqj[j - 1]
            align_i += "-"
            j -= 1
        elif p == UP:
            align_j += "-"
            align_i += seqi[i - 1]
            i -= 1
        else:
            raise Exception('wtf!')

    return align_j[::-1], align_i[::-1]


def get_consensus(sequences, ignore_gaps=False, min_gap_proportion=0):
    """if *ignore_gaps* is False, only report a gap as the consensus if it
    constitutes at least *min_gap_proportion* at sites.
    """
    positions = list(zip(*sequences))
    counts_df = pandas.DataFrame([collections.Counter(x) for x in positions]).fillna(0)
    if "-" in counts_df.columns:
        if ignore_gaps:
            counts_df.drop(["-"], axis=1, inplace=True)
        else:
            # reset all gap-counts where gap is less than min_gap_proportion
            if min_gap_proportion > 0:
                threshold = min_gap_proportion * float(len(sequences))
                counts_df.loc[counts_df["-"] < threshold, ["-"]] = 0

    consensus = "".join(counts_df.idxmax(axis=1)).upper()
    return consensus


def get_anchor_consensus(sequences):
    # filter sequences, removing anything with non-majority length
    lengths = [len(x) for x in sequences if len(x) > 0]
    median_length = numpy.median(lengths)
    s = [x for x in sequences if len(x) == median_length]
    return get_consensus(s)


def iterate_bed(bed_file, merge_intervals):
    if merge_intervals:
        contig, start, end = None, None, None
        for bed in bed_file.fetch(parser=pysam.asBed()):
            if contig != bed.contig:
                if contig is not None:
                    yield contig, start, end
                contig = bed.contig
                start, end = bed.start, bed.end
            end = bed.end
        yield contig, start, end
    else:
        for bed in bed_file.fetch(parser=pysam.asBed()):
            yield bed.contig, bed.start, bed.end


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-e", "--input-bed-file", dest="input_bed_file", type=str,
        help="input file with intervals. Tab-delimited file of intervals "
        "in bed format to restrict analysis to. ")

    parser.add_argument(
        "-m", "--merge-intervals", dest="merge_intervals", action="store_true",
        help="merge intervals in bed file. Useful if you have a site bed-file "
        )

    parser.add_argument(
        "-f", "--reference-fasta-file", dest="reference_fasta_file",
        help="reference genomic sequence in fasta format. "
        )

    parser.add_argument(
        "-c", "--barcode-fasta-file", dest="barcode_fasta_file",
        help="barcode sequence in fasta format. Variable positions "
        "should be marked by N "
        )

    parser.set_defaults(
        reference_fasta_file=None,
        barcode_fasta_file=None,
        merge_intervals=False,
        input_bed_file=None,
        anchor=5,
    )

    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv, add_output_options=True)

    if args.stdin != sys.stdin:
        bamfile = args.stdin.name
    elif args:
        if len(args) > 1:
            raise ValueError("multiple bam files provided in arguments")
        bamfile = args[0]
    else:
        bamfile = "-"

    if args.barcode_fasta_file:
        with pysam.FastxFile(args.barcode_fasta_file) as inf:
            barcode_sequence = next(inf).sequence
    else:
        barcode_sequence = None

    if not os.path.exists(args.reference_fasta_file):
        raise OSError("reference fasta file {} does not exist".format(
            args.reference_fasta_file))

    if not os.path.exists(args.input_bed_file):
        raise OSError("input bed file {} does not exist".format(
            args.input_bed_file))

    bed_in = pysam.TabixFile(args.input_bed_file)
    pysam_in = pysam.AlignmentFile(bamfile)
    anchor = args.anchor

    for region_idx, vals in enumerate(iterate_bed(bed_in, args.merge_intervals)):

        if region_idx > 0:
            raise NotImplementedError("output for multiple regions not yet implemented")

        contig, region_start, region_end = vals
        upstream_anchors, downstream_anchors = [], []
        counter = E.Counter()

        unaligned_fn = E.get_output_file("unaligned_{}.fasta".format(region_idx))
        with iotools.open_file(unaligned_fn, "w") as outf:
            for read in pysam_in.fetch(contig, region_start, region_end):
                counter.overlapping_reads += 1
                try:
                    pairs = read.get_aligned_pairs(with_seq=True)
                except ValueError:
                    counter.no_md_tag += 1
                    continue

                map_ref2read_pos = dict((x[1], x[0]) for x in pairs if x[0] is not None)
                map_ref2ref_base = dict((x[1], x[2]) for x in pairs if x[0] is not None)

                upstream_anchor = "".join(
                    map_ref2ref_base.get(x, "") for x in range(region_start - anchor, region_start))

                downstream_anchor = "".join(
                    map_ref2ref_base.get(x, "") for x in range(region_end, region_end + anchor))

                # check if at least one anchor is aligned
                upstream_matches = sum([x.isupper() for x in upstream_anchor])
                downstream_matches = sum([x.isupper() for x in downstream_anchor])

                if upstream_matches < anchor and downstream_matches < anchor:
                    counter.no_anchor += 1
                    continue
                seq = read.query_alignment_sequence

                # collect full length anchors
                upstream_anchor_start, upstream_anchor_end = region_start - anchor, region_start
                downstream_anchor_start, downstream_anchor_end = region_end, region_end + anchor

                if upstream_anchor_start in map_ref2read_pos and upstream_anchor_end in map_ref2read_pos:
                    upstream_anchors.append(
                        seq[map_ref2read_pos[upstream_anchor_start]:map_ref2read_pos[upstream_anchor_end]])
                if downstream_anchor_start in map_ref2read_pos and downstream_anchor_end in map_ref2read_pos:
                    downstream_anchors.append(
                        seq[map_ref2read_pos[downstream_anchor_start]:map_ref2read_pos[downstream_anchor_end]])

                # get region to align
                read_start = min((map_ref2read_pos.get(x, len(seq))
                                  for x in range(region_start - anchor, region_start)))
                if read_start == len(seq):
                    read_start = 0
                read_end = max((map_ref2read_pos.get(x, 0) + 1
                                for x in range(region_end, region_end + anchor)))
                if read_end == 1:
                    read_end = len(seq)
                counter.collected_reads += 1
                outf.write(">{}/{}-{}\n{}\n".format(read.query_name,
                                                    read_start, read_end,
                                                    seq[read_start:read_end]))
        counter.downstream_anchors = len(downstream_anchors)
        counter.upstream_anchors = len(upstream_anchors)

        E.info(counter)

        if counter.overlapping_reads == 0:
            E.warn("no sequences overlapping region")
            continue

        if counter.downstream_anchors == 0 or counter.upstream_anchors == 0:
            E.warn("at least one anchor undefined")
            continue

        if counter.collected_reads == 1:
            E.warn("only single sequence, multiple aligment skipped")
            with iotools.open_file(unaligned_fn) as inf:
                stdout = inf.read()
        else:
            # G-INS-i -> global alignment algorithm
            E.info("starting mafft multiple alignment")
            stdout = E.run("mafft --globalpair --maxiterate 100 --quiet --op 2 --ep 0.5 {}".format(unaligned_fn),
                           return_stdout=True)

        aligned_fn = E.get_output_file("aligned_{}.fasta".format(region_idx))
        with iotools.open_file(aligned_fn, "w") as outf:
            outf.write(stdout)

        mali = stdout.splitlines()
        identifiers = [mali[x] for x in range(0, len(mali), 2)]
        sequences = [mali[x].upper() for x in range(1, len(mali), 2)]
        consensus = get_consensus(sequences)

        E.info("after alignment: consensus={}".format(consensus))

        # gap filtering -> remove highly gappy columns
        consensus = get_consensus(sequences, min_gap_proportion=0.9)

        E.info("after anchor trimming: consensus={}".format(consensus))

        take = [idx for idx, x in enumerate(consensus) if x != "-"]
        sequences = ["".join([s[x] for x in take]) for s in sequences]
        consensus = get_consensus(sequences, min_gap_proportion=0.9)

        E.info("after gap filtering: consensus={}".format(consensus))

        # get anchor consensus and chop it off
        consensus = get_consensus(sequences, ignore_gaps=True)
        upstream_anchor = get_anchor_consensus(upstream_anchors)
        downstream_anchor = get_anchor_consensus(downstream_anchors)

        upstream_anchor_start = consensus.find(upstream_anchor)
        downstream_anchor_start = consensus.rfind(downstream_anchor)

        E.info("anchor consensus (no gaps)={}, upstream={}, downstream={}, upstream_idx={}, downstream_idx={}".format(
            consensus, upstream_anchor, downstream_anchor, upstream_anchor_start, downstream_anchor_start))

        if upstream_anchor_start < 0 or downstream_anchor_start < 0:
            E.warn("can't locate anchor, no output produced")
            continue

        upstream_anchor_end = upstream_anchor_start + len(upstream_anchor)
        if upstream_anchor_end >= downstream_anchor_start:
            E.warn("anchor not in correct order, no output produced")
            continue

        sequences = [x[upstream_anchor_end:downstream_anchor_start] for x in sequences]
        consensus = get_consensus(sequences)

        E.info("after anchor trimming: consensus={}".format(consensus))

        truncated_fn = E.get_output_file("aligned_truncated_{}.fasta".format(region_idx))
        with iotools.open_file(truncated_fn, "w") as outf:
            outf.write("\n".join(
                "{}\n{}\n".format(x, y) for x, y in zip(identifiers, sequences)))

        positions = list(zip(*sequences))
        bases = ["A", "C", "G", "T"]
        df = pandas.DataFrame([collections.Counter(x) for x in positions]).fillna(0)
        for missing_base in [x for x in bases if x not in df.columns]:
            df[missing_base] = 0
        df["gapped_depth"] = df.sum(axis=1)
        df["depth"] = df[bases].sum(axis=1)
        df["consensus"] = df[bases].idxmax(axis=1)
        df["consensus_counts"] = df.lookup(df.index, df.consensus)
        df["consensus_support"] = df.consensus_counts / df.depth
        df["offconsensus_counts"] = df.depth - df.consensus_counts
        df.loc[df.consensus_counts == 0, "consensus"] = "N"
        df["region_id"] = region_idx

        # replace "gap" consensus positions with + character
        alignment = global_align(re.sub("-", "+", consensus), barcode_sequence)
        E.info("alignment: consensus {}".format(alignment[0]))
        E.info("alignment: barcode   {}".format(alignment[1]))

        barcode_idx = 0
        deleted_barcode_bases = []
        rows = []
        for c, b in zip(*alignment):
            if c == "-":
                deleted_barcode_bases.append(barcode_idx)
                barcode_idx += 1
            elif b == "N":
                rows.append((barcode_idx, "variable"))
                barcode_idx += 1
            elif b == "-":
                rows.append(("", "insertion"))
            elif b == c:
                rows.append((barcode_idx, "fixed-match"))
                barcode_idx += 1
            else:
                rows.append((barcode_idx, "fixed-mismatch"))
                barcode_idx += 1

        alignment_df = pandas.DataFrame.from_records(rows,
                                                     columns=["barcode_pos", "barcode_class"])

        assert len(alignment_df) == len(df)
        df = pandas.concat([df, alignment_df], axis=1)
        with E.open_output_file("pileup") as outf:
            df.to_csv(outf, sep="\t", index=True, index_label="position")

        observed_barcode_sequence = "".join(df[df.barcode_class == "variable"].consensus)
        headers = df.consensus_support.describe().index
        eval_df = df.loc[df.barcode_class.isin(("variable", "fixed-match", "fixed-mismatch")), ]
        median_consensus_depth = eval_df.consensus_counts.median()
        # zero stuff out if depth is low
        if median_consensus_depth <= 2:
            deleted_barcode_bases = []

        outf = args.stdout
        # modules to recover partial bar-codes
        outf.write("\t".join(map(str,
                                 ["barcode", "ndeleted_barcode_bases", "deleted_barcode_bases"] +
                                 ["support_{}".format(x) for x in headers] +
                                 ["counts_{}".format(x) for x in headers] +
                                 ["offcounts_{}".format(x) for x in headers])) + "\n")

        outf.write("\t".join(map(str, [
            observed_barcode_sequence,
            len(deleted_barcode_bases),
            ",".join(map(str, deleted_barcode_bases))] +
                                 eval_df.consensus_support.describe().tolist() +
                                 eval_df.consensus_counts.describe().tolist() +
                                 eval_df.offconsensus_counts.describe().tolist())) + "\n")

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
