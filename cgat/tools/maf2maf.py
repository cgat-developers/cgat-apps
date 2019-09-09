"""

This script currently only processes pairwise MAF alignments.

"""

import sys
import re
import collections
import cgatcore.experiment as E
import cgatcore.iotools as iotools
from cgat.Genomics import parse_region_string


def iterate_maf_blocks(infile):
    block = []
    for line in infile:
        if line.startswith("a"):
            if block:
                yield block
            block = []
        block.append(line)
    yield block

# The following is derived from bx.python (MIT licence)


def format_tabular(rows, align=None):
    if len(rows) == 0:
        return ""
    lengths = [len(col) for col in rows[0]]
    for row in rows[1:]:
        for i in range(0, len(row)):
            lengths[i] = max(lengths[i], len(row[i]))
    rval = ""
    for row in rows:
        for i in range(0, len(row)):
            if align and align[i] == "l":
                rval += row[i].ljust(lengths[i])
            else:
                rval += row[i].rjust(lengths[i])
            rval += " "
        rval += "\n"
    return rval


def parse_block(block):
    RECORD = collections.namedtuple(
        "RECORD",
        ("key", "src", "start", "size", "strand", "srcsize", "text"))

    key, src, start, size, strand, srcsize, text = re.split("\s+", block[1].strip())
    query = RECORD(key, src, int(start), int(size), strand, int(srcsize), text)
    key, src, start, size, strand, srcsize, text = re.split("\s+", block[2].strip())
    target = RECORD(key, src, int(start), int(size), strand, int(srcsize), text)
    try:
        key, text = re.split("\s+", block[3].strip())
    except ValueError:
        qual = None
    return block[0], query, target, qual


def main(argv=sys.argv):

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-f", "--input-filter-tsv",
        dest="input_filter_tsv", type=str,
        help="list with identifiers to remove. "
        )

    parser.add_argument(
        "--set-prefix", dest="set_prefix", type=str,
        help="set sequence prefix ")

    parser.add_argument(
        "--min-length", dest="min_length", type=int,
        help="minimum alignment length ")

    parser.add_argument(
        "--method", dest="methods", action="append",
        choices=("shift-region", ),
        help="methods to apply ")

    parser.set_defaults(
        input_maf_file=None,
        input_filter_tsv=None,
        set_prefix=None,
        min_length=0,
        methods=[],
    )

    (args) = E.start(parser, argv)

    if args.input_filter_tsv:
        with iotools.open_file(args.input_filter_tsv) as inf:
            skip_id = set([x[:-1] for x in inf])
    else:
        skip_id = False

    counter = E.Counter()

    if args.set_prefix:
        prefix = "s {}".format(args.set_prefix)
    else:
        prefix = None

    for block in iterate_maf_blocks(args.stdin):
        counter.blocks_input += 1
        if skip_id:
            if block[2].startswith("s "):
                id = re.match("s (\S+)", block[2]).groups()[0]
                if id in skip_id:
                    counter.blocks_skipped_id += 1
                    continue

        if args.min_length:
            if block[2].startswith("s "):
                id, pos, length = re.match("s (\S+)\s+(\d+)\s+(\d+)", block[2]).groups()
                if int(length) <= args.min_length:
                    counter.blocks_skipped_length += 1
                    continue

        if prefix:
            block[2] = prefix + block[2][4:]

        if block[2].startswith("s "):
            header, ali1, ali2, qual = parse_block(block)
            if "shift-region" in args.methods:
                rows = []
                contig, start, end = parse_region_string(ali1.src)
                ali1 = ali1._replace(src=contig, start=start + ali1.start)
                rows.append(list(map(str, ali1)))
                rows.append(list(map(str, ali2)))
                if qual:
                    rows.append(list(map(str, qual)))
                lines = [header]
                lines.append(format_tabular(rows, "llrrrrl"))
                lines.append("\n")
                block = lines
        counter.blocks_output += 1
        args.stdout.write("".join(block))

    E.info(counter)
    E.stop()
