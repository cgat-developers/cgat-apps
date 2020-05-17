"""Compare alignments in BAM Files
====================================

This script compares the alignments in a :term:`BAM` formatted file
against alignments in a reference (also :term:`BAM` formatted).

The script outputs a summary of missing alignments as well as a detailed
table containing for each read the following information:

.. csv-table::
   :header: "column_name", "description"

   "read", "name of the read"
   "length", "length of the read"
   "status", "status, one of (unique|multi)_(mapped|mismapped)"
   "overlap", "physical overlap between reads"
   "comp_contig", "contig that read maps to in test set",
   "comp_start", "chromosomal position that reads maps to"
   "comp_end", "chromosomal position that reads maps to"
   "ref_contig", "contig that read maps to in reference set"
   "ref_start", "chromosomal position that reads maps to"
   "ref_end", "chromosomal position that reads maps to"
   "shared_misaligned", "misaligned residue pairs between test and reference"
   "shared_aligned", "correctly aligned residue pairs shared between test and reference"
   "shared_insertion", "correctly aligned insertions"
   "shared_deletion", "correctly aligned deletions"
   "comp_aligned", "aligned pairs in test read"
   "comp_insertion", "inserted bases in test read"
   "comp_deletion", "deleted bases in test read"
   "ref_aligned", "aligned pairs in reference read"
   "ref_insertion", "inserted bases in reference read"
   "ref_deletion", "deleted bases in reference read"

.. note::

   This script assumes single-end data (for example, ONT).

"""
import sys
import os
import pysam
import re
import cgatcore.experiment as E


def count_pairs(s):
    insertions = len([x for x in s if x[1] is None])
    deletions = len([x for x in s if x[0] is None])
    paired = len(s) - insertions - deletions
    return paired, insertions, deletions


def iterate_read_pairs(stream1, stream2, qname_fn=None):

    read1 = next(stream1)
    read2 = next(stream2)

    def passthrough(x):
        return x

    if qname_fn is None:
        qname_fn = passthrough

    try:
        n1 = qname_fn(read1.query_name)
        n2 = qname_fn(read2.query_name)
        while 1:
            if n1 == n2:
                yield(read1, read2)
                read1 = next(stream1)
                read2 = next(stream2)
                n1 = qname_fn(read1.query_name)
                n2 = qname_fn(read2.query_name)
            elif n1 < n2:
                yield(read1, None)
                read1 = next(stream1)
                n1 = qname_fn(read1.query_name)
            else:
                yield(None, read2)
                read2 = next(stream2)
                n2 = qname_fn(read2.query_name)
    except StopIteration:
        pass


def group_pairs(stream):

    reads1, reads2 = [], None
    last_name = None

    for read1, read2 in stream:
        if read1 is None:
            name = read2.query_name
            if last_name == name:
                reads2 = read2
            else:
                if reads2 is not None:
                    yield reads1, reads2
                reads1, reads2 = [], read2
                last_name = name
        elif read2 is None:
            name = read1.query_name
            assert last_name == name
            reads1.append(read1)
        else:
            name = read1.query_name
            if name != last_name:
                if reads2 is not None:
                    yield reads1, reads2
                reads1, reads2 = [read1], read2
            else:
                assert False, "this should not happen"
            last_name = name

    yield reads1, reads2


def main(argv=sys.argv):

    parser = E.ArgumentParser(version="%prog version: $Id$",
                              usage=globals()["__doc__"])

    parser.add_argument(
        "-i", "--input-bam", dest="input_bam_file", type=str,
        help="input bam file")

    parser.add_argument(
        "-f", "--reference-bam", dest="reference_bam_file", type=str,
        help="reference BAM file ")

    parser.add_argument(
        "-q", "--query-name-regex", dest="query_name_regex", type=str,
        help="regular expression to apply on query name. "
        "Potentially required to match samtools sort order and should "
        "evaluate to an integer ")

    parser.set_defaults(
        input_bam_file=None,
        reference_bam_file=None,
        query_name_regex=None,
    )

    (options, args) = E.start(parser, argv, add_output_options=True)

    if len(args) == 2:
        options.input_bam_file = args[0]
        options.reference_bam_file = args[1]

    if options.input_bam_file is None:
        raise ValueError("please supply a BAM file as input")

    if options.reference_bam_file is None:
        raise ValueError("please supply a BAM file as reference")

    # update paths to absolute
    options.input_bam_file = os.path.abspath(options.input_bam_file)
    options.reference_bam_file = os.path.abspath(options.reference_bam_file)

    if not os.path.exists(options.input_bam_file):
        raise OSError("input bam file {} does not exist".format(
            options.input_bam_file))

    if not os.path.exists(options.reference_bam_file):
        raise OSError("reference bam file {} does not exist".format(
            options.reference_bam_file))

    bam_in = pysam.AlignmentFile(options.input_bam_file)
    ref_in = pysam.AlignmentFile(options.reference_bam_file)

    outf_mapped = E.open_output_file("mapped")
    outf_mapped.write("\t".join(
        ["read",
         "length",
         "status",
         "overlap",
         "comp_contig",
         "comp_start",
         "comp_end",
         "ref_contig",
         "ref_start",
         "ref_end",
         "shared_misaligned",
         "shared_aligned",
         "shared_insertion",
         "shared_deletion",
         "comp_aligned",
         "comp_insertion",
         "comp_deletion",
         "ref_aligned",
         "ref_insertion",
         "ref_deletion"]) + "\n")

    outf_missing = E.open_output_file("missing")
    outf_missing.write("\t".join(
        ["read", "length", "status", "aligned",
         "insertion", "deletion"]) + "\n")

    counter = E.Counter()

    if options.query_name_regex:
        rx = re.compile(options.query_name_regex)

    def extract_query(x):
        return int(rx.search(x).groups()[0])

    qname_fn = None
    if options.query_name_regex:
        qname_fn = extract_query

    for reads_cmp, read_ref in group_pairs(iterate_read_pairs(
            bam_in.fetch(until_eof=True),
            ref_in.fetch(until_eof=True),
            qname_fn=qname_fn)):

        if len(reads_cmp) == 0:
            counter.missing += 1
            pairs_ref = set(read_ref.get_aligned_pairs())
            outf_missing.write("\t".join(
                map(str, (
                    read_ref.query_name,
                    read_ref.query_length,
                    "missing") +
                    count_pairs(pairs_ref))) + "\n")
            continue

        if len(reads_cmp) > 1:
            # multiple matches
            counter.multi_mapping += 1
            prefix = "multi_"
        else:
            counter.unique_mapping += 1
            prefix = "unique_"

        is_mapped = False
        for read_cmp in reads_cmp:

            counter.paired += 1

            if read_cmp.is_unmapped:
                counter.unmapped += 1
                pairs_ref = set(read_ref.get_aligned_pairs())
                outf_missing.write("\t".join(
                    map(str, (
                        read_ref.query_name,
                        read_ref.query_length,
                        "unmapped") +
                        count_pairs(pairs_ref))) + "\n")
                continue

            overlap = max(0, (min(read_cmp.reference_end,
                                  read_ref.reference_end) -
                              max(read_cmp.reference_start,
                                  read_ref.reference_start)))

            pairs_cmp = set(read_cmp.get_aligned_pairs())
            pairs_ref = set(read_ref.get_aligned_pairs())
            shared_cmp = pairs_cmp.intersection(pairs_ref)
            unique_cmp = pairs_cmp.difference(pairs_ref)
            missaligned = len([x for x, y in unique_cmp
                               if x is not None and y is not None])

            if read_cmp.reference_name != read_ref.reference_name or \
               overlap == 0:
                status = "mismapped"
            else:
                counter.overlap += 1
                status = "mapped"
                is_mapped = True

            outf_mapped.write("\t".join(
                map(str, (read_cmp.query_name,
                          read_cmp.query_length,
                          prefix + status,
                          overlap,
                          read_cmp.reference_name,
                          read_cmp.reference_start,
                          read_cmp.reference_end,
                          read_ref.reference_name,
                          read_ref.reference_start,
                          read_ref.reference_end,
                          missaligned) +
                    count_pairs(shared_cmp) +
                    count_pairs(pairs_cmp) +
                    count_pairs(pairs_ref))) + "\n")
        else:
            if is_mapped:
                status = "mapped"
            else:
                status = "mismapped"

            counter[prefix + status] += 1

    with E.open_output_file("summary") as outf:
        outf.write("category\tcounts\n")
        outf.write(counter.asTable() + "\n")

    E.stop()

if __name__ == "__main__":
    sys.exit(main())
