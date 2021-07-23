'''
gff2fasta.py - output sequences from genomic features
=====================================================

:Tags: Genomics Intervals Sequences GFF Fasta Transformation

Purpose
-------

This script outputs the genomic sequences for intervals within
a :term:`gff` or :term: `gtf` formatted file.

The ouput can be optionally masked and filtered.

Usage
-----

If you want to convert a ``features.gff`` file with intervals information
into a :term:`fasta` file containing the sequence of each interval, use this
script as follows::

   python gff2fasta.py --genome-file=hg19 < features.gff > features.fasta

The input can also be a :term:`gtf` formatted file. In that case, use the
``--is-gtf`` option::

   python gff2fasta.py --genome-file=hg19 --is-gtf < features.gtf >\
 features.fasta

If you want to add a polyA tail onto each transcript you can use the `extend`
options:

   python gff2fasta.py --genome-file=hg19 --is-gtf
   --extend-at=3 --extend-by=125 --extend-with=A
   < features.gtf > features.fasta

If you want to merge the sequence of similar features together, please use
``--merge-overlapping``::

   python gff2fasta.py --genome-file=hg19 --merge-overlapping < features.gff >\
 features.fasta

It is possible to filter the output by selecting a minimum or maximum number
of nucleotides in the resultant fasta sequence with ``--max-length`` or
``--min-interval-length`` respectively::

   python gff2fasta.py --genome-file=hg19 --max-length=100\
 < features.gff > features.fasta

Or you can also filter the output by features name with the ``--feature``
option::

   python gff2fasta.py --genome-file=hg19 --feature=exon < features.gff\
 > features.fasta

On the other hand, low-complexity regions can be masked with the ``--masker``
option and a given :term:`gff` formatted file::

   python gff2fasta.py --genome-file=hg19 --masker=dust\
 --maskregions-bed-file=intervals.gff < features.gff > features.fasta

where ``--masker`` can take the following values: ``dust``, ``dustmasker``,
and ``softmask``.

Options
-------

``--is-gtf``
  Tells the script to expect a :term:`gtf` format file

``--genome-file``
  PATH to Fasta file of genome build to use

``--merge-overlapping``
  Merge features in :term:`gtf`/:term:`gff` file that are adjacent and share
  attributes

``--method=filter --filter-method``
  Filter on a :term:`gff` feature such as ``exon`` or ``CDS``

``--maskregions-bed-file``
  Mask sequences in intervals in :term:`gff` file

``--remove-masked-regions``
  Remove sequences in intervals in :term:`gff` file rather than masking them

``--min-interval-length``
  Minimum output sequence length

``--max-length``
  Maximum output sequence length

``--extend-at``
  Extend sequence at 3', 5' or both end.  Optionally '3only' or '5only' will
  return only the 3' or 5' extended sequence

``--extend-by``
  Used in conjunction with ``--extend-at``, the number of nucleotides to extend
  by

``--extend-with``
  Optional. Used in conjunction with ``--extend-at`` and ``--extend-by``.
  Instead of extending by the genomic sequence, extend by this string repeated
  n times, where n is --entend-by


``--masker``
  Masker type to use: dust, dustmasker, soft or none

``--fold-at``
  Fold the fasta sequence every n bases

``--naming-attribute``
  Use this attribute to name the fasta entries

Command line options
--------------------
'''

import sys
import quicksect
import cgatcore.experiment as E
import cgat.GTF as GTF
import cgat.Genomics as Genomics
import cgatcore.iotools as iotools
import cgat.IndexedFasta as IndexedFasta
import cgat.Intervals as Intervals
import cgat.Masker as Masker


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("--is-gtf", dest="is_gtf", action="store_true",
                        help="input is gtf instead of gff.")

    parser.add_argument("-g", "--genome-file", dest="genome_file", type=str,
                        help="filename with genome.")

    parser.add_argument(
        "-m", "--merge-adjacent", dest="merge", action="store_true",
        help="merge adjacent intervals with the same attributes.")

    parser.add_argument(
        "-e", "--feature", dest="feature", type=str,
        help="filter by a feature, for example 'exon', 'CDS'."
        " If set to the empty string, all entries are output "
        ".")

    parser.add_argument(
        "-f", "--maskregions-bed-file", dest="filename_masks",
        type=str, metavar="gff",
        help="mask sequences with regions given in gff file "
        ".")

    parser.add_argument(
        "--remove-masked-regions", dest="remove_masked_regions",
        action="store_true",
        help="remove regions instead of masking .")

    parser.add_argument(
        "--min-interval-length", dest="min_length", type=int,
        help="set minimum length for sequences output "
        )

    parser.add_argument(
        "--max-length", dest="max_length", type=int,
        help="set maximum length for sequences output "
        )

    parser.add_argument(
        "--extend-at", dest="extend_at", type=str,
        choices=("none", "3", "5", "both", "3only", "5only"),
        help="extend at no end, 3', 5' or both ends. If "
        "3only or 5only are set, only the added sequence "
        "is returned")

    parser.add_argument(
        "--header-attributes", dest="header_attr",
        action="store_true",
        help="add GFF entry attributes to the FASTA record"
        " header section")

    parser.add_argument(
        "--extend-by", dest="extend_by", type=int,
        help="extend by # bases")

    parser.add_argument(
        "--extend-with", dest="extend_with", type=str,
        help="extend using base")

    parser.add_argument(
        "--masker", dest="masker", type=str,
        choices=("dust", "dustmasker", "softmask", "none"),
        help="apply masker .")

    parser.add_argument(
        "--fold-at", dest="fold_at", type=int,
        help="fold sequence every n bases.")

    parser.add_argument(
        "--fasta-name-attribute", dest="naming_attribute", type=str,
        help="use attribute to name fasta entry. Currently only compatable"
        " with gff format .")

    parser.set_defaults(
        is_gtf=False,
        genome_file=None,
        merge=False,
        feature=None,
        filename_masks=None,
        remove_masked_regions=False,
        min_length=0,
        max_length=0,
        extend_at=None,
        extend_by=100,
        extend_with=None,
        masker=None,
        fold_at=None,
        naming_attribute=False,
        header_attr=False,
    )

    (args) = E.start(parser)

    if args.genome_file:
        fasta = IndexedFasta.IndexedFasta(args.genome_file)
        contigs = fasta.getContigSizes()

    if args.is_gtf:
        iterator = GTF.transcript_iterator(GTF.iterator(args.stdin))
    else:
        gffs = GTF.iterator(args.stdin)
        if args.merge:
            iterator = GTF.joined_iterator(gffs)
        else:
            iterator = GTF.chunk_iterator(gffs)

    masks = None
    if args.filename_masks:
        masks = {}
        with iotools.open_file(args.filename_masks, "r") as infile:
            e = GTF.readAsIntervals(GTF.iterator(infile))

        # convert intervals to intersectors
        for contig in list(e.keys()):
            intersector = quicksect.IntervalTree()
            for start, end in e[contig]:
                intersector.add(start, end)
            masks[contig] = intersector

    ninput, noutput, nmasked, nskipped_masked = 0, 0, 0, 0
    nskipped_length = 0
    nskipped_noexons = 0

    feature = args.feature

    # iterator is a list containing groups (lists) of features.
    # Each group of features have in common the same transcript ID, in case of
    # GTF files.
    for ichunk in iterator:

        ninput += 1

        if feature:
            chunk = [x for x in ichunk if x.feature == feature]
        else:
            chunk = ichunk

        if len(chunk) == 0:
            nskipped_noexons += 1
            E.info("no features in entry from "
                   "%s:%i..%i - %s" % (ichunk[0].contig,
                                       ichunk[0].start,
                                       ichunk[0].end,
                                       str(ichunk[0])))
            continue

        contig, strand = chunk[0].contig, chunk[0].strand

        if args.is_gtf:
            name = chunk[0].transcript_id
        else:
            if args.naming_attribute:
                attr_dict = {x.split("=")[0]: x.split("=")[1]
                             for x in chunk[0].attributes.split(";")}
                name = attr_dict[args.naming_attribute]
            else:
                name = str(chunk[0].attributes)

        lcontig = contigs[contig]
        positive = Genomics.IsPositiveStrand(strand)
        intervals = [(x.start, x.end) for x in chunk]
        intervals.sort()

        if masks:
            if contig in masks:
                masked_regions = []
                for start, end in intervals:
                    masked_regions += [(x.start, x.end)
                                       for x in masks[contig].find(quicksect.Interval(start, end))]

                masked_regions = Intervals.combine(masked_regions)
                if len(masked_regions):
                    nmasked += 1

                if args.remove_masked_regions:
                    intervals = Intervals.truncate(intervals, masked_regions)
                else:
                    raise NotImplementedError("unimplemented")

                if len(intervals) == 0:
                    nskipped_masked += 1
                    if args.loglevel >= 1:
                        args.stdlog.write("# skipped because fully masked: "
                                          "%s: regions=%s masks=%s\n" %
                                          (name,
                                           str([(x.start,
                                                 x.end) for x in chunk]),
                                           masked_regions))
                    continue

        out = intervals

        if args.extend_at and not args.extend_with:
            if args.extend_at == "5only":
                intervals = [(max(0, intervals[0][0] - args.extend_by),
                              intervals[0][0])]
            elif args.extend_at == "3only":
                intervals = [(intervals[-1][1],
                              min(lcontig,
                                  intervals[-1][1] + args.extend_by))]
            else:
                if args.extend_at in ("5", "both"):
                    intervals[0] = (max(0,
                                        intervals[0][0] - args.extend_by),
                                    intervals[0][1])
                if args.extend_at in ("3", "both"):
                    intervals[-1] = (intervals[-1][0],
                                     min(lcontig,
                                         intervals[-1][1] + args.extend_by))

        if not positive:
            intervals = [(lcontig - x[1], lcontig - x[0])
                         for x in intervals[::-1]]
            out.reverse()

        s = [fasta.getSequence(contig, strand, start, end)
             for start, end in intervals]
        # IMS: allow for masking of sequences
        s = Masker.maskSequences(s, args.masker)
        l = sum([len(x) for x in s])
        if (l < args.min_length or
                (args.max_length and l > args.max_length)):
            nskipped_length += 1
            if args.loglevel >= 1:
                args.stdlog.write("# skipped because length out of bounds "
                                  "%s: regions=%s len=%i\n" %
                                  (name, str(intervals), l))
                continue

        if args.extend_at and args.extend_with:
            extension = "".join((args.extend_with,) * args.extend_by)

            if args.extend_at in ("5", "both"):
                s[1] = extension + s[1]
            if args.extend_at in ("3", "both"):
                s[-1] = s[-1] + extension

        if args.fold_at:
            n = args.fold_at
            s = "".join(s)
            seq = "\n".join([s[i:i+n] for i in range(0, len(s), n)])
        else:
            seq = "\n".join(s)

        if args.header_attr:
            attributes = " ".join([":".join([ax, ay]) for ax, ay in chunk[0].asDict().items()])
            args.stdout.write(">%s %s:%s:%s feature:%s %s\n%s\n" % (name,
                                                                    contig,
                                                                    strand,
                                                                    ";".join(["%i-%i" %
                                                                              x for x in out]),
                                                                    chunk[0].feature,
                                                                    attributes,
                                                                    seq))
        else:
            args.stdout.write(">%s %s:%s:%s\n%s\n" % (name,
                                                      contig,
                                                      strand,
                                                      ";".join(["%i-%i" %
                                                                x for x in out]),
                                                      seq))

        noutput += 1

    E.info("ninput=%i, noutput=%i, nmasked=%i, nskipped_noexons=%i, "
           "nskipped_masked=%i, nskipped_length=%i" %
           (ninput, noutput, nmasked, nskipped_noexons,
            nskipped_masked, nskipped_length))

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
