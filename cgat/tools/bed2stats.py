'''bed2stats.py - summary of bed file contents
============================================

:Tags: Genomics Intervals Summary BED

Purpose
-------

This script takes a :term:`bed`-formatted file as input and outputs the number
of intervals and bases in the bed file. Counts can be subdivided by setting
the ``--aggregate-by`` command line option:

``contig``
   output counts per contig (column 1)

``name``
   output counts grouped by the name field in the :term:`bed` formatted
   file (column 4)

``track``
   output counts per track in the :term:`bed` formatted file.

Note that a count of bases usually makes only sense if the intervals
submitted are non-overlapping.

If the option --add-percent is given, an additional column will output
the percent of the genome covered by intervals. This requires a
--genome-file to be given as well.

Usage
-----

To count the number of intervals, type::

   cgat bed2table < in.bed

+-----+--------+----------+------+
|track|ncontigs|nintervals|nbases|
+-----+--------+----------+------+
|all  |23      |556       |27800 |
+-----+--------+----------+------+

To count per contig::

   cgat bed2table --aggregate=contig < in.bed

+-----+--------+----------+------+
|track|ncontigs|nintervals|nbases|
+-----+--------+----------+------+
|chrX |1       |11        |550   |
+-----+--------+----------+------+
|chr13|1       |12        |600   |
+-----+--------+----------+------+
|chr12|1       |37        |1850  |
+-----+--------+----------+------+
|...  |...     |...       |...   |
+-----+--------+----------+------+

Type::

   cgat bed2table --help

for command line help.

Command line options
--------------------

'''
import sys
import collections
import cgat.Bed as Bed
import cgatcore.experiment as E
import cgat.IndexedFasta as IndexedFasta


class Counter:

    headers = ["ncontigs", "nintervals", "nbases"]
    headers_percent = ["ncontigs", "nintervals", "nbases", "pbases"]

    def __init__(self):
        self.intervals_per_contig = collections.defaultdict(int)
        self.bases_per_contig = collections.defaultdict(int)
        self.size = None

    def setSize(self, size):
        self.size = size

    def add(self, bed):
        self.intervals_per_contig[bed.contig] += 1
        self.bases_per_contig[bed.contig] += bed.end - bed.start

    def __str__(self):
        bases = sum(self.bases_per_contig.values())
        if self.size is None:
            return "%i\t%i\t%i" % (len(self.intervals_per_contig),
                                   sum(self.intervals_per_contig.values()),
                                   bases,
                                   )
        else:
            return "%i\t%i\t%i\t%5.2f" % (len(self.intervals_per_contig),
                                          sum(self.intervals_per_contig.values(
                                          )),
                                          sum(self.bases_per_contig.values()),
                                          100.0 * bases / self.size
                                          )


def main(argv=None):

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument(
        "-g", "--genome-file", dest="genome_file", type=str,
        help="filename with genome.")

    parser.add_argument(
        "-a", "--aggregate-by", dest="aggregate", type=str,
        choices=("name", "contig", "track", "none"),
        help="aggregate counts by feature.")

    parser.add_argument(
        "-p", "--add-percent", dest="add_percent", action="store_true",
        help="add percentages.")

    parser.set_defaults(
        genome_file=None,
        aggregate="none",
        add_percent=False,
    )

    (args) = E.start(parser, argv)

    # get files
    if args.genome_file:
        fasta = IndexedFasta.IndexedFasta(args.genome_file)
    else:
        if args.add_percent:
            raise ValueError("--add-percent option requires --genome-file")
        fasta = None

    if args.add_percent and not args.aggregate == "contig":
        raise NotImplementedError(
            "--add-percent option requires --aggregate=contig")

    counts = collections.defaultdict(Counter)
    total = Counter()
    output_totals = True

    if args.aggregate == "track":
        keyf = lambda x: x.track
    elif args.aggregate == "name":
        keyf = lambda x: x.name
    elif args.aggregate == "contig":
        keyf = lambda x: x.contig
    else:
        keyf = lambda x: "all"
        output_totals = False

    for bed in Bed.iterator(args.stdin):
        counts[keyf(bed)].add(bed)
        total.add(bed)

    outf = args.stdout

    key = "track"
    if args.add_percent:
        outf.write("%s\t%s\n" % (key, "\t".join(Counter.headers_percent)))
    else:
        outf.write("%s\t%s\n" % (key, "\t".join(Counter.headers)))

    total_bases = 0
    for key, count in sorted(counts.items()):
        if args.add_percent:
            total_bases += fasta.getLength(key)
            count.setSize(fasta.getLength(key))

        outf.write("%s\t%s\n" % (key, str(count)))

    if output_totals:
        if args.add_percent:
            count.setSize(total_bases)
        outf.write("%s\t%s\n" % ("total", str(total)))
    E.stop()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
