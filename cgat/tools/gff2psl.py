'''
gff2psl.py - convert from gff to psl
====================================

:Tags: Genomics Intervals GFF PSL Conversion

Note: This script is scheduled for deprecation in May 2024.

Purpose
-------

This scripts converts from a :term:`gff` formatted
file to a :term:`psl` formatted file.
The output can be modified by the following command line options:

--allow-duplicates
    keep duplicate entries from gff/gtf input file

--genome-file
    restrict output to gff/gtf entries with contigs in fasta file

--queries-tsv-file
    restrict output to queries in fasta file

Usage
-----

Example::

   python gff2psl.py < in.gff > out.psl

Type::

   python gff2psl.py --help

for command line help.
genome-file

Command line options
--------------------
'''

import sys
import cgatcore.experiment as E
import cgat.IndexedFasta as IndexedFasta
import cgat.Blat as Blat
import cgat.GTF as GTF
import alignlib_lite
import cgat.Intervals as Intervals


def print_deprecation_warning():
    warning_message = ("""WARNING: 'gff2psl.py' is deprecated and will be removed in May 2024.""")
    print(warning_message, file=sys.stderr)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    print_deprecation_warning()

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("--is-gtf", dest="is_gtf", action="store_true",
                        help="input is gtf.")

    parser.add_argument("--no-header", dest="with_header", action="store_false",
                        help="do not output BLAT header.")

    parser.add_argument("-g", "--genome-file", dest="genome_file", type=str,
                        help="filename with genome.")

    parser.add_argument("--queries-tsv-file", dest="input_filename_queries", type=str,
                        help="fasta filename with queries.")

    parser.add_argument("--allow-duplicates", dest="allow_duplicates", action="store_true",
                        help="""permit duplicate entries. Adjacent exons of a transcript will still be merged.""")

    parser.set_defaults(is_gtf=False,
                        genome_file=None,
                        with_header=True,
                        allow_duplicates=False,
                        test=None)

    (args) = E.start(parser, add_pipe_options=True)

    if args.genome_file:
        genome_fasta = IndexedFasta.IndexedFasta(args.genome_file)
    else:
        genome_fasta = None

    if args.input_filename_queries:
        queries_fasta = IndexedFasta.IndexedFasta(
            args.input_filename_queries)
    else:
        queries_fasta = None

    ninput, noutput, nskipped = 0, 0, 0

    if args.is_gtf:
        iterator = GTF.transcript_iterator(GTF.iterator_filtered(GTF.iterator(sys.stdin),
                                                                 feature="exon"),
                                           strict=not args.allow_duplicates)
    else:
        iterator = GTF.joined_iterator(GTF.iterator(sys.stdin))

    if args.with_header:
        args.stdout.write(Blat.Match().getHeader() + "\n")

    for gffs in iterator:

        if args.test and ninput >= args.test:
            break

        ninput += 1

        result = alignlib_lite.py_makeAlignmentBlocks()

        xstart = 0

        intervals = Intervals.combine([(gff.start, gff.end) for gff in gffs])

        for start, end in intervals:
            xend = xstart + end - start

            result.addDiagonal(xstart, xend,
                               start - xstart)
            xstart = xend

        entry = Blat.Match()
        entry.mQueryId = gffs[0].transcript_id
        entry.mSbjctId = gffs[0].contig
        entry.strand = gffs[0].strand

        if genome_fasta:
            if entry.mSbjctId in genome_fasta:
                entry.mSbjctLength = genome_fasta.getLength(entry.mSbjctId)
            else:
                entry.mSbjctLength = result.getColTo()

        if queries_fasta:
            if entry.mQueryId in queries_fasta:
                entry.mQueryLength = queries_fasta.getLength(entry.mQueryId)
        else:
            entry.mQueryLength = result.getRowTo()

        entry.fromMap(result)

        args.stdout.write(str(entry) + "\n")
        noutput += 1

    E.info("ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped))

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
