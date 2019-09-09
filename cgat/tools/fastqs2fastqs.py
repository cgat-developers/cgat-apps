'''fastqs2fastqs.py - manipulate (merge/reconcile) fastq files
=============================================================

:Tags: Genomics NGS FASTQ FASTQ Manipulation

Purpose
-------

This script manipulates multiple fastq files and outputs
new fastq files. Currently only the method ``reconcile``
is implemented.

reconcile
+++++++++

Reconcile reads from a pair of fastq files.

This method takes two fastq files and outputs two fastq files such
that all reads in the output are present in both output files.

The typical use case is that two fastq files containing the first and
second part of a read pair have been independently filtered, for
example by quality scores, truncation, etc. As a consequence some
reads might be missing from one file but not the other. The reconcile
method will output two files containing only reads that are common to
both files.

The two files must be sorted by read identifier.

Example input, read2 and read3 are only present in either of the
files:

   # File1        # File 2

   @read1         @read1
   AAA            AAA
   +              +
   !!!            !!!
   @read2         @read3
   CCC            TTT
   +              +
   !!!            !!!
   @read4         @read4
   GGG            GGG
   +              +
   !!!            !!!

Example output, only the reads common to both files are output::

   # File1        # File 2

   @read1         @read1
   AAA            AAA
   +              +
   !!!            !!!
   @read4         @read4
   GGG            GGG
   +              +
   !!!            !!!

Usage
-----

Example::

   python fastqs2fastqs.py \
            --method=reconcile \
            --output-filename-pattern=myReads_reconciled.%s.fastq \
            myReads.1.fastq.gz myReads.2.fastq.gz

In this example we take a pair of fastq files, reconcile by read
identifier and output 2 new fastq files named
``myReads_reconciled.1.fastq.gz`` and
``myReads_reconciled.2.fastq.gz``.

Type::

   python fastqs2fastqs.py --help

for command line help.

Command line options
--------------------

'''

import sys
import re
import pysam

import cgatcore.iotools as iotools
import cgatcore.experiment as E
import cgat.FastqTools as fastqtools


class PatternGetter:

    def __init__(self, pattern):
        self.pattern = re.compile(pattern)

    def __call__(self, id):
        return self.pattern.search(id).groups()[0]


def plain_getter(id):
    return id


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("-m", "--method", dest="method", type=str,
                        choices=('reconcile', 'filter-by-sequence'),
                        help="method to apply.")

    parser.add_argument(
        "-c", "--chop-identifier", dest="chop", action="store_true",
        help="whether or not to trim last character of the  "
        "sequence name. For example sometimes ids in the first "
        "file in the pair will end with \1 and the second "
        "with \2. If --chop-identifier is not specified "
        "then the results will be wrong.")

    parser.add_argument(
        "-u", "--unpaired", dest="unpaired", action="store_true",
        help="whether or not to write out unpaired reads "
        "to a separate file")

    parser.add_argument(
        "--id-pattern-1", dest="id_pattern_1",
        help="If specified will use the first group from the"
        "pattern to determine the ID for the first read")

    parser.add_argument(
        "--id-pattern-2", dest="id_pattern_2",
        help="As above but for read 2")

    parser.add_argument(
        "--input-filename-fasta",
        dest="input_filename_fasta", type=str,
        help="input filename of FASTA formatted sequence "
        "for method 'filter-by-sequence'.")

    parser.add_argument(
        "--filtering-kmer-size",
        dest="filtering_kmer_size", type=int,
        help="kmer size for method 'filter-by-sequence'.")

    parser.add_argument(
        "--filtering-min-kmer-matches",
        dest="filtering_min_kmer_matches", type=int,
        help="minimum number of matches 'filter-by-sequence'.")

    parser.set_defaults(
        method="reconcile",
        chop=False,
        unpaired=False,
        input_filename_fasta=None,
        filtering_kmer_size=10,
        filtering_min_kmer_matches=20
    )

    # add common options (-h/--help, ...) and parse command line
    (args, unknown) = E.start(parser,
                              argv=argv,
                              add_output_options=True,
                              unknowns=True)

    if len(unknown) != 2:
        raise ValueError(
            "please supply at least two fastq files on the commandline")

    fn1, fn2 = unknown
    counter = E.Counter()

    if args.id_pattern_1:
        id1_getter = PatternGetter(args.id_pattern_1)
    else:
        id1_getter = plain_getter

    if args.id_pattern_2:
        id2_getter = PatternGetter(args.id_pattern_2)
    else:
        id2_getter = plain_getter

    if args.method == "reconcile":

        # IMS: switching to no store second set of read names and only use
        # lazily. Since generators don't have a size must keep track
        id_lengths = {fn1: 0, fn2: 0}

        def getIds(infile, id_getter=plain_getter):
            '''return ids in infile.'''
            aread = infile.readline
            while True:
                l = [aread().rstrip("\r\n") for i in range(4)]
                if not l[0]:
                    break
                r = id_getter(l[0].split()[0])
                # decide if to chop read number off
                id_lengths[infile.name] += 1
                if args.chop:
                    yield r[:-1]
                else:
                    yield r

        def write(outfile, infile, take, unpaired_file=None,
                  id_getter=plain_getter):
            '''filter fastq files with ids in take.'''
            aread = infile.readline
            while True:
                l = [aread().rstrip("\r\n") for i in range(4)]
                if not l[0]:
                    break
                r = id_getter(l[0].split()[0])
                if args.chop:
                    r = r[:-1]
                if r not in take:
                    if unpaired_file is None:
                        continue
                    else:
                        unpaired_file.write("\n".join(l) + "\n")
                else:
                    outfile.write("\n".join(l) + "\n")

        E.info("reading first in pair")
        inf1 = iotools.open_file(fn1)
        ids1 = set(getIds(inf1, id1_getter))

        E.info("reading second in pair")
        inf2 = iotools.open_file(fn2)
        # IMS: No longer keep as a set, but lazily evaluate into intersection
        # leads to large memory saving for large inf2, particularly if
        # inf1 is small.
        ids2 = getIds(inf2, id2_getter)
        take = ids1.intersection(ids2)

        E.info("first pair: %i reads, second pair: %i reads, "
               "shared: %i reads" %
               (id_lengths[fn1],
                id_lengths[fn2],
                len(take)))

        if args.unpaired:
            unpaired_filename = E.open_output_file(
                "unpaired.fastq.gz", "w")
        else:
            unpaired_filename = None

        with E.open_output_file("1", "w") as outf:
            inf = iotools.open_file(fn1)
            E.info("writing first in pair")
            write(outf, inf, take, unpaired_filename, id1_getter)

        with E.open_output_file("2", "w") as outf:
            inf = iotools.open_file(fn2)
            E.info("writing second in pair")
            write(outf, inf, take, unpaired_filename, id2_getter)

        counter.output = len(take)

        if args.unpaired:
            unpaired_filename.close()

    elif args.method == "filter-by-sequence":

        with pysam.FastxFile(args.input_filename_fasta) as inf:
            for record in inf:
                query_sequence = record.sequence
                break

        with pysam.FastxFile(fn1, persist=False) as inf1, \
                pysam.FastxFile(fn2, persist=False) as inf2, \
                E.open_output_file("matched.fastq.1.gz", "w") as outf_matched1, \
                E.open_output_file("matched.fastq.2.gz", "w") as outf_matched2, \
                E.open_output_file("unmatched.fastq.1.gz", "w") as outf_unmatched1, \
                E.open_output_file("unmatched.fastq.2.gz", "w") as outf_unmatched2:
            counter = fastqtools.filter_by_sequence(
                query_sequence,
                inf1,
                inf2,
                outf_matched1,
                outf_matched2,
                outf_unmatched1,
                outf_unmatched2,
                kmer_size=args.filtering_kmer_size,
                min_kmer_matches=args.filtering_min_kmer_matches)
        args.stdout.write(
            "\t".join(("input", "matched", "unmatched", "percent_matched")) + "\n")

        args.stdout.write(
            "\t".join(map(str, (
                counter.input, counter.matched, counter.unmatched,
                100.0 * counter.matched / counter.input))) + "\n")

    E.info(str(counter))
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
