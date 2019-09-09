'''
fastq2fastq.py - manipulate fastq files
=============================================

:Tags: Genomics NGS Sequences FASTQ Manipulation

Purpose
-------

This script performs manipulations on :term:`fastq` formatted
files. For example it can be used to change the quality score format
or sample a subset of reads.

The script predominantly is used for manipulation of single fastq
files. However, for some of its functionality it will take paired data
using the ``--pair-fastq-file`` and ``--output-filename-pattern`` options.
This applies to the ``sample`` and ``sort`` methods.

Usage
-----

Example::
  In this example we randomly sample 50% of reads from paired data provided in
  two :term:`fastq` files.

   head in.fastq.1

   @SRR111956.1 HWUSI-EAS618:7:1:27:1582 length=36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   +SRR111956.1 HWUSI-EAS618:7:1:27:1582 length=36
   =@A@9@BAB@;@BABA?=;@@BB<A@9@;@2>@;??
   @SRR111956.2 HWUSI-EAS618:7:1:29:1664 length=36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCC
   +SRR111956.2 HWUSI-EAS618:7:1:29:1664 length=36
   =B@9@0>A<B=B=AAA?;*(@A>(@<=*9=9@BA>7
   @SRR111956.3 HWUSI-EAS618:7:1:38:878 length=36
   AGTGAGCAGGGAAACAATGTCTGTCTAAGAATTTGA

   head in.fastq.2

   +SRR111956.3 HWUSI-EAS618:7:1:38:878 length=36
   <?@BA?;A=@BA>;@@7###################
   @SRR111956.4 HWUSI-EAS618:7:1:38:1783 length=36
   ATTAGTATTATCCATTTATATAATCAATAAAAATGT
   +SRR111956.4 HWUSI-EAS618:7:1:38:1783 length=36
   ?ABBA2CCBBB2?=BB@C>=AAC@A=CBB#######
   @SRR111956.5 HWUSI-EAS618:7:1:39:1305 length=36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   +SRR111956.5 HWUSI-EAS618:7:1:39:1305 length=36
   AA>5;A>*91?=AAA@@BBA<B=?ABA>2>?A<BB@

   command-line::
     cat in.fastq.1 | python fastq2fastq.py
                      --method=sample --sample-size 0.5
                      --pair-fastq-file in.fastq.2
                      --output-filename-pattern out.fastq.2
                      > out.fastq.1

   head out.fastq.1
   @SRR111956.1 HWUSI-EAS618:7:1:27:1582 length=36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   +
   =@A@9@BAB@;@BABA?=;@@BB<A@9@;@2>@;??
   @SRR111956.2 HWUSI-EAS618:7:1:29:1664 length=36
   CCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCCCC
   +
   =B@9@0>A<B=B=AAA?;*(@A>(@<=*9=9@BA>7
   @SRR111956.3 HWUSI-EAS618:7:1:38:878 length=36
   AGTGAGCAGGGAAACAATGTCTGTCTAAGAATTTGA
   +
   <?@BA?;A=@BA>;@@7###################
   @SRR111956.4 HWUSI-EAS618:7:1:38:1783 length=36
   ATTAGTATTATCCATTTATATAATCAATAAAAATGT
   +
   ?ABBA2CCBBB2?=BB@C>=AAC@A=CBB#######

Options
-------

The following methods are implemented (``--method``).

``change-format``

    change the quality format to new format given as
    target-format. Options are ``sanger``,
  ``solexa``, ``phred64``, ``integer`` and ``illumina-1.8``

``sample``

    Sub-sample a fastq file. The size of the sample is set by
    --sample-size

``unique``

    Remove duplicate reads based on read name

``trim3``

    Trim a fixed number of nucleotides from the 3' end of reads.
    (see ``--num-bases``). Note that there are better tools for
   trimming.

``trim5``

    Trim a fixed number of nucleotides from the 5' end of reads.
    (see ``--num-bases``). Note that there are better tools for
   trimming.

``sort``

    Sort the fastq file by read name.

``renumber-reads``

    Rename the reads based on pattern given in ``--pattern-identifier``
    e.g. ``--pattern-identifier="read_%010i"``

Type::

   python fastq2fastq.py --help

for command line help.


Command line options
--------------------

'''
import collections
import sys
import os
import re
import random
import pysam
import numpy
import cgatcore.experiment as E
import cgatcore.iotools as iotools
import cgat.Fastq as Fastq
import cgat.Genomics as Genomics


def process_cgat(options):

    c = E.Counter()

    assert options.input_fastq_file == "-"

    if options.method == "change-format":
        for record in Fastq.iterate_convert(options.stdin,
                                            format=options.target_format,
                                            guess=options.guess_format):
            c.input += 1
            options.stdout.write("%s\n" % record)
            c.output += 1

    elif options.method == "grep":
        for record in Fastq.iterate(options.stdin):
            if re.match(options.grep_pattern, record.seq):
                options.stdout.write("%s\n" % record)

    elif options.method == "reverse-complement":
        for record in Fastq.iterate(options.stdin):
            record.seq = Genomics.complement(record.seq)
            record.quals = record.quals[::-1]
            options.stdout.write("%s\n" % record)

    elif options.method == "sample":
        sample_threshold = min(1.0, options.sample_size)

        random.seed(options.seed)

        if options.pair:
            if not options.output_filename_pattern:
                raise ValueError(
                    "please specify output filename pattern for "
                    "second pair (--output-filename-pattern)")

            outfile1 = options.stdout
            outfile2 = iotools.open_file(options.output_filename_pattern, "w")

            for record1, record2 in zip(
                    Fastq.iterate(options.stdin),
                    Fastq.iterate(iotools.open_file(options.pair))):
                c.input += 1
                if random.random() <= sample_threshold:
                    c.output += 1
                    outfile1.write("%s\n" % record1)
                    outfile2.write("%s\n" % record2)
        else:
            for record in Fastq.iterate(options.stdin):
                c.input += 1
                if random.random() <= sample_threshold:
                    c.output += 1
                    options.stdout.write("%s\n" % record)

    elif options.method == "apply":
        ids = set(iotools.read_list(iotools.open_file(options.apply)))

        for record in Fastq.iterate(options.stdin):
            c.input += 1
            if re.sub(" .*", "", record.identifier).strip() in ids:
                c.output += 1
                options.stdout.write("%s\n" % record)

    elif options.method == "trim3":
        trim3 = options.nbases
        for record in Fastq.iterate(options.stdin):
            c.input += 1
            record.trim(trim3)
            options.stdout.write("%s\n" % record)
            c.output += 1

    elif options.method == "trim5":
        trim5 = options.nbases
        for record in Fastq.iterate(options.stdin):
            c.input += 1
            record.trim5(trim5)
            options.stdout.write("%s\n" % record)
            c.output += 1

    elif options.method == "unique":
        keys = set()
        for record in Fastq.iterate(options.stdin):
            c.input += 1
            if record.identifier in keys:
                continue
            else:
                keys.add(record.identifier)
            options.stdout.write("%s\n" % record)
            c.output += 1

    # Need to change this to incorporate both pairs
    elif options.method == "sort":
        if not options.pair:
            # This is quicker for a single fastq file
            statement = "paste - - - - | sort -k1,1 -t ' ' | tr '\t' '\n'"
            os.system(statement)
        else:
            if not options.output_filename_pattern:
                raise ValueError(
                    "please specify output filename for second pair "
                    "(--output-filename-pattern)")
            E.warn(
                "consider sorting individual fastq files - "
                "this is memory intensive")
            entries1 = {}
            entries2 = {}

            for record1, record2 in zip(
                    Fastq.iterate(options.stdin),
                    Fastq.iterate(iotools.open_file(options.pair))):
                entries1[
                    record1.identifier[:-2]] = (record1.seq, record1.quals)
                entries2[
                    record2.identifier[:-2]] = (record2.seq, record2.quals)

            outfile1 = options.stdout
            outfile2 = iotools.open_file(options.output_filename_pattern, "w")
            assert len(set(entries1.keys()).intersection(
                set(entries2.keys()))) == len(entries1),\
                "paired files do not contain the same reads "\
                "need to reconcile files"

            for entry in sorted(entries1):
                outfile1.write("@%s/1\n%s\n+\n%s\n" %
                               (entry, entries1[entry][0], entries1[entry][1]))
                outfile2.write("@%s/2\n%s\n+\n%s\n" %
                               (entry, entries2[entry][0], entries2[entry][1]))

    elif options.method == "renumber-reads":
        id_count = 1
        for record in Fastq.iterate(options.stdin):
            record.identifier = options.renumber_pattern % id_count
            id_count += 1
            options.stdout.write("@%s\n%s\n+\n%s\n" %
                                 (record.identifier, record.seq, record.quals))
    return c


def process_daisy(options):

    filter_n = "filter-N" in options.methods

    filter_ont = "filter-ONT" in options.methods

    if "filter-identifier" in options.methods:
        if options.input_filter_tsv is None:
            raise ValueError("please set --input-filter-tsv for method filter-identifier")
        with iotools.open_file(options.input_filter_tsv) as inf:
            filter_identifier = set([x.split()[0].strip() for x in inf.readlines()])
    else:
        filter_identifier = False

    if options.output_removed_tsv:
        outf_removed_tsv = iotools.open_file(options.output_removed_tsv, "w")
    else:
        outf_removed_tsv = None

    if options.output_removed_fastq:
        outf_removed_fastq = iotools.open_file(options.output_removed_fastq, "w")
    else:
        outf_removed_fastq = None

    if options.set_prefix:
        prefix = "{}".format(options.set_prefix)
    else:
        prefix = None

    quality_offset = options.quality_offset
    counter = E.Counter()

    with pysam.FastxFile(options.input_fastq_file) as inf:
        for read in inf:
            counter.input += 1
            remove = False
            if filter_n:
                chars = collections.Counter(read.sequence)
                if "N" in chars and \
                   100.0 * chars["N"] / len(read.sequence) > options.max_percent_N:
                    remove = True
                    counter.filter_n += 1

            if filter_identifier:
                if read.name not in filter_identifier:
                    counter.filter_identifier += 1
                    remove = True

            if filter_ont:
                quals = read.get_quality_array()
                n = len(quals)
                if n < options.min_sequence_length or \
                        float(sum(quals)) / n < options.min_average_quality:
                    counter.remove_ont += 1
                    remove = True

            if remove:
                counter.removed += 1
                if outf_removed_tsv:
                    outf_removed_tsv.write(read.name + "\n")
                if outf_removed_fastq:
                    outf_removed_fastq.write(str(read) + "\n")
                continue

            if prefix:
                read.name = prefix + read.name[2:]

            if quality_offset:
                quals = numpy.array(read.get_quality_array())
                quals += quality_offset
                quals[quals < 0] = 0
                quals += 33
                # pysam fastq is read-only, so fudge it:
                # Note: not outputting description
                read = "@{}\n{}\n+\n{}".format(
                    read.name,
                    read.sequence,
                    "".join([chr(x) for x in quals]))

            counter.output += 1

            options.stdout.write(str(read) + "\n")

    if outf_removed_tsv:
        outf_removed_tsv.close()

    if outf_removed_fastq:
        outf_removed_fastq.close()

    if options.output_stats_tsv:
        with iotools.open_file(options.output_stats_tsv, "w") as outf:
            outf.write(counter.asTable(as_rows=False) + "\n")

    return counter


def main(argv=sys.argv):

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument(
        "-i", "--input-fastq-file", dest="input_fastq_file", type=str,
        help="input fastq file. "
        )

    parser.add_argument(
        "--output-removed-tsv", dest="output_removed_tsv", type=str,
        help="if given, sequence identifiers of removed sequences will "
        "be stored in this file ")

    parser.add_argument(
        "--output-stats-tsv", dest="output_stats_tsv", type=str,
        help="if given, output statistics will be written to this file. "
        )

    parser.add_argument(
        "--output-removed-fastq", dest="output_removed_fastq", type=str,
        help="if given, removed fastq records will "
        "be stored in this file ")

    parser.add_argument(
        "-m", "--method", dest="methods", action="append", type=str,
        choices=("filter-N",
                 "filter-identifier",
                 "filter-ONT",
                 "offset-quality",
                 "apply",
                 "change-format",
                 "renumber-reads",
                 "sample",
                 "sort",
                 "trim3",
                 "trim5",
                 "unique",
                 "reverse-complement",
                 "grep"),
        help="methods to apply ")

    parser.add_argument(
        "--set-prefix", dest="set_prefix", type=str,
        help="set sequence prefix ")

    parser.add_argument(
        "--input-filter-tsv", dest="input_filter_tsv", type=str,
        help="list of sequence ides to filter ")

    parser.add_argument(
        "--min-average-quality", dest="min_average_quality", type=float,
        help="minimum average quality ")

    parser.add_argument(
        "--min-sequence-length", dest="min_sequence_length", type=int,
        help="minimum sequence length ")

    parser.add_argument(
        "--quality-offset", dest="quality_offset", type=int,
        help="offset to modify quality values with ")

    parser.add_argument(
        "--target-format", dest="target_format", type=str,
        choices=('sanger', 'solexa', 'phred64', 'integer', 'illumina-1.8'),
        help="guess quality score format and set quality scores "
        "to format.")

    parser.add_argument(
        "--guess-format", dest="guess_format", type=str,
        choices=('sanger', 'solexa', 'phred64', 'integer', 'illumina-1.8'),
        help="quality score format to assume if ambiguous.")

    parser.add_argument(
        "--sample-size", dest="sample_size", type=float,
        help="proportion of reads to sample. "
        "Provide a proportion of reads to sample, e.g. 0.1 for 10%, "
        "0.5 for 50%, etc.")

    parser.add_argument(
        "--pair-fastq-file", dest="pair", type=str,
        help="if data is paired, filename with second pair. "
        "Implemented for sampling.")

    parser.add_argument(
        "--map-tsv-file", dest="map_tsv_file", type=str,
        help="filename with tab-separated identifiers mapping for "
        "method apply.")

    parser.add_argument(
        "--num-bases", dest="nbases", type=int,
        help="number of bases to trim.")

    parser.add_argument(
        "--seed", dest="seed", type=int,
        help="seed for random number generator.")

    parser.add_argument(
        "--pattern-identifier", dest="renumber_pattern", type=str,
        help="rename reads in file by pattern")

    parser.add_argument(
        "--grep-pattern", dest="grep_pattern", type=str,
        help="subset to reads matching pattern")

    parser.set_defaults(
        input_fastq_file="-",
        methods=[],
        change_format=None,
        guess_format=None,
        sample_size=0.1,
        nbases=0,
        pair=None,
        apply=None,
        seed=None,
        renumber_pattern="read_%010i",
        grep_pattern=".*",
        max_percent_N=10.0,
        set_prefix=None,
        output_removed_tsv=None,
        output_removed_fastq=None,
        output_stats_tsv=None,
        input_filter_tsv=None,
        min_average_quality=0,
        min_sequence_length=0,
        quality_offset=0,
    )

    (args, unknown) = E.start(parser,
                              argv,
                              add_output_options=True,
                              unknowns=True)

    if len(unknown) == 1:
        args.input_fastq_file = unknown[0]

    if len(args.methods) == 0:
        raise ValueError("no method specified, please use --method")

    # this script combines two scripts with different functionalities
    # TODO: to be sanitized
    if args.methods[0] in ["apply",
                           "change-format",
                           "renumber-reads",
                           "sample",
                           "sort",
                           "trim3",
                           "trim5",
                           "unique",
                           "reverse-complement",
                           "grep"]:
        args.method = args.methods[0]
        counter = process_cgat(args)
    else:
        counter = process_daisy(args)

    E.info(counter)
    E.stop()
