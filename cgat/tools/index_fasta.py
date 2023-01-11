'''index_fasta.py - Index fasta formatted files
============================================

:Tags: Genomics Sequences FASTA Manipulation

Purpose
-------

This script indexes one or more :term:`fasta` formatted files into a
database that can be used by other scripts in the cgat code collection
and :mod:`IndexedFasta` for quick access to a particular part of a sequence.
This is very useful for large genomic sequences.

By default, the database is itself a :term:`fasta` formatted file in
which all line breaks and other white space characters have been
removed.  Compression methods are available to conserve disk space,
though they do come at a performance penalty.

The script implements several indexing and compression methods.  The
default method uses no compression and builds a simple random access
index based on a table of absolute file positions.  The sequence is
stored in a plain fasta file with one line per sequence allowing to
extract a sequence segment by direct file positioning.

Alternatively, the sequence can be block-compressed using different
compression methods (gzip, lzo, bzip). These are mostly for research
purposes.

See also http://pypi.python.org/pypi/pyfasta for another
implementation.  Samtools provides similar functionality with the
``samtools faidx`` command and block compression has been implemented
in the `bgzip http://samtools.sourceforge.net/tabix.shtml>`_ tool.

The script permits supplying synonyms to the database index. For
example, setting ``--synonyms=chrM=chrMT`` will ensure that the
mitochondrial genome sequence is returned both for the keys ``chrM``
and ``chrMT``.

Examples
--------

Index a collection of fasta files in a compressed archive::

   python index_fasta.py oa_ornAna1_softmasked ornAna1.fa.gz > oa_ornAna1_softmasked.log

To retrieve a segment::

   python index_fasta.py --extract=chr5:1000:2000 oa_ornAna1_softmasked

Indexing from a tar file is possible::

   python index_fasta.py oa_ornAna1_softmasked ornAna1.tar.gz > oa_ornAna1_softmasked.log

Indexing from stdin requires to use the ``-`` place-holder::

   zcat ornAna1.fa.gz | python index_fasta.py oa_ornAna1_softmasked - > oa_ornAna1_softmasked.log

Usage
-----

Type::

   cgat index_genome DATABASE [SOURCE...|-] [OPTIONS]
   cgat index_genome DATABASE [SOURCE...|-] --compression=COMPRESSION --random-access-points=100000

To create indexed DATABASE from SOURCE. Supply - as SOURCE to read from stdin.
If the output is to be compressed, a spacing for the random access points must
be supplied.

Type::

   cgat index_genome DATABASE --extract=CONTIG:[STRAND]:START:END

To extract the bases on the STRAND strand, between START to END from
entry CONTIG, from DATABASE.

Command line options
--------------------

'''
import cgat.IndexedFasta as IndexedFasta
import cgatcore.experiment as E
import sys


def main(argv=None):

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument(
        "-e", "--extract", dest="extract", type=str,
        help="extract region for testing purposes. Format is "
        "contig:strand:from:to. "
        "The default coordinates are 0-based "
        "open/closed coordinates on both strands, but can be changed "
        "by --input-format. "
        "For example, 'chr1:+:10:12' will return "
        "bases 11 and 12 on chr1. Elements from the end of the "
        "string can be omitted. For example, 'chr1' will return "
        "all of chromosome 'chr1'.")

    input_format_choices = ("one-forward-open", "zero-both-open")
    parser.add_argument("-i", "--input-format", dest="input_format",
                        type=str,
                        choices=input_format_choices,
                        help="coordinate format of input. Valid choices are "
                        "%s. See --extract." %
                        ", ".join(input_format_choices))

    parser.add_argument(
        "-s", "--synonyms", dest="synonyms", type=str,
        help="list of synonyms. This is a comma separated with list "
        "of equivalence relations. For example, chrM=chrMT "
        "means that chrMT will refer to chrM and either "
        "can be used to retrieve a sequence ")

    group = parser.add_argument_group("Bencharking options")

    group.add_argument("-b", "--benchmark", dest="benchmark",
                       action="store_true",
                       help="benchmark time for read access ")
    group.add_argument("--benchmark-num-iterations",
                       dest="benchmark_num_iterations",
                       type=int,
                       help="number of iterations for benchmark ")
    group.add_argument("--benchmark-fragment-size",
                       dest="benchmark_fragment_size",
                       type=int,
                       help="benchmark: fragment size.")
    parser.add_argument_group(group)

    group = parser.add_argument_group("Validation options")

    group.add_argument("--verify", dest="verify", type=str,
                       help="verify against other database.")

    group.add_argument("--verify-iterations", dest="verify_num_iterations",
                       type=int,
                       help="number of iterations for verification ")
    parser.add_argument_group(group)

    file_format_choices = ("fasta", "auto", "fasta.gz", "tar", "tar.gz")
    parser.add_argument("--file-format", dest="file_format", type=str,
                        choices=file_format_choices,
                        help="file format of input. Supply if data comes "
                        "from stdin "
                        "Valid choices are \\%s." %
                        ", ".join(file_format_choices))

    parser.add_argument("-a", "--clean-sequence", dest="clean_sequence",
                        action="store_true",
                        help="remove X/x from DNA sequences - they cause "
                        "errors in exonerate.")

    parser.add_argument("--allow-duplicates", dest="allow_duplicates",
                        action="store_true",
                        help="allow duplicate identifiers. Further occurances "
                        "of an identifier are suffixed by an '_%%i' ")

    parser.add_argument("--regex-identifier", dest="regex_identifier",
                        type=str,
                        help="regular expression for extracting the "
                        "identifier from fasta description line.")

    parser.add_argument("--force-output", dest="force", action="store_true",
                        help="force overwriting of existing files ")

    translator_choices = ("solexa", "phred", "bytes", "range200")
    parser.add_argument("-t", "--translator", dest="translator", type=str,
                        choices=translator_choices,
                        help="translate numerical quality scores. "
                        "Valid choices are \\%s." %
                        ", ".join(translator_choices))

    group = parser.add_argument_group("Compression options")

    compression_choices = ("lzo", "zlib", "gzip", "dictzip", "bzip2", "debug")
    group.add_argument("-c", "--compression", dest="compression", type=str,
                       choices=compression_choices,
                       help="compress database, using specified compression "
                       "method. "
                       "Valid choices are %s, but depend on availability on the "
                       "system " % ", ".join(compression_choices))

    group.add_argument("--random-access-points", dest="random_access_points",
                       type=int,
                       help="set random access points every # number "
                       "of nucleotides for block compression schemes ")

    group.add_argument(
        "--compress-index", dest="compress_index",
        action="store_true",
        help="compress index. The default is to use a plain-text, "
        "human-readable index.")

    parser.add_argument_group(group)

    parser.set_defaults(
        extract=None,
        input_format="zero-both-open",
        benchmark_fragment_size=1000,
        benchmark_num_iterations=1000000,
        benchmark=False,
        compression=None,
        random_access_points=0,
        synonyms=None,
        verify=None,
        verify_num_iterations=100000,
        verify_fragment_size=100,
        clean_sequence=False,
        allow_duplicates=False,
        regex_identifier=None,
        compress_index=False,
        file_format="auto",
        force=False,
        translator=None)

    (args, unknown) = E.start(parser,
                              unknowns=True)

    if args.synonyms:
        synonyms = {}
        for x in args.synonyms.split(","):
            a, b = x.split("=")
            a = a.strip()
            b = b.strip()
            if a not in synonyms:
                synonyms[a] = []
            synonyms[a].append(b)
    else:
        synonyms = None

    if args.translator:
        if args.translator == "phred":
            args.translator = IndexedFasta.TranslatorPhred()
        elif args.translator == "solexa":
            args.translator = IndexedFasta.TranslatorSolexa()
        elif args.translator == "bytes":
            args.translator = IndexedFasta.TranslatorBytes()
        elif args.translator == "range200":
            args.translator = IndexedFasta.TranslatorRange200()
        else:
            raise ValueError("unknown translator %s" % args.translator)

    if args.extract:
        fasta = IndexedFasta.IndexedFasta(unknown[0])
        fasta.setTranslator(args.translator)
        converter = IndexedFasta.getConverter(args.input_format)

        contig, strand, start, end = IndexedFasta.parseCoordinates(
            args.extract)
        sequence = fasta.getSequence(contig, strand,
                                     start, end,
                                     converter=converter)
        args.stdout.write(">%s\n%s\n" %
                          (args.extract, sequence))

    elif args.benchmark:
        import timeit
        timer = timeit.Timer(
            stmt="IndexedFasta.benchmarkRandomFragment(fasta=fasta, size=%i)" %
            (args.benchmark_fragment_size),
            setup="from cgat import IndexedFasta\n"
            "fasta=IndexedFasta.IndexedFasta('%s')" % (unknown[0]))

        t = timer.timeit(number=args.benchmark_num_iterations)
        args.stdout.write("iter\tsize\ttime\n")
        args.stdout.write("%i\t%i\t%i\n" % (
            args.benchmark_num_iterations,
            args.benchmark_fragment_size, t))

    elif args.verify:
        fasta1 = IndexedFasta.IndexedFasta(unknown[0])
        fasta2 = IndexedFasta.IndexedFasta(args.verify)
        nerrors1 = IndexedFasta.verify(fasta1, fasta2,
                                       args.verify_num_iterations,
                                       args.verify_fragment_size,
                                       stdout=args.stdout)
        args.stdout.write("errors=%i\n" % (nerrors1))
        nerrors2 = IndexedFasta.verify(fasta2, fasta1,
                                       args.verify_num_iterations,
                                       args.verify_fragment_size,
                                       stdout=args.stdout)
        args.stdout.write("errors=%i\n" % (nerrors2))
    elif args.compress_index:
        fasta = IndexedFasta.IndexedFasta(unknown[0])
        fasta.compressIndex()
    else:
        if args.loglevel >= 1:
            args.stdlog.write("# creating database %s\n" % unknown[0])
            args.stdlog.write("# indexing the following files: \n# %s\n" %
                              (" \n# ".join(unknown[1:])))
            args.stdlog.flush()

            if synonyms:
                args.stdlog.write("# Applying the following synonyms:\n")
                for k, v in list(synonyms.items()):
                    args.stdlog.write("# %s=%s\n" % (k, ",".join(v)))
                args.stdlog.flush()
        if len(unknown) < 2:
            print(globals()["__doc__"])
            sys.exit(1)

        iterator = IndexedFasta.MultipleFastaIterator(
            unknown[1:],
            regex_identifier=args.regex_identifier,
            format=args.file_format)

        IndexedFasta.createDatabase(
            unknown[0],
            iterator,
            synonyms=synonyms,
            random_access_points=args.random_access_points,
            compression=args.compression,
            clean_sequence=args.clean_sequence,
            allow_duplicates=args.allow_duplicates,
            translator=args.translator,
            force=args.force)

    E.stop()


if __name__ == "__main__":
    sys.exit(main())
