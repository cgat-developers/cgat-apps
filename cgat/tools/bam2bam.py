'''bam2bam.py - modify bam files
=============================

Purpose
-------

This script reads a :term:`bam` formatted file from stdin, performs an
action (see methods below) then outputs a modified :term:`bam`
formatted file on stdout.

.. note::
   You need to redirect logging information to a file (via -L) or turn it off
   via -v 0 in order to get a valid sam/bam file.

Documentation
-------------

The script implements the following methods:

``set-nh``

   set the NH flag. Some tools (bowtie, bwa) do not set the NH flag.
   If set, this option will set the NH flag (for mapped reads).
   This option requires the bam/sam file to be sorted by read name.

``unset-unmapped_mapq``

   some tools set the mapping quality of unmapped reads. This
   causes a violation in the Picard tools.

``filter``

   remove alignments based on a variety of flags. The filtering method
   is determined by the ``--filter-method`` option. These may be
   ``unique``, ``non-unique``, ``mapped``, ``NM`` or ``CM``.  If
   ``unique`` is set, only uniquely mapping reads will be output. If
   ``non-unique`` is set then only multi-mapping reads will be
   output. This method first checks for the NH flag - if set, a unique
   match should have at most NH=1 hits.  If not set, the method checks
   for BWA flags. Currently it checks if X0 is set (X0=Number of best
   hits found by BWA).  If ``mapped`` is given, unmapped reads will be
   removed. If ``NM`` or ``CM`` is set, the alignment of reads in two
   sam files (input and reference) is compared and only reads with a
   lower number of mismatches in the input compared to the reference
   sam file will be kept. If ``CM`` is set, the colourspace mismatch
   tag (for ABI Solid reads) will be used to count differences to the
   reference sam file. By default, the ``NM`` (number of mismatches)
   tag is used. The tag that is used needs to present in both input
   sam file and the reference sam file. If ``unique`` is given this
   wil NOT remove any unmapped reads.  This can be achieved by
   providing the ``filter`` option twice, once each with ``mapped``
   and ``unique``.

   .. note::

      The filter methods can't currently combined with any of
      the other methods - this is work in progress.

``strip-sequence``

   remove the sequence from all reads in a bam-file. Note that
   stripping the sequence will also remove the quality scores.
   Stripping is not reversible if the read names are not unique.

``strip-quality``

   remove the quality scores from all reads in a bam-file.
   Stripping is not reversible if the read names are not unique.

``set-sequence``

   set the sequence and quality scores in the bam file to some dummy
   values ('A' for sequence, 'F' for quality which is a valid score in
   most fastq encodings. Necessary for some tools that can not work
   with bam-files without sequence.

``unstrip``

   add sequence and quality scores back to a bam file. Requires a
   :term:`fastq` formatted file with the sequences and quality scores
   to insert.

``unset-unmapped-mapq``

   sets the mapping quality of unmapped reads to 0.

``keep-first-base``

   keep only the first base of reads so that read counting tools will
   only consider the first base in the counts

``downsample-single``

   generates a downsampled :term:`bam` file by randomly subsampling
   reads from a single ended :term:`bam` file. The downsmpling
   retains multimapping reads. The use of this requires downsampling
   parameter to be set and optionally randomseed.

``downsample-paired``

   generates a downsampled :term:`bam` file by randomly subsampling
   reads from a paired ended :term:`bam` file. The downsampling
   retains multimapping reads. The use of this requires downsampling
   parameter to be set and optionally randomseed.

``add-sequence-error``

   add a certain amount of random error to read sequences. This method
   picks a certain proportion of positions within a read's sequence
   and alters the nucleotide to a randomly chosen alternative. The
   model is naive and applies uniform probabilities for positions and
   nucleotides. The method does not update base qualities, the
   alignment and the NM flag. As a result, error rates that are
   computed via the NM flag will be unaffected. The error rate is set
   by --error-rate.

By default, the script works from stdin and outputs to stdout.

Usage
-----

For example::

   cgat bam2bam --method=filter --filter-method=mapped < in.bam > out.bam

will remove all unmapped reads from the bam-file.

Example for running downsample::

cgat bam2bam --method=downsample-paired --downsample=30000
--randomseed=1 -L out.log < Paired.bam > out.bam

Type::

   cgat bam2bam --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import tempfile
import shutil
import random
import pysam
import cgatcore.experiment as E
import cgatcore.iotools as iotools
import itertools
import math

from cgat.BamTools.bamtools import bam2bam_filter_bam, SetNH


class SubsetBam(object):

    ''' base class for performing downsampling on single and
    paired bam file

    A dictionary of the read names is made and then an
    array matching the numbers of reads in the dictionary will be
    created. This array contains 1s that match the number of reads
    to be downsampled to and 0s to fill out the rest of the array.

    The read is kept if the value in the dictionary matches 1. The
    read is yielded and iterated over to produce the output bam.

    The script will handle multimapping reads.
    '''

    def __init__(self, infile, downsample, paired_end=None,
                 single_end=None, random_seed=None):

        self.pysam_in1, self.pysam_in2 = itertools.tee(infile)
        self.downsample = downsample
        self.paired_end = paired_end
        self.single_end = single_end
        self.random_seed = random_seed

    def list_of_reads(self, paired=None):

        '''
        This will create a dictionary of uniqe reads in the bam
        '''
        read_list = []

        for read in self.pysam_in1:
            if paired is True:
                if read.is_proper_pair:
                    read_list.append(read.qname)
            else:
                read_list.append(read.qname)

        return sorted(set(read_list))

    def downsample_paired(self):

        '''
        This function will downsample a paired bam file.
        It will retain multimapping reads if they have not been
        pre-filtered
        '''
        if self.random_seed is not None:
            random.seed(self.random_seed)

        collect_list = self.list_of_reads(paired=True)
        read_list = random.sample(collect_list, self.downsample)

        if self.downsample == len(collect_list):
            E.warn('''The downsample reads is equal to the
            number of unique reads''')
            for read in self.pysam_in2:
                yield read

        # yield read if it is in read_list, if the read is multimapped
        # then all multimapping reads will be yielded
        else:
            E.warn('''Multimaping reads have been detected and these will
            be output to the final bam file''')
            for read in self.pysam_in2:
                if read.qname in read_list:
                    yield read

    def downsample_single(self):

        '''
        This function will downsample a single bam file.
        It will retain multimapping reads if not pre-filtered
        '''
        if self.random_seed is not None:
            random.seed(self.random_seed)

        collect_list = self.list_of_reads(paired=False)
        read_list = random.sample(collect_list, self.downsample)

        if self.downsample == len(collect_list):
            E.warn('''The downsample reads is equal to the
            number of unique reads''')
            for read in self.pysam_in2:
                yield read

        # yield read if it is in read_list, if the reads is multimapped
        # then all multimapping reads will be yielded
        else:
            E.warn('''Multimaping reads have been detected and these will
            be output to the final bam file''')
            for read in self.pysam_in2:
                if read.qname in read_list:
                    yield read


def process_bam(infile, outfile, options):

    if "filter" in options.methods:
        if "remove-list" in options.filter_methods or "keep-list" in options.filter_methods:

            it = infile.fetch(until_eof=True)
            c = E.Counter()
            if "remove-list" in options.filter_methods:
                for read in it:
                    c.input += 1
                    if read.query_name in filter_query_names:
                        c.skipped += 1
                        continue
                    outfile.write(read)
                    c.output += 1
            elif "keep-list" in options.filter_methods:
                for read in it:
                    c.input += 1
                    if read.query_name not in filter_query_names:
                        c.skipped += 1
                        continue
                    outfile.write(read)
                    c.output += 1

            E.info("category\tcounts\n%s\n" % c.asTable())
        else:
            remove_mismatches, colour_mismatches = False, False

            if "NM" in options.filter_methods:
                remove_mismatches = True

            elif "CM" in options.filter_methods:
                remove_mismatches = True
                colour_mismatches = True

            if "min-length" in options.filter_methods and options.minimum_read_length == 0:
                raise ValueError("please specify --minimum-read-length when using "
                                 "--filter-method=min-read-length")

            if "min-average-base-quality" in options.filter_methods and options.minimum_average_base_quality == 0:
                raise ValueError("please specify --min-average-base-quality when "
                                 "using --filter-method=min-average-base-quality")

            if remove_mismatches:
                if not options.reference_bam:
                    raise ValueError(
                        "requiring reference bam file for removing by "
                        "mismatches")

                pysam_ref = pysam.AlignmentFile(options.reference_bam, "rb")
            else:
                pysam_ref = None

            # filter and flags are the opposite way around
            c = bam2bam_filter_bam(
                infile, outfile, pysam_ref,
                remove_nonunique="unique" in options.filter_methods,
                remove_unique="non-unique" in options.filter_methods,
                remove_contigs=None,
                remove_unmapped="mapped" in options.filter_methods,
                remove_mismatches=remove_mismatches,
                filter_error_rate=options.error_rate,
                colour_mismatches=colour_mismatches,
                minimum_read_length=options.minimum_read_length,
                minimum_average_base_quality=options.minimum_average_base_quality)

            options.stdlog.write("category\tcounts\n%s\n" % c.asTable())
    else:

        # set up the modifying iterators
        it = infile.fetch(until_eof=True)

        def nop(x):
            return None
        # function to check if processing should start
        pre_check_f = nop

        if "unset-unmapped-mapq" in options.methods:
            def unset_unmapped_mapq(i):
                for read in i:
                    if read.is_unmapped:
                        read.mapq = 0
                    yield read
            it = unset_unmapped_mapq(it)

        if "set-sequence" in options.methods:
            def set_sequence(i):
                for read in i:
                    # can't get at length of unmapped reads
                    if read.is_unmapped:
                        read.seq = "A"
                        read.qual = "F"
                    else:
                        read.seq = "A" * read.inferred_length
                        read.qual = "F" * read.inferred_length

                    yield read
            it = set_sequence(it)

        if "strip-sequence" in options.methods or "strip-quality" in \
           options.methods:
            def strip_sequence(i):
                for read in i:
                    read.seq = None
                    yield read

            def check_sequence(reads):
                if reads[0].seq is None:
                    return 'no sequence present'
                return None

            def strip_quality(i):
                for read in i:
                    read.qual = None
                    yield read

            def check_quality(reads):
                if reads[0].qual is None:
                    return 'no quality information present'
                return None

            def strip_match(i):
                for read in i:
                    try:
                        nm = read.opt('NM')
                    except KeyError:
                        nm = 1
                    if nm == 0:
                        read.seq = None
                    yield read

            if options.strip_method == "all":
                if "strip-sequence" in options.methods:
                    it = strip_sequence(it)
                    pre_check_f = check_sequence
                elif "strip-quality" in options.methods:
                    it = strip_quality(it)
                    pre_check_f = check_quality
            elif options.strip_method == "match":
                it = strip_match(it)

        if "unstrip" in options.methods:
            def buildReadDictionary(filename):
                if not os.path.exists(filename):
                    raise OSError("file not found: %s" % filename)
                fastqfile = pysam.FastxFile(filename)
                fastq2sequence = {}
                for x in fastqfile:
                    if x.name in fastq2sequence:
                        raise ValueError(
                            "read %s duplicate - can not unstrip" % x.name)

                    fastq2sequence[x.name] = (x.sequence, x.quality)
                return fastq2sequence

            if not options.fastq_pair1:
                raise ValueError(
                    "please supply fastq file(s) for unstripping")
            fastq2sequence1 = buildReadDictionary(options.fastq_pair1)
            if options.fastq_pair2:
                fastq2sequence2 = buildReadDictionary(options.fastq_pair2)

            def unstrip_unpaired(i):
                for read in i:
                    read.seq, read.qual = fastq2sequence1[read.qname]
                    yield read

            def unstrip_pair(i):
                for read in i:
                    if read.is_read1:
                        read.seq, read.qual = fastq2sequence1[read.qname]
                    else:
                        read.seq, read.qual = fastq2sequence2[read.qname]
                    yield read

            if options.fastq_pair2:
                it = unstrip_pair(it)
            else:
                it = unstrip_unpaired(it)

        if "set-nh" in options.methods:
            it = SetNH(it)

        # keep first base of reads by changing the cigarstring to
        # '1M' and, in reads mapping to the reverse strand,
        # changes the pos to aend - 1
        # Needs to be refactored to make it more general
        # (last base, midpoint, ..)
        if "keep_first_base" in options.methods:
            def keep_first_base(i):
                for read in i:
                    if read.is_reverse:
                        read.pos = read.aend - 1
                        read.cigarstring = '1M'
                    elif not read.is_unmapped:
                        read.cigarstring = '1M'
                    yield read
            it = keep_first_base(it)

        # read first read and check if processing should continue
        # only possible when not working from stdin
        # Refactoring: use cache to also do a pre-check for
        # stdin input.
        if not infile.is_stream:
            # get first read for checking pre-conditions
            first_reads = list(infile.head(1))

            msg = pre_check_f(first_reads)
            if msg is not None:
                if options.force:
                    E.warn('proccessing continues, though: %s' % msg)
                else:
                    E.warn('processing not started: %s' % msg)
                    return

        if "downsample-single" in options.methods:

            if not options.downsample:
                raise ValueError("Please provide downsample size")

            else:
                down = SubsetBam(infile=it,
                                 downsample=options.downsample,
                                 paired_end=None,
                                 single_end=True,
                                 random_seed=options.random_seed)
                it = down.downsample_single()

        if "downsample-paired" in options.methods:

            if not options.downsample:
                raise ValueError("Please provide downsample size")

            else:
                down = SubsetBam(infile=it,
                                 downsample=options.downsample,
                                 paired_end=True,
                                 single_end=None,
                                 random_seed=options.random_seed)
                it = down.downsample_paired()

        if "add-sequence-error" in options.methods:
            def add_sequence_error(i):
                error_rate = options.error_rate
                map_nuc2var = {"A": "CGT",
                               "C": "AGT",
                               "G": "ACT",
                               "T": "ACG"}
                for read in i:
                    sequence = list(read.query_sequence)
                    quals = read.query_qualities
                    npos = int(math.floor(len(sequence) * error_rate))
                    positions = random.sample(range(len(sequence)), npos)
                    for pos in positions:
                        try:
                            alt = map_nuc2var[sequence[pos]]
                        except KeyError:
                            continue
                        sequence[pos] = alt[random.randint(0, len(alt) - 1)]

                    read.query_sequence = "".join(sequence)
                    read.query_qualities = quals
                    yield read

            it = add_sequence_error(it)

        # continue processing till end
        for read in it:
            outfile.write(read)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("-m", "--methods", dest="methods", type=str,
                        action="append",
                        choices=("filter",
                                 "keep-first-base",
                                 "set-nh",
                                 "set-sequence",
                                 "strip-sequence",
                                 "strip-quality",
                                 "unstrip",
                                 "unset-unmapped-mapq",
                                 "downsample-single",
                                 "downsample-paired",
                                 "add-sequence-error"),
                        help="methods to apply ")

    parser.add_argument("--strip-method", dest="strip_method", type=str,
                        choices=("all", "match"),
                        help="define which sequences/qualities to strip. "
                        "match means that stripping only applies to entries "
                        "without mismatches (requires NM tag to be present). "
                        )

    parser.add_argument("--filter-method", dest="filter_methods",
                        action="append", type=str,
                        choices=('NM', 'CM',
                                 "mapped", "unique", "non-unique",
                                 "remove-list",
                                 "keep-list",
                                 "error-rate",
                                 "min-read-length",
                                 "min-average-base-quality"),
                        help="filter method to apply to remove alignments "
                        "from a bam file. Multiple methods can be supplied "
                        )

    parser.add_argument("--reference-bam-file", dest="reference_bam",
                        type=str,
                        help="bam-file to filter with ")

    parser.add_argument("--force-output", dest="force", action="store_true",
                        help="force processing. Some methods such "
                        "as strip/unstrip will stop processing if "
                        "they think it not necessary "
                        )

    parser.add_argument("--output-sam", dest="output_sam", action="store_true",
                        help="output in sam format ")

    parser.add_argument(
        "--first-fastq-file", "-1", dest="fastq_pair1", type=str,
        help="fastq file with read information for first "
        "in pair or unpaired. Used for unstripping sequence "
        "and quality scores ")

    parser.add_argument(
        "--second-fastq-file", "-2", dest="fastq_pair2", type=str,
        help="fastq file with read information for second "
        "in pair. Used for unstripping sequence "
        "and quality scores  ")

    parser.add_argument(
        "--downsample", dest="downsample",
        type=int,
        help="Number of reads to downsample to")

    parser.add_argument(
        "--filename-read-list", dest="filename_read_list",
        type=str,
        help="Filename with list of reads to filter if 'keep-list' or 'remove-list' "
        "filter method is chosen ")

    parser.add_argument(
        "--error-rate", dest="error_rate",
        type=float,
        help="error rate to use as filter. Reads with an error rate "
        "higher than the threshold will be removed ")

    parser.add_argument(
        "--minimum-read-length", dest="minimum_read_length",
        type=int,
        help="minimum read length when filtering ")

    parser.add_argument(
        "--minimum-average-base-quality", dest="minimum_average_base_quality",
        type=float,
        help="minimum average base quality when filtering ")

    parser.set_defaults(
        methods=[],
        output_sam=False,
        reference_bam=None,
        filter_methods=[],
        strip_method="all",
        force=False,
        fastq_pair1=None,
        fastq_pair2=None,
        downsample=None,
        random_seed=None,
        filename_read_list=None,
        error_rate=None,
        minimum_read_length=0,
        minimum_average_base_quality=0,
    )

    # add common options (-h/--help, ...) and parse command line
    (args, unknown) = E.start(parser, argv=argv, unknowns=True)

    if args.stdin != sys.stdin:
        bamfile = args.stdin.name
    elif unknown:
        bamfile = unknown[0]
        if len(unknown) > 1:
            raise ValueError("multiple bam files provided in arguments")
    else:
        bamfile = "-"

    if "remove-list" in args.filter_methods or "keep-list" in args.filter_methods:
        if "remove-list" in args.filter_methods and "keep-list" in args.filter_methods:
            raise ValueError("it is not possible to specify remove-list and keep-list")

        with iotools.open_file(args.filename_read_list) as inf:
            filter_query_names = set([x.strip() for x in inf.readlines() if not x.startswith("#")])
        E.info("read query_sequence filter list with {} read names".format(len(filter_query_names)))

    if "error-rate" in args.filter_methods and not args.error_rate:
        raise ValueError("filtering by error-rate requires --error-rate to be set")

    if "add-sequence-error" in args.methods and not args.error_rate:
        raise ValueError("--add-error-rate requires --error-rate to be set")

    E.info('processing %s' % bamfile)
    if bamfile != "-" and iotools.is_empty(bamfile):
        E.warn('ignoring empty file %s' % bamfile)
        E.stop()
        return

    if args.stdout != sys.stdout:
        output_bamfile = args.stdout.name
    else:
        output_bamfile = "-"
        if args.stdlog == sys.stdout:
            raise ValueError("redirect log-stream to file (--log) if outputting to stdout")

    if args.output_sam:
        output_mode = "wh"
    else:
        output_mode = "wb"

    # reading bam from stdin does not work with only the "r" tag
    with pysam.AlignmentFile(bamfile, "rb") as pysam_in:
        with pysam.AlignmentFile(output_bamfile, output_mode,
                                 template=pysam_in) as pysam_out:
            process_bam(pysam_in, pysam_out, args)

    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
