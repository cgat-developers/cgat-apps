'''gff2gff.py - manipulate gff files
=================================

:Tags: Genomics Intervals GFF Manipulation

Note: The script can parse AGP files but this functionality
is scheduled for deprecation in the next release.

Purpose
-------

This scripts reads a :term:`gff` formatted file, applies a
transformation and outputs the new intervals in :term:`gff` format.
The type of transformation chosen is given through the `--method``
option. Below is a list of available transformations:

``complement-groups``

    output the complenent intervals for the features in the file, for
    example to output introns from exons. The option ``--group-field``
    sets field/attribute to group by, e.g gene_id, transcript_id, feature,
    source.

``combine-groups``

    combine all features in a group into a single interval.  The
    option ``--group-field`` sets field/attribute to group by, see
    alse ``complement-groups``.

``to-forward-coordinates``

    translate all features forward coordinates.

``to-forward-strand``

    convert to forward strand

``add-upstream-flank/add-downstream-flank/add-flank``

   add an upstream/downstream flanking segment to first/last exon of a group.
   The amount added is given through the options ``--extension-upstream`` and
   ``--extension-downstream``. If ``--flank-method`` is ``extend``, the
   first/last exon will be extended, otherwise a new feature will be added.

``crop``

   crop features according to features in a separate gff file.
   If a feature falls in the middle of another, two entries will be
   output.""" )

``crop-unique``

   remove non-unique features from gff file.

``merge-features``

   merge consecutive features.

``join-features``

   group consecutive features.

``filter-range``
   extract features overlapping a chromosomal range. The range can be
   set by the ``--filter-range`` option.

``sanitize``
   reconcile chromosome names between ENSEMBL/UCSC or with an indexed
   genomic fasta file (see :doc:`index_fasta`). Raises an exception if
   an unknown contig is found (unless ``--skip-missing`` is set). The
   method to sanitize is specified by ``--sanitize-method``.The
   method to sanitize is specified by ``--sanitize-method``. Options for
   ```--sanitize-method``` include "ucsc", "ensembl", "genome".
   A pattern of contigs to remove can be given in the option
   ``--contig-pattern``.
   If ``--sanitize-method`` is set to ``ucsc`` or ``ensembl``, the option
   ``--assembly-report`` is required to allow for accurate mapping
   between UCSC and Ensembl. If not found in the assembly report the
   contig names are forced into the desired convention, either by removing
   or prepending ``chr``, this is useful for :term:`gff` files with custom
   contigs. The Assembly Report can be found on the NCBI assembly page
   under the link "Download the full sequence report".
   If ``--sanitize-method`` is set to ``genome``, the genome file has to be
   provided via the option ``--genome-file`` or ``--contigs-tsv-file``

``skip-missing``

   skip entries on missing contigs. This prevents exception from being raised

``filename-agp``

    agp file to map coordinates from contigs to scaffolds

``rename-chr``

    Renames chromosome names. Source and target names are supplied as a file
    with two columns. Examples are available at:
    https://github.com/dpryan79/ChromosomeMappings
    Note that unmapped chromosomes are dropped from the output file.


Usage
-----

For many downstream applications it is helpful to make sure
that a :term:`gff` formatted file contains only features on
placed chromosomes.

As an example, to sanitise hg38 chromosome names and remove
chromosome matching the regular expression patterns
"ChrUn", "_alt" or "_random", use the following:

   cat in.gff
   | gff2gff.py --method=sanitize --sanitize-method=ucsc
                --assembly-report=/path/to/file --skip-missing
   | gff2gff.py --remove-contigs="chrUn,_random,_alt" > gff.out

The "--skip-missing" option prevents an exception being
raised if entries are found on missing chromosomes

Another example, to rename UCSC chromosomes to ENSEMBL.

    cat ucsc.gff
    | gff2gff.py --method=rename-chr
                 --rename-chr-file=ucsc2ensembl.txt > ensembl.gff

Type::

   cgat gff2gff --help

for command line help.

Command line options
--------------------

'''

import sys
import re
import collections
import numpy
import quicksect
import pandas as pd
import cgatcore.experiment as E
import cgatcore.iotools as iotools
import cgat.GTF as GTF
import cgat.AGP as AGP
import cgat.Genomics as Genomics
import cgat.IndexedFasta as IndexedFasta
import cgat.Intervals as Intervals
import csv


def combineGFF(gffs,
               min_distance,
               max_distance,
               min_features,
               max_features,
               merge=True,
               output_format="%06i"):
    """join intervals in gff file.

    Note: strandedness is ignored
    """

    E.info("joining features: min distance=%i, max_distance=%i, "
           "at least %i and at most %i features." %
           (min_distance, max_distance, min_features, max_features))

    def iterate_chunks(gffs):

        last = next(gffs)
        to_join = [last]

        for gff in gffs:
            d = gff.start - last.end
            if gff.contig == last.contig:
                assert gff.start >= last.start, "input file should be sorted by contig and position: d=%i:\n%s\n%s\n" % (
                    d, last, gff)

            if gff.contig != last.contig or \
                    (max_distance and d > max_distance) or \
                    (min_distance and d < min_distance) or \
                    (max_features and len(to_join) >= max_features):

                if min_features or len(to_join) >= min_features:
                    yield to_join
                to_join = []

            last = gff
            to_join.append(gff)

        if len(to_join) >= min_features:
            yield to_join
        raise StopIteration

    id = 1
    ninput, noutput, nfeatures = 0, 0, 0

    if merge:
        for to_join in iterate_chunks(gffs):

            ninput += 1
            y = GTF.Entry()
            t = output_format % id
            y.fromGTF(to_join[0], t, t)
            y.start = to_join[0].start
            y.end = to_join[-1].end

            yield(y)
            nfeatures += 1

            noutput += 1
            id += 1
    else:

        for to_join in iterate_chunks(gffs):

            ninput += 1
            for x in to_join:
                y = GTF.Entry()
                t = output_format % id
                y.fromGTF(x, t, t)
                yield(y)
                nfeatures += 1

            noutput += 1
            id += 1

    E.info("ninput=%i, noutput=%i, nfeatures=%i" %
           (ninput, noutput, nfeatures))


def cropGFFUnique(gffs, ignore_strand=True):
    """crop intervals in gff file.

    Only unique regions are kept. This method ignores the feature field.

    If ignore_strand is set, strand information for cropping
    is ignored.
    """

    # read regions to crop with and convert intervals to intersectors
    E.info("reading gff for cropping: started.")
    gffs = list(gffs)
    if len(gffs) == 0:
        return

    def _cmp_without_strand(this, last):
        return this.contig != last.contig

    def _cmp_with_strand(this, last):
        return this.contig != last.contig or this.strand != last.strand

    if ignore_strand:
        gffs.sort(key=lambda x: (x.contig, x.start))
        comp = _cmp_without_strand
    else:
        gffs.sort(key=lambda x: (x.contig, x.strand, x.start))
        comp = _cmp_with_strand

    last = gffs[0]
    c = E.Counter()

    c.input = len(gffs)

    for this in gffs[1:]:
        if comp(this, last) or last.end <= this.start:
            # no overlap
            if last.start < last.end:
                c.output += 1
                yield(last)
            last = this

        elif this.end <= last.start:
            # this ends before last (happens in multiway overlaps)
            # nothing happens
            pass

        elif last.end <= this.end:
            # last ends before this
            l = last.end
            last.end = this.start
            this.start = l
            if last.start < last.end:
                E.info("overlap")
                E.info(str(this))
                E.info(str(last))
                c.overlaps += 1
                c.output += 1
                yield(last)
            last = this

        elif last.end > this.end:
            # this within last - split last
            l = last.end
            last.end = this.start
            if last.start < last.end:
                c.output += 1
                c.splits += 1
                yield(last)
            last.start = this.end
            last.end = l

    if last.start < last.end:
        c.output += 1
        yield(last)

    E.info("cropping finished: %s" % str(c))


def cropGFF(gffs, filename_gff):
    """crop intervals in gff file."""

    # read regions to crop with and convert intervals to intersectors
    E.info("reading gff for cropping: started.")

    other_gffs = GTF.iterator(iotools.open_file(filename_gff, "r"))

    cropper = GTF.readAsIntervals(other_gffs)

    ntotal = 0
    for contig in list(cropper.keys()):
        intersector = quicksect.IntervalTree()
        for start, end in cropper[contig]:
            intersector.add(start, end)
            ntotal += 1
        cropper[contig] = intersector

    E.info("reading gff for cropping: finished.")
    E.info("reading gff for cropping: %i contigs with %i intervals." %
           (len(cropper), ntotal))

    ninput, noutput, ncropped, ndeleted = 0, 0, 0, 0

    # do the actual cropping
    for gff in gffs:

        ninput += 1

        if gff.contig in cropper:

            start, end = gff.start, gff.end
            overlaps = cropper[gff.contig].find(quicksect.Interval(start, end))

            if overlaps:

                l = end - start
                a = numpy.ones(l)
                for i in overlaps:
                    s = max(0, i.start - start)
                    e = min(l, i.end - start)
                    a[s:e] = 0

                segments = Intervals.fromArray(a)
                if len(segments) == 0:
                    ndeleted += 1
                else:
                    ncropped += 1

                for s, e in segments:
                    gff.start, gff.end = s + start, e + start
                    noutput += 1
                    yield(gff)

                continue

        noutput += 1

        yield(gff)

    E.info("ninput=%i, noutput=%i, ncropped=%i, ndeleted=%i" %
           (ninput, noutput, ncropped, ndeleted))


def renameChromosomes(gffs, chr_map):

    ninput, noutput, nskipped = 0, 0, 0

    for gff in gffs:

        ninput += 1

        if gff.contig in chr_map.keys():
            gff.contig = chr_map[gff.contig]
        else:
            nskipped += 1
            continue

        noutput += 1
        yield gff

    E.info("ninput = %i, noutput=%i, nskipped=%i" %
           (ninput, noutput, nskipped))


def main(argv=None):

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("-m", "--method", dest="method", type=str,
                        choices=("add-flank",
                                 "add-upstream-flank",
                                 "add-downstream-flank",
                                 "crop",
                                 "crop-unique",
                                 "complement-groups",
                                 "combine-groups",
                                 "filter-range",
                                 "join-features",
                                 "merge-features",
                                 "sanitize",
                                 "to-forward-coordinates",
                                 "to-forward-strand",
                                 "rename-chr"),
                        help="method to apply ")

    parser.add_argument(
        "--ignore-strand", dest="ignore_strand",
        help="ignore strand information.", action="store_true")

    parser.add_argument("--is-gtf", dest="is_gtf", action="store_true",
                        help="input will be treated as gtf.")

    parser.add_argument(
        "-c", "--contigs-tsv-file", dest="input_filename_contigs",
        type=str,
        help="filename with contig lengths.")

    parser.add_argument(
        "--agp-file", dest="input_filename_agp", type=str,
        help="agp file to map coordinates from contigs to scaffolds.")

    parser.add_argument(
        "-g", "--genome-file", dest="genome_file", type=str,
        help="filename with genome.")

    parser.add_argument(
        "--crop-gff-file", dest="filename_crop_gff", type=str,
        help="GFF/GTF file to crop against.")

    parser.add_argument(
        "--group-field", dest="group_field", type=str,
        help="""gff field/attribute to group by such as gene_id, "
        "transcript_id, ... .""")

    parser.add_argument(
        "--filter-range", dest="filter_range", type=str,
        help="extract all elements overlapping a range. A range is "
        "specified by eithor 'contig:from..to', 'contig:+:from..to', "
        "or 'from,to' .")

    parser.add_argument(
        "--sanitize-method", dest="sanitize_method", type=str,
        choices=("ucsc", "ensembl", "genome"),
        help="method to use for sanitizing chromosome names. "
        ".")

    parser.add_argument(
        "--flank-method", dest="flank_method", type=str,
        choices=("add", "extend"),
        help="method to use for adding flanks. ``extend`` will "
        "extend existing features, while ``add`` will add new features. "
        ".")

    parser.add_argument(
        "--skip-missing", dest="skip_missing", action="store_true",
        help="skip entries on missing contigs. Otherwise an "
        "exception is raised .")

    parser.add_argument(
        "--contig-pattern", dest="contig_pattern", type=str,
        help="a comma separated list of regular expressions specifying "
        "contigs to be removed when running method sanitize .")

    parser.add_argument(
        "--assembly-report", dest="assembly_report", type=str,
        help="path to assembly report file which allows mapping of "
        "ensembl to ucsc contigs when running method sanitize .")

    parser.add_argument(
        "--assembly-report-hasids",
        dest="assembly_report_hasIDs", type=int,
        help="path to assembly report file which allows mapping of "
        "ensembl to ucsc contigs when running method sanitize .")

    parser.add_argument(
        "--assembly-report-ucsccol", dest="assembly_report_ucsccol",
        type=int,
        help="column in the assembly report containing ucsc contig ids"
        ".")

    parser.add_argument(
        "--assembly-report-ensemblcol", dest="assembly_report_ensemblcol",
        type=int,
        help="column in the assembly report containing ensembl contig ids")

    parser.add_argument(
        "--assembly-extras", dest="assembly_extras",
        type=str,
        help="additional mismatches between gtf and fasta to fix when"
        "sanitizing the genome .")

    parser.add_argument(
        "--extension-upstream", dest="extension_upstream", type=float,
        help="extension for upstream end .")

    parser.add_argument(
        "--extension-downstream", dest="extension_downstream", type=float,
        help="extension for downstream end .")

    parser.add_argument(
        "--min-distance", dest="min_distance", type=int,
        help="minimum distance of features to merge/join .")

    parser.add_argument(
        "--max-distance", dest="max_distance", type=int,
        help="maximum distance of features to merge/join .")

    parser.add_argument(
        "--min-features", dest="min_features", type=int,
        help="minimum number of features to merge/join .")

    parser.add_argument(
        "--max-features", dest="max_features", type=int,
        help="maximum number of features to merge/join .")

    parser.add_argument(
        "--rename-chr-file", dest="rename_chr_file", type=str,
        help="mapping table between old and new chromosome names."
        "TAB separated 2-column file.")

    parser.set_defaults(
        input_filename_contigs=False,
        filename_crop_gff=None,
        input_filename_agp=False,
        genome_file=None,
        rename_chr_file=None,
        add_up_flank=None,
        add_down_flank=None,
        complement_groups=False,
        crop=None,
        crop_unique=False,
        ignore_strand=False,
        filter_range=None,
        min_distance=0,
        max_distance=0,
        min_features=1,
        max_features=0,
        extension_upstream=1000,
        extension_downstream=1000,
        sanitize_method="ucsc",
        flank_method="add",
        output_format="%06i",
        skip_missing=False,
        is_gtf=False,
        group_field=None,
        contig_pattern=None,
        assembly_report=None,
        assembly_report_hasIDs=1,
        assembly_report_ensemblcol=4,
        assembly_report_ucsccol=9,
        assembly_extras=None
    )

    (args) = E.start(parser, argv=argv)

    contigs = None
    genome_fasta = None
    chr_map = None

    if args.input_filename_contigs:
        contigs = Genomics.readContigSizes(
            iotools.open_file(args.input_filename_contigs, "r"))

    if args.genome_file:
        genome_fasta = IndexedFasta.IndexedFasta(args.genome_file)
        contigs = genome_fasta.getContigSizes()

    if args.rename_chr_file:
        chr_map = {}
        with open(args.rename_chr_file, 'r') as filein:
            reader = csv.reader(filein, delimiter='\t')
            for row in reader:
                if len(row) != 2:
                    raise ValueError(
                        "Mapping table must have exactly two columns")
                chr_map[row[0]] = row[1]
        if not len(chr_map.keys()) > 0:
            raise ValueError("Empty mapping dictionnary")

    if args.assembly_report:
        df = pd.read_csv(args.assembly_report, comment="#",
                         header=None, sep="\t")
        # fixes naming inconsistency in assembly report: ensembl chromosome
        # contigs found in columnn 0, ensembl unassigned contigs found in
        # column 4.
        if args.assembly_report_hasIDs == 1:
            ucsccol = args.assembly_report_ucsccol
            ensemblcol = args.assembly_report_ensemblcol
            df.loc[df[1] == "assembled-molecule", ensemblcol] = df.loc[
                df[1] == "assembled-molecule", 0]
            if args.sanitize_method == "ucsc":
                assembly_dict = df.set_index(ensemblcol)[ucsccol].to_dict()
            elif args.sanitize_method == "ensembl":
                assembly_dict = df.set_index(ucsccol)[ensemblcol].to_dict()
            else:
                raise ValueError(''' When using assembly report,
                please specify sanitize method as either
                "ucsc" or "ensembl" to specify direction of conversion
                ''')
        else:
            assembly_dict = {}
        if args.assembly_extras is not None:
            assembly_extras = args.assembly_extras.split(",")
            for item in assembly_extras:
                item = item.split("-")
                assembly_dict[item[0]] = item[1]

    if args.method in ("forward_coordinates", "forward_strand",
                       "add-flank", "add-upstream-flank",
                       "add-downstream-flank") \
       and not contigs:
        raise ValueError("inverting coordinates requires genome file")

    if args.input_filename_agp:
        agp = AGP.AGP()
        agp.readFromFile(iotools.open_file(args.input_filename_agp, "r"))
    else:
        agp = None

    gffs = GTF.iterator(args.stdin)

    if args.method in ("add-upstream-flank",
                       "add-downstream-flank",
                       "add-flank"):

        add_upstream_flank = "add-upstream-flank" == args.method
        add_downstream_flank = "add-downstream-flank" == args.method
        if args.method == "add-flank":
            add_upstream_flank = add_downstream_flank = True

        upstream_flank = int(args.extension_upstream)
        downstream_flank = int(args.extension_downstream)
        extend_flank = args.flank_method == "extend"

        if args.is_gtf:
            iterator = GTF.flat_gene_iterator(gffs)
        else:
            iterator = GTF.joined_iterator(gffs, args.group_field)

        for chunk in iterator:
            is_positive = Genomics.IsPositiveStrand(chunk[0].strand)
            chunk.sort(key=lambda x: (x.contig, x.start))
            lcontig = contigs[chunk[0].contig]

            if extend_flank:
                if add_upstream_flank:
                    if is_positive:
                        chunk[0].start = max(
                            0, chunk[0].start - upstream_flank)
                    else:
                        chunk[-1].end = min(
                            lcontig,
                            chunk[-1].end + upstream_flank)
                if add_downstream_flank:
                    if is_positive:
                        chunk[-1].end = min(lcontig,
                                            chunk[-1].end + downstream_flank)
                    else:
                        chunk[0].start = max(
                            0, chunk[0].start - downstream_flank)
            else:
                if add_upstream_flank:
                    gff = GTF.Entry()
                    if is_positive:
                        gff.copy(chunk[0])
                        gff.end = gff.start
                        gff.start = max(0, gff.start - upstream_flank)
                        chunk.insert(0, gff)
                    else:
                        gff.copy(chunk[-1])
                        gff.start = gff.end
                        gff.end = min(lcontig, gff.end + upstream_flank)
                        chunk.append(gff)
                    gff.feature = "5-Flank"
                    gff.mMethod = "gff2gff"
                if add_downstream_flank:
                    gff = GTF.Entry()
                    if is_positive:
                        gff.copy(chunk[-1])
                        gff.start = gff.end
                        gff.end = min(lcontig, gff.end + downstream_flank)
                        chunk.append(gff)
                    else:
                        gff.copy(chunk[0])
                        gff.end = gff.start
                        gff.start = max(0, gff.start - downstream_flank)
                        chunk.insert(0, gff)
                    gff.feature = "3-Flank"
                    gff.mMethod = "gff2gff"

            if not is_positive:
                chunk.reverse()

            for gff in chunk:
                args.stdout.write(str(gff) + "\n")

    elif args.method == "complement-groups":

        iterator = GTF.joined_iterator(gffs,
                                       group_field=args.group_field)

        for chunk in iterator:
            if args.is_gtf:
                chunk = [x for x in chunk if x.feature == "exon"]
                if len(chunk) == 0:
                    continue
            chunk.sort(key=lambda x: (x.contig, x.start))
            x = GTF.Entry()
            x.copy(chunk[0])
            x.start = x.end
            x.feature = "intron"
            for c in chunk[1:]:
                x.end = c.start
                args.stdout.write(str(x) + "\n")
                x.start = c.end

    elif args.method == "combine-groups":

        iterator = GTF.joined_iterator(gffs,
                                       group_field=args.group_field)

        for chunk in iterator:
            chunk.sort(key=lambda x: (x.contig, x.start))
            x = GTF.Entry()
            x.copy(chunk[0])
            x.end = chunk[-1].end
            x.feature = "segment"
            args.stdout.write(str(x) + "\n")

    elif args.method == "join-features":
        for gff in combineGFF(gffs,
                              min_distance=args.min_distance,
                              max_distance=args.max_distance,
                              min_features=args.min_features,
                              max_features=args.max_features,
                              merge=False,
                              output_format=args.output_format):
            args.stdout.write(str(gff) + "\n")

    elif args.method == "merge-features":
        for gff in combineGFF(gffs,
                              min_distance=args.min_distance,
                              max_distance=args.max_distance,
                              min_features=args.min_features,
                              max_features=args.max_features,
                              merge=True,
                              output_format=args.output_format):
            args.stdout.write(str(gff) + "\n")

    elif args.method == "crop":
        for gff in cropGFF(gffs, args.filename_crop_gff):
            args.stdout.write(str(gff) + "\n")

    elif args.method == "crop-unique":
        for gff in cropGFFUnique(gffs):
            args.stdout.write(str(gff) + "\n")

    elif args.method == "filter-range":

        contig, strand, interval = None, None, None
        try:
            contig, strand, start, sep, end = re.match(
                "(\S+):(\S+):(\d+)(\.\.|-)(\d+)",
                args.filter_range).groups()
        except AttributeError:
            pass

        if not contig:
            try:
                contig, start, sep, end = re.match(
                    "(\S+):(\d+)(\.\.|-)(\d+)", args.filter_range).groups()
                strand = None
            except AttributeError:
                pass

        if not contig:
            try:
                start, end = re.match(
                    "(\d+)(\.\.|\,|\-)(\d+)", args.filter_range).groups()
            except AttributeError:
                raise "can not parse range %s" % args.filter_range
            contig = None
            strand = None

        if start:
            interval = (int(start), int(end))
        else:
            interval = None

        E.debug("filter: contig=%s, strand=%s, interval=%s" %
                (str(contig), str(strand), str(interval)))

        for gff in GTF.iterator_filtered(gffs, contig=contig,
                                         strand=strand,
                                         interval=interval):
            args.stdout.write(str(gff) + "\n")

    elif args.method == "sanitize":

        def assemblyReport(id):
            if id in assembly_dict.keys():
                id = assembly_dict[id]
            # if not in dict, the contig name is forced
            # into the desired convention, this is helpful user
            # modified gff files that contain additional contigs
            elif args.sanitize_method == "ucsc":
                if not id.startswith("contig") and not id.startswith("chr"):
                    id = "chr%s" % id
            elif args.sanitize_method == "ensembl":
                if id.startswith("contig"):
                    return id[len("contig"):]
                elif id.startswith("chr"):
                    return id[len("chr"):]
            return id

        if args.sanitize_method == "genome":
            if genome_fasta is None:
                raise ValueError(
                    "please specify --genome-file= when using "
                    "--sanitize-method=genome")
            f = genome_fasta.getToken
        else:
            if args.assembly_report is None:
                raise ValueError(
                    "please specify --assembly-report= when using "
                    "--sanitize-method=ucsc or ensembl")
            f = assemblyReport

        skipped_contigs = collections.defaultdict(int)
        outofrange_contigs = collections.defaultdict(int)
        filtered_contigs = collections.defaultdict(int)

        for gff in gffs:
            try:
                gff.contig = f(gff.contig)
            except KeyError:
                if args.skip_missing:
                    skipped_contigs[gff.contig] += 1
                    continue
                else:
                    raise

            if genome_fasta:
                lcontig = genome_fasta.getLength(gff.contig)
                if lcontig < gff.end:
                    outofrange_contigs[gff.contig] += 1
                    continue

            if args.contig_pattern:
                to_remove = [re.compile(x)
                             for x in args.contig_pattern.split(",")]
                if any([x.search(gff.contig) for x in to_remove]):
                    filtered_contigs[gff.contig] += 1
                    continue

            args.stdout.write(str(gff) + "\n")

        if skipped_contigs:
            E.info("skipped %i entries on %i contigs: %s" %
                   (sum(skipped_contigs.values()),
                    len(list(skipped_contigs.keys(
                    ))),
                    str(skipped_contigs)))

        if outofrange_contigs:
            E.warn("skipped %i entries on %i contigs because they are out of range: %s" %
                   (sum(outofrange_contigs.values()),
                    len(list(outofrange_contigs.keys())),
                    str(outofrange_contigs)))

        if filtered_contigs:
            E.info("filtered out %i entries on %i contigs: %s" %
                   (sum(filtered_contigs.values()),
                    len(list(filtered_contigs.keys())),
                    str(filtered_contigs)))

    elif args.method == "rename-chr":
        if not chr_map:
            raise ValueError("please supply mapping file")

        for gff in renameChromosomes(gffs, chr_map):
            args.stdout.write(str(gff) + "\n")

    else:

        for gff in gffs:

            if args.method == "forward_coordinates":
                gff.invert(contigs[gff.contig])

            if args.method == "forward_strand":
                gff.invert(contigs[gff.contig])
                gff.strand = "+"

            if agp:
                # note: this works only with forward coordinates
                gff.contig, gff.start, gff.end = agp.mapLocation(
                    gff.contig, gff.start, gff.end)

            args.stdout.write(str(gff) + "\n")

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
