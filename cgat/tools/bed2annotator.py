'''
bed2annotator.py - convert bed to annotator format
==================================================

:Tags: Python

Purpose
-------

This script converts a bed file into annotator compatible regions. Depending on the option --section
this script will create:

   segments
      a segments file

   annotations
      a file with annotations. Each bed track is a separate annotation.

   workspace
      a file with a workspace

Usage
-----

Example::

   python bed2annotator2tsv.py --help

Type::

   python bed2annotator2tsv.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import itertools
import collections

import cgatcore.experiment as E
import cgat.Bed as Bed
import cgat.IndexedFasta as IndexedFasta


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("-g", "--genome-file", dest="genome_file", type=str,
                        help="filename with genome.")

    parser.add_argument("-f", "--features", dest="features", type=str,
                        help="feature to collect.")

    parser.add_argument("-i", "--files", dest="files", action="append",
                        help="use multiple annotations.")

    parser.add_argument("-a", "--annotations", dest="annotations", type=str,
                        help="aggregate name for annotations if only single file is provided from STDIN.")

    parser.add_argument("--map-tsv-file", dest="input_filename_map", type=str,
                        help="filename with a map of gene_ids to categories.")

    parser.add_argument("-l", "--max-length", dest="max_length", type=str,
                        help="maximum segment length.")

    parser.add_argument("-m", "--merge-overlapping", dest="merge", action="store_true",
                        help="merge overlapping bed segments.")

    parser.add_argument("-s", "--section", dest="section", type=str,
                        choices=("segments", "annotations", "workspace"),
                        help="annotator section.")

    parser.add_argument("--subset", dest="subsets", type=str, action="append",
                        help="add filenames to delimit subsets within the gff files. The syntax is filename.gff,label,filename.ids.")

    parser.set_defaults(
        genome_file=None,
        feature=None,
        remove_random=True,
        section="segments",
        annotations="annotations",
        max_length=100000,
        files=[],
        subsets=[],
        input_filename_map=None,
        merge=False,
    )

    (args, unknown) = E.start(parser,
                              unknowns=True)

    args.files += unknown
    if len(args.files) == 0:
        args.files.append("-")
    args.files = list(
        itertools.chain(*[re.split("[,; ]+", x) for x in args.files]))

    if args.subsets:
        subsets = collections.defaultdict(list)
        for s in args.subsets:
            filename_gff, label, filename_ids = s.split(",")
            subsets[filename_gff].append((label, filename_ids))
        args.subsets = subsets

    if args.genome_file:
        fasta = IndexedFasta.IndexedFasta(args.genome_file)
    else:
        fasta = None

    if args.section == "segments":
        prefix = "##Segs"
    elif args.section == "annotations":
        prefix = "##Id"
    elif args.section == "workspace":
        prefix = "##Work"
    else:
        raise ValueError("unknown section %s" % args.section)

    if args.max_length:
        max_length = args.max_length
    else:
        max_length = 0

    ninput, ntracks, ncontigs, nsegments, ndiscarded = 0, 0, 0, 0, 0

    if args.section in ("annotations"):
        contigs = set()
        it = itertools.groupby(
            Bed.iterator(args.stdin), key=lambda x: x.track["name"])

        map_track2segments = {}
        for track, beds in it:
            ntracks += 1
            map_track2segments[track] = []
            first_segment = nsegments

            beds = list(beds)

            if args.merge:
                beds = Bed.merge(beds)

            for bed in beds:
                contig, start, end = bed.contig, bed.start, bed.end

                if args.remove_random and "random" in contig:
                    continue

                if max_length > 0 and end - start > max_length:
                    ndiscarded += 1
                    continue

                contigs.add(contig)
                map_track2segments[track].append(nsegments)
                args.stdout.write(
                    "%s\t%i\t%s\t(%i,%i)\n" % (prefix, nsegments, contig, start, end))
                nsegments += 1

            args.stdout.write("##Ann\t%s\t%s\n" % (
                track, "\t".join(["%i" % x for x in range(first_segment, nsegments)])))
            E.info("track %s: annotated with %i segments" %
                   (track, nsegments - first_segment))

        ncontigs = len(contigs)
        E.info("ninput=%i, ntracks=%i, ncontigs=%i, nsegments=%i, ndiscarded=%i" %
               (ninput, ntracks, ncontigs, nsegments, ndiscarded))

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
