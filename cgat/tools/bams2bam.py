'''bams2bam.py - merge genomic and transcriptome mapped bamfiles
====================================================================

:Tags: Genomics NGS Geneset BAM Manipulation

Purpose
-------

This script takes as input two BAM files from an RNASeq experiment.
The first bam file (:file:`bamG`) should contain reads mapped against
the genome using a mapper permitting splicing (e.g. tophat). The
second bam file (:file:`bamT`) should contain reads mapped against
known transcripts. This script will write a new bam file that removes
reads from :term:`bamG` that map to regions that are conflicting with
those in :term:`bamT`.

.. note::
   Note that if junctions are supplied, the resultant bam files will not
   be sorted by position.

.. glossary::

   bamG
      :term:`bam` formatted file with reads mapped against the genome

   bamT
      :term:`bam` formatted file with reads mapped against transcripts

Usage
-----

Example::

   python bams2bam.py bamT.bam bamG.bam

Type::

   python bams2bam.py --help

for command line help.

Documentation
-------------

The script needs to look-up reads via their names. It thus builds an
index of reads mapping

This script requires the NM attributes to be set. If it is not set,
you will need to set a policy.

Command line options
--------------------

'''

import os
import sys
import pysam

import cgatcore.experiment as E
import cgat.GTF as GTF
import cgatcore.iotools as iotools
import cgat.Bed as Bed
import cgat.IndexedGenome as IndexedGenome
from cgat.BamTools.bamtools import bams2bam_filter


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(version="%prog version: $Id$",
                              usage=globals()["__doc__"])

    parser.add_argument(
        "-g", "--gtf-file", dest="filename_gtf", type=str,
        help="filename with gene models in gtf format ")

    parser.add_argument(
        "-m", "--filename-mismapped", dest="filename_mismapped", type=str,
        help="output bam file for mismapped reads ")

    parser.add_argument(
        "-j", "--junctions-bed-file", dest="filename_junctions", type=str,
        help="bam file with reads mapped across junctions ")

    parser.add_argument(
        "-r", "--filename-regions", dest="filename_regions", type=str,
        help="filename with regions to remove in bed format ")

    parser.add_argument(
        "-t", "--transcripts-gtf-file", dest="filename_transcriptome",
        type=str,
        help="bam file with reads mapped against transcripts ")

    parser.add_argument(
        "-p", "--map-tsv-file", dest="filename_map", type=str,
        help="filename mapping transcript numbers (used by "
        "--filename-transciptome) to transcript names "
        "(used by --filename-gtf) ")

    parser.add_argument(
        "-s", "--filename-stats", dest="filename_stats", type=str,
        help="filename to output stats to ")

    parser.add_argument(
        "-o", "--colour",
        dest="colour_mismatches", action="store_true",
        help="mismatches will use colour differences (CM tag) ")

    parser.add_argument(
        "-i", "--ignore-mismatches",
        dest="ignore_mismatches", action="store_true",
        help="ignore mismatches ")

    parser.add_argument(
        "-c", "--remove-contigs", dest="remove_contigs", type=str,
        help="','-separated list of contigs to remove ")

    parser.add_argument(
        "-f", "--force-output", dest="force", action="store_true",
        help="force overwriting of existing files ")

    parser.add_argument("-u", "--unique", dest="unique", action="store_true",
                        help="remove reads not matching uniquely ")

    parser.add_argument("--output-sam", dest="output_sam", action="store_true",
                        help="output in sam format ")

    parser.set_defaults(
        filename_gtf=None,
        filename_mismapped=None,
        filename_junctions=None,
        filename_transcriptome=None,
        filename_map=None,
        remove_contigs=None,
        force=False,
        unique=False,
        colour_mismatches=False,
        ignore_mismatches=False,
        output_sam=False,
        filename_table=None,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)

    if len(args) != 1:
        raise ValueError("please supply one bam file")

    bamfile_genome = args[0]
    genome_samfile = pysam.AlignmentFile(bamfile_genome, "rb")

    if options.remove_contigs:
        options.remove_contigs = options.remove_contigs.split(",")

    if options.filename_map:
        E.info("reading map")
        id_map = iotools.read_map(
            iotools.open_file(options.filename_map), has_header=True)
        id_map = dict([(y, x) for x, y in id_map.items()])
    else:
        id_map = None

    transcripts = {}
    if options.filename_gtf:
        E.info("indexing geneset")
        mapped, missed = 0, 0
        for gtf in GTF.transcript_iterator(
                GTF.iterator(iotools.open_file(options.filename_gtf))):
            gtf.sort(key=lambda x: x.start)
            transcript_id = gtf[0].transcript_id
            if id_map:
                try:
                    transcript_id = id_map[transcript_id]
                    mapped += 1
                except KeyError:
                    missed += 1
                    continue
            transcripts[transcript_id] = gtf

        E.info("read %i transcripts from geneset (%i mapped, %i missed)" %
               (len(transcripts), mapped, missed))

    regions_to_remove = None
    if options.filename_regions:
        E.info("indexing regions")
        regions_to_remove = IndexedGenome.Simple()
        for bed in Bed.iterator(iotools.open_file(options.filename_regions)):
            regions_to_remove.add(bed.contig, bed.start, bed.end)
        E.info("read %i regions" % len(regions_to_remove))

    if options.filename_transcriptome:
        transcripts_samfile = pysam.AlignmentFile(options.filename_transcriptome,
                                                  "rb")
    else:
        transcripts_samfile = None

    if options.output_sam:
        output_samfile = pysam.AlignmentFile("-", "wh", template=genome_samfile)
    else:
        output_samfile = pysam.AlignmentFile("-", "wb", template=genome_samfile)

    if options.filename_mismapped:
        if not options.force and os.path.exists(options.filename_mismapped):
            raise IOError("output file %s already exists" %
                          options.filename_mismapped)
        output_mismapped = pysam.AlignmentFile(options.filename_mismapped,
                                               "wb",
                                               template=genome_samfile)
    else:
        output_mismapped = None

    if options.filename_junctions:
        junctions_samfile = pysam.AlignmentFile(options.filename_junctions,
                                                "rb")
    else:
        junctions_samfile = None

    c = bams2bam_filter(genome_samfile,
                        output_samfile,
                        output_mismapped,
                        transcripts_samfile,
                        junctions_samfile,
                        transcripts,
                        regions=regions_to_remove,
                        unique=options.unique,
                        remove_contigs=options.remove_contigs,
                        colour_mismatches=options.colour_mismatches,
                        ignore_mismatches=options.ignore_mismatches,
                        ignore_transcripts=transcripts_samfile is None,
                        ignore_junctions=junctions_samfile is None)

    if options.filename_stats:
        outf = iotools.open_file(options.filename_stats, "w")
        outf.write("category\tcounts\n%s\n" % c.asTable())
        outf.close()

    if options.filename_transcriptome:
        transcripts_samfile.close()

    genome_samfile.close()
    output_samfile.close()
    if output_mismapped:
        output_mismapped.close()

    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
