'''
bed2fasta.py - get sequences from bed file
==========================================

:Tags: Genomics Intervals Sequences Conversion BED FASTA



Purpose
-------

This script outputs nucleotide sequences for intervals within
a :term:`bed` formatted file using a corresponding genome file.

Usage
-----

A required input to bed2fasta.py is a cgat indexed genome. To obtain an
idexed human reference genome we would type

Example::
   cat hg19.fasta | index_fasta.py hg19 > hg19.log

This file would then serve as the --genome-file when we wish to extract
sequences from a :term:`bed` formatted file.


For example we could now type::

   cat in.bed | python bed2fasta.py --genome-file hg19 > out.fasta

Where we take a set of genomic intervals (e.g. from a human ChIP-seq experiment)
and output their respective nucleotide sequences.


Type::

   python bed2fasta.py --help

for command line help.

Command line options
--------------------

'''
import sys
import cgatcore.experiment as E
import cgat.Bed as Bed
import cgat.IndexedFasta as IndexedFasta
import cgat.Masker as Masker


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("-g", "--genome-file", dest="genome_file", type=str,
                        help="filename with genomic sequence to retrieve "
                        "sequences from.")

    parser.add_argument("-m", "--masker", dest="masker", type=str,
                        choices=("dust", "dustmasker", "softmask", "none"),
                        help="apply masker to mask output sequences ")

    parser.add_argument("--output-mode", dest="output_mode", type=str,
                        choices=("intervals", "leftright", "segments"),
                        help="what to output. "
                        "'intervals' generates a single sequence for "
                        "each bed interval. 'leftright' generates two "
                        "sequences, one in each direction, for each bed "
                        "interval. 'segments' can be used to output "
                        "sequence from bed12 files so that sequence only covers "
                        "the segments ")

    parser.add_argument("--min-sequence-length", dest="min_length", type=int,
                        help="require a minimum sequence length ")

    parser.add_argument("--max-sequence-length", dest="max_length", type=int,
                        help="require a maximum sequence length ")

    parser.add_argument(
        "--extend-at", dest="extend_at", type=str,
        choices=("none", "3", "5", "both", "3only", "5only"),
        help="extend at 3', 5' or both or no ends. If 3only or 5only "
        "are set, only the added sequence is returned ")

    parser.add_argument(
        "--extend-by", dest="extend_by", type=int,
        help="extend by # bases")

    parser.add_argument(
        "--use-strand", dest="ignore_strand",
        action="store_false",
        help="use strand information and return reverse complement "
        "on intervals located on the negative strand. ")

    parser.set_defaults(
        genome_file=None,
        masker=None,
        output_mode="intervals",
        min_length=0,
        max_length=0,
        extend_at=None,
        extend_by=100,
        ignore_strand=True,
    )

    (args) = E.start(parser)

    if args.genome_file:
        fasta = IndexedFasta.IndexedFasta(args.genome_file)
        contigs = fasta.getContigSizes()
        fasta.setConverter(IndexedFasta.getConverter("zero-both-open"))

    counter = E.Counter()
    ids, seqs = [], []

    E.info("collecting sequences")
    for bed in Bed.setName(Bed.iterator(args.stdin)):
        counter.input += 1

        lcontig = fasta.getLength(bed.contig)

        if args.ignore_strand:
            strand = "+"
        else:
            strand = bed.strand

        # Apply the --extend-by logic
        if args.extend_by:
            if args.extend_at in ("both", "3", "3only"):
                bed.end = min(bed.end + args.extend_by, lcontig)
            if args.extend_at in ("both", "5", "5only"):
                bed.start = max(0, bed.start - args.extend_by)

        if args.output_mode == "segments" and bed.columns == 12:
            ids.append("%s %s:%i..%i (%s) %s %s" %
                       (bed.name, bed.contig, bed.start, bed.end, strand,
                        bed["blockSizes"], bed["blockStarts"]))
            seg_seqs = [fasta.getSequence(bed.contig, strand, start, end)
                        for start, end in bed.toIntervals()]
            seqs.append("".join(seg_seqs))

        elif (args.output_mode == "intervals" or
              args.output_mode == "segments"):
            ids.append("%s %s:%i..%i (%s)" %
                       (bed.name, bed.contig, bed.start, bed.end, strand))
            seqs.append(
                fasta.getSequence(bed.contig, strand, bed.start, bed.end))

        elif args.output_mode == "leftright":
            l = bed.end - bed.start

            start, end = max(0, bed.start - l), bed.end - l
            ids.append("%s_l %s:%i..%i (%s)" %
                       (bed.name, bed.contig, start, end, strand))
            seqs.append(fasta.getSequence(bed.contig, strand, start, end))

            start, end = bed.start + l, min(lcontig, bed.end + l)
            ids.append("%s_r %s:%i..%i (%s)" %
                       (bed.name, bed.contig, start, end, strand))
            seqs.append(fasta.getSequence(bed.contig, strand, start, end))

    E.info("collected %i sequences" % len(seqs))

    masked = Masker.maskSequences(seqs, args.masker)
    args.stdout.write(
        "\n".join([">%s\n%s" % (x, y) for x, y in zip(ids, masked)]) + "\n")

    E.info("masked %i sequences" % len(seqs))

    counter.output = len(seqs)

    E.info("%s" % counter)

    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
