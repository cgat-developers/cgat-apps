'''fasta2fasta.py - operate on sequences
=====================================

:Tags: Sequences

Purpose
-------

perform operations (masking, renaming) on a stream of fasta formatted sequences.

Available edit operations are:

translate
   translate sequences using the standard genetic code.

translate-to-stop
   translate until first stop codon

truncate-at-stop
   truncate sequence at first stop codon

back-translate
   convert nucleotide sequence to peptide sequence
   Requires parameter of second fasta file with peptide sequences.

mark-codons
   adds a space after each codon

apply-map
   rename sequence identifiers from a given map Requires parameter
   with filename of a map. The map is a tab-separated file mapping old
   to new names.

build-map
   rename sequence identifiers numerically and save output in a
   tab-separated file.  Requires parameter with filename of a map. The
   map is a tab-separated file mapping new to old names and will be
   newly created. Any exiting file of the same name will be
   overwritten.

pseudo-codons
   translate, but keep register with codons

interleaved-codons
   mix amino acids and codons

filter
   remove sequence according to certain criteria. For example,
   --method=filter --filter-method=min-length=5  --filter-method=max-length=10

map-codons:

remove-gaps
   remove all gaps in the sequence

mask-stops
   mask all stop codons

mask-seg
   mask sequence by running seg

mask-bias
   mask sequence by running bias

mask-codons
   mask codon sequence given a masked amino acid sequence.
   Requires parameter with masked amino acids in fasta format.

mask-incomplete-codons
   mask codons that are partially masked or gapped

mask-soft
   combine hard-masked (NNN) sequences with unmasked sequences to generate
   soft masked sequence (masked regions in lower case)

remove-stops
   remove stop codons

upper
   convert sequence to upper case

lower
   convert sequence to lower case

reverse-complement
   build the reverse complement

shuffle
   shuffle each sequence

sample
   select a certain proportion of sequences

Parameters are given to the option ``parameters`` in a comma-separated
list in the order that the edit operations are called upon.

Exclusion/inclusion is tested before applying any id mapping.

Usage
-----

Example::

   python fasta2fasta.py --method=translate < in.fasta > out.fasta

Type::

   python fasta2fasta.py --help

for command line help.

Command line options
---------------------

'''
import sys
import string
import re
import random
from itertools import zip_longest
import pysam

import cgatcore.experiment as E
import cgatcore.iotools as iotools
import cgat.Genomics as Genomics
import cgat.FastaIterator as FastaIterator
import cgat.Masker as Masker


def getCodons(sequence, gap_chars="-."):
    """get codons in sequence."""
    codons, codon = [], []
    full_codons, full_codon = [], []
    for x in range(len(sequence)):
        c = sequence[x]
        full_codon.append(c)
        if c not in gap_chars:
            codon.append(c)
        if len(codon) % 3 == 0:
            codons.append(codon)
            full_codons.append(full_codon)
            codon = []
            full_codon = []

    if full_codon:
        full_codons.append(full_codon)
        codons.append(codon)

    return full_codons, codons


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument(
        "-m", "--method", dest="methods", type=str, action="append",
        choices=("translate",
                 "translate-to-stop",
                 "truncate-at-stop",
                 "back-translate",
                 "mark-codons",
                 "apply-map",
                 "build-map",
                 "pseudo-codons",
                 "filter",
                 "interleaved-codons",
                 "map-codons",
                 "remove-gaps",
                 "mask-seg",
                 "mask-bias",
                 "mask-codons",
                 "mask-incomplete-codons",
                 "mask-stops",
                 "mask-soft",
                 "map-identifier",
                 "nop",
                 "remove-stops",
                 "upper",
                 "lower",
                 "reverse-complement",
                 "sample",
                 "shuffle"),
        help="method to apply to sequences.")

    parser.add_argument(
        "-p", "--parameters", dest="parameters", type=str,
        help="parameter stack for methods that require one ")

    parser.add_argument(
        "-x", "--ignore-errors", dest="ignore_errors", action="store_true",
        help="ignore errors.")

    parser.add_argument("--sample-proportion", dest="sample_proportion",
                        type=float,
                        help="sample proportion.")

    parser.add_argument(
        "--exclude-pattern", dest="exclude_pattern", type=str,
        help="exclude all sequences with ids matching pattern ")

    parser.add_argument(
        "--include-pattern", dest="include_pattern", type=str,
        help="include only sequences with ids matching pattern ")

    parser.add_argument(
        "--filter-method", dest="filter_methods", type=str,
        action="append",
        help="filtering methods to apply ")

    parser.add_argument(
        "-t", "--sequence-type", dest="type", type=str,
        choices=("aa", "na"),
        help="sequence type (aa or na) . This option determines "
        "which characters to use for masking.")

    parser.add_argument(
        "-l", "--template-identifier", dest="template_identifier",
        type=str,
        help="template for numerical identifier"
        "for the operation --build-map. A %%i is replaced by the position "
        "of the sequence in the file.")

    parser.add_argument(
        "--map-tsv-file", dest="map_tsv_file",
        type=str,
        help="input filename with map for identifiers. The first row is a header")

    parser.add_argument(
        "--fold-width", dest="fold_width", type=int,
        help="fold width for sequence output. 0 is unfolded ")

    parser.set_defaults(
        methods=[],
        parameters="",
        type="na",
        aa_mask_chars="xX",
        aa_mask_char="x",
        na_mask_chars="nN",
        na_mask_char="n",
        gap_chars="-.",
        gap_char="-",
        template_identifier="ID%06i",
        ignore_errors=False,
        exclude_pattern=None,
        include_pattern=None,
        sample_proportion=None,
        filter_methods=[],
        input_filename_fasta="-",
        input_filename_map=None,
        fold_width=80
    )

    (args, unknown) = E.start(parser,
                              unknowns=True)

    if len(unknown) > 0:
        args.input_filename_fasta = unknown[0]

    args.parameters = args.parameters.split(",")

    rx_include, rx_exclude = None, None
    if args.include_pattern:
        rx_include = re.compile(args.include_pattern)
    if args.exclude_pattern:
        rx_exclude = re.compile(args.exclude_pattern)

    iterator = FastaIterator.FastaIterator(args.stdin)

    nseq = 0

    map_seq2nid = {}

    map_identifier = ("apply-map" in args.methods or
                      "map-identifier" in args.methods)
    if map_identifier:
        if args.input_filename_map is None:
            raise ValueError("for method=map-identifier use --map-tsv-file")
        with iotools.open_file(args.input_filename_map) as infile:
            map_identifier = iotools.read_map(infile, has_header=True)

    if args.type == "na":
        mask_chars = args.na_mask_chars
        mask_char = args.na_mask_char
    else:
        mask_chars = args.aa_mask_chars
        mask_char = args.aa_mask_char

    if "map-codons" in args.methods:
        map_codon2code = iotools.ReadMap(open(args.parameters[0], "r"))
        del args.parameters[0]

    if "mask-soft" in args.methods:
        f = args.parameters[0]
        del args.parameters[0]
        hard_masked_iterator = FastaIterator.FastaIterator(open(f, "r"))

    if "mask-codons" in args.methods or "back-translate" in args.methods:

        # open a second stream to read sequences from
        f = args.parameters[0]
        del args.parameters[0]

        other_iterator = FastaIterator.FastaIterator(open(f, "r"))

    if "sample" in args.methods:
        if not args.sample_proportion:
            raise ValueError("specify a sample proportion")
        sample_proportion = args.sample_proportion
    else:
        sample_proportion = None

    filter_min_sequence_length = None
    filter_max_sequence_length = None
    filter_id_list = None
    for f in args.filter_methods:
        if f.startswith("min-length"):
            filter_min_sequence_length = int(f.split("=")[1])
        elif f.startswith("max-length"):
            filter_max_sequence_length = int(f.split("=")[1])
        elif f.startswith("id-file"):
            filter_id_list = [line[:-1] for line in iotools.open_file(f.split("=")[1])]

    def raiseIfNotCodon(l, title):
        '''raise ValueError if sequence length l is not divisible by
        3'''

        if l % 3 != 0:
            raise ValueError(
                "length of sequence %s not divisible by 3" % (title))

    iterator = pysam.FastxFile(args.input_filename_fasta)

    c = E.Counter()

    fold_width = args.fold_width

    def fold(s, w):
        return "\n".join([s[x:x+w] for x in range(0, len(s), w)])

    for record in iterator:
        c.nseq += 1
        c.input += 1

        sequence = re.sub(" ", "", record.sequence)
        l = len(sequence)

        if rx_include and not rx_include.search(record.name):
            c.skipped += 1
            continue

        if rx_exclude and rx_exclude.search(record.name):
            c.skipped += 1
            continue

        if sample_proportion:
            if random.random() > sample_proportion:
                continue

        if not (filter_id_list is None or record.name in filter_id_list):
            c.skipped += 1
            continue

        for method in args.methods:

            if method == "translate":
                # translate such that gaps are preserved
                seq = []

                ls = len(re.sub('[%s]' % args.gap_chars, sequence, ""))

                if ls % 3 != 0:
                    msg = "length of sequence %s (%i) not divisible by 3" % (
                        record.name, ls)
                    c.errors += 1
                    if args.ignore_errors:
                        E.warn(msg)
                        continue
                    else:
                        raise ValueError(msg)

                for codon in [sequence[x:x + 3] for x in range(0, l, 3)]:
                    aa = Genomics.MapCodon2AA(codon)
                    seq.append(aa)

                sequence = "".join(seq)

            elif method == "back-translate":
                # translate from an amino acid alignment to codon alignment
                seq = []

                try:
                    other_record = next(other_iterator)
                except StopIteration:
                    raise ValueError("run out of sequences")

                if record.name != other_record.title:
                    raise "sequence titles don't match: %s %s" % (
                        record.name, other_record.title)

                other_sequence = re.sub(
                    "[ %s]" % args.gap_chars, "", other_record.sequence)

                if len(other_sequence) % 3 != 0:
                    raise ValueError(
                        "length of sequence %s not divisible by 3" %
                        (other_record.title))

                r = re.sub("[%s]" % args.gap_chars, "", sequence)
                if len(other_sequence) != len(r) * 3:
                    raise ValueError(
                        "length of sequences do not match: %i vs %i" %
                        (len(other_sequence), len(r)))

                x = 0
                for aa in sequence:
                    if aa in args.gap_chars:
                        c = args.gap_char * 3
                    else:
                        c = other_sequence[x:x + 3]
                        x += 3
                    seq.append(c)

                sequence = "".join(seq)

            elif method == "pseudo-codons":
                raiseIfNotCodon(l, record.name)
                seq = []

                for codon in [sequence[x:x + 3] for x in range(0, l, 3)]:

                    aa = Genomics.MapCodon2AA(codon)
                    seq.append(aa)

                sequence = "   ".join(seq)

            elif method == "reverse-complement":
                sequence = sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

            elif method in ("mask-stops", "remove-stops"):
                c = []
                codon = []
                new_sequence = []

                if method == "mask-stops":
                    char = args.na_mask_char
                elif method == "remove-stops":
                    char = args.gap_char

                for x in sequence:

                    if x not in args.gap_chars:
                        codon.append(x.upper())

                    c.append(x)

                    if len(codon) == 3:
                        codon = "".join(codon).upper()
                        # mask all non-gaps
                        if Genomics.IsStopCodon(codon):

                            for x in c:
                                if x in args.gap_chars:
                                    new_sequence.append(x)
                                else:
                                    new_sequence.append(char)
                        else:
                            new_sequence += c

                        c = []
                        codon = []

                new_sequence += c

                sequence = "".join(new_sequence)

            elif method == "mask-soft":
                # Get next hard masked record and extract sequence and length
                try:
                    cur_hm_record = next(hard_masked_iterator)
                except StopIteration:
                    break
                hm_sequence = re.sub(" ", "", cur_hm_record.sequence)
                lhm = len(hm_sequence)
                new_sequence = []

                # Check lengths of unmasked and soft masked sequences the same
                if l != lhm:
                    raise ValueError(
                        "length of unmasked and hard masked sequences not "
                        "identical for record %s" %
                        (record.name))

                # Check if hard masked seq contains repeat (N), if so replace N
                # with lowercase sequence from unmasked version
                if sequence == hm_sequence:
                    pass
                else:
                    for x, y in zip_longest(sequence, hm_sequence):
                        if y == "N":
                            new_sequence += x.lower()
                        else:
                            new_sequence += x.upper()
                sequence = "".join(new_sequence)

            elif method == "map-codons":
                raiseIfNotCodon(l, record.name)
                seq = []

                for codon in (sequence[x:x + 3].upper()
                              for x in range(0, l, 3)):

                    if codon not in map_codon2code:
                        aa = "X"
                    else:
                        aa = map_codon2code[codon]
                    seq.append(aa)

                sequence = "".join(seq)

            elif method == "interleaved-codons":
                raiseIfNotCodon(l, record.name)
                seq = []

                for codon in [sequence[x:x + 3] for x in range(0, l, 3)]:

                    aa = Genomics.MapCodon2AA(codon)
                    seq.append("%s:%s" % (aa, codon))

                sequence = " ".join(seq)

            elif method == "translate-to-stop":
                seq = []

                for codon in [sequence[x:x + 3] for x in range(0, l, 3)]:

                    if Genomics.IsStopCodon(codon):
                        break

                    aa = Genomics.MapCodon2AA(codon)
                    seq.append(aa)

                sequence = "".join(seq)

            elif method == "truncate-at-stop":
                seq = []

                for codon in [sequence[x:x + 3] for x in range(0, l, 3)]:

                    if Genomics.IsStopCodon(codon):
                        break
                    seq.append(codon)

                sequence = "".join(seq)

            elif method == "remove-gaps":

                seq = []
                for s in sequence:
                    if s in args.gap_chars:
                        continue
                    seq.append(s)

                sequence = "".join(seq)

            elif method == "upper":
                sequence = sequence.upper()

            elif method == "lower":
                sequence = sequence.lower()

            elif method == "mark-codons":
                raiseIfNotCodon(l, record.name)
                seq = []

                sequence = " ".join([sequence[x:x + 3]
                                     for x in range(0, l, 3)])

            elif method == "apply-map":
                id = re.match("^(\S+)", record.name).groups()[0]
                if id in map_seq2nid:
                    rest = record.name[len(id):]
                    record.name = map_seq2nid[id] + rest

            elif method == "build-map":
                # build a map of identifiers
                id = re.match("^(\S+)", record.name).groups()[0]
                new_id = args.template_identifier % nseq
                if id in map_seq2nid:
                    raise "duplicate fasta entries - can't map those: %s" % id
                map_seq2nid[id] = new_id
                record.name = new_id

            elif method == "mask-bias":
                masker = Masker.MaskerBias()
                sequence = masker(sequence)

            elif method == "mask-seg":
                masker = Masker.MaskerSeg()
                sequence = masker(sequence)

            elif method == "shuffle":
                s = list(sequence)
                random.shuffle(s)
                sequence = "".join(s)

            elif method == "mask-incomplete-codons":
                seq = list(sequence)
                for x in range(0, l, 3):
                    nm = len([x for x in seq[x:x + 3] if x in mask_chars])
                    if 0 < nm < 3:
                        seq[x:x + 3] = [mask_char] * 3
                sequence = "".join(seq)

            elif method == "mask-codons":
                # mask codons based on amino acids given as reference
                # sequences.
                other_record = next(other_iterator)

                if other_record is None:
                    raise ValueError("run out of sequences.")

                if record.name != other_record.title:
                    raise ValueError(
                        "sequence titles don't match: %s %s" %
                        (record.name, other_record.title))

                other_sequence = re.sub(" ", "", other_record.sequence)

                if len(other_sequence) * 3 != len(sequence):
                    raise ValueError(
                        "sequences for %s don't have matching lengths %i - %i" %
                        (record.name, len(other_sequence) * 3,
                         len(sequence)))

                seq = list(sequence)
                c = 0
                for x in other_sequence:
                    if x in args.aa_mask_chars:
                        if x.isupper():
                            seq[c:c + 3] = [args.na_mask_char.upper()] * 3
                        else:
                            seq[c:c + 3] = [args.na_mask_char.lower()] * 3
                    c += 3

                sequence = "".join(seq)

        l = len(sequence)
        if filter_min_sequence_length is not None and \
           l < filter_min_sequence_length:
            c.skipped += 1

        if filter_max_sequence_length is not None and \
           l > filter_max_sequence_length:
            c.skipped += 1
            continue

        record.sequence = sequence
        if fold_width >= 0:
            if record.comment:
                args.stdout.write(">{} {}\n{}\n".format(
                    record.name,
                    record.comment,
                    fold(record.sequence, fold_width)))
            else:
                args.stdout.write(">{}\n{}\n".format(
                    record.name,
                    fold(record.sequence, fold_width)))
        else:
            args.stdout.write(str(record) + "\n")

        c.output += 1

    if "build-map" in args.methods:
        p = args.parameters[0]
        if p:
            outfile = iotools.open_file(p, "w")
        else:
            outfile = args.stdout

        outfile.write("old\tnew\n")
        for old_id, new_id in list(map_seq2nid.items()):
            outfile.write("%s\t%s\n" % (old_id, new_id))
        if p:
            outfile.close()

    E.info(c)
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
