'''
fastas2fasta.py - concatenate sequences from multiple fasta files
==========================================================================

:Tags: Genomics Sequences MultipleAlignments FASTA Manipulation

Purpose
-------

This script reads sequences from two or more :term:`fasta` formatted
files and outputs a new file with the sequences concatenated per
entry.

All files must have the same number of sequences and the id of
the first file is output.

Usage
-----

Example::

   python fastas2fasta.py a.fasta b.fasta > c.fasta

If a.fasta is::

  >1
  AAACC
  >2
  CCCAA

and b.fasta is::

  >a
  GGGGTTT
  >b
  TTTTGGG

then the output will be::

  >1
  AAACCGGGGTTT
  >2
  CCCAATTTTGGG


Type::

   python fastas2fasta.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re

import cgatcore.experiment as E
import cgatcore.iotools as iotools
import cgat.FastaIterator as FastaIterator


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    (args, unknown) = E.start(parser, unknowns=True)

    if len(unknown) < 2:
        raise ValueError(
            "please supply at least two filenames to concatenate.")

    iterators = []
    for a in unknown:
        iterators.append(FastaIterator.FastaIterator(iotools.open_file(a, "r")))

    ninput, noutput, nerrors = 0, 0, 0

    while 1:

        sequences = []
        ids = []

        for iterator in iterators:
            try:
                cur_record = next(iterator)
            except StopIteration:
                break

            sequences.append(re.sub(" ", "", cur_record.sequence))
            ids.append(cur_record.title)

        if not sequences:
            break
        ninput += 1

        if len(sequences) != len(iterators):
            raise ValueError("unequal number of sequences in files")

        noutput += 1

        args.stdout.write(">%s\n%s\n" % (ids[0],
                                         "".join(sequences)))

    E.info("ninput=%i, noutput=%i, nerrors=%i" % (ninput, noutput, nerrors))

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
