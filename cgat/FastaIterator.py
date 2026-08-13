"""FastaIterator.py - Iteration over fasta files
================================================

This module provides a simple iterator of Fasta formatted files.  The
difference to the biopython iterator is that the iterators in this
module skip over comment lines starting with "#".

.. note::
   Another way to access the information in :term:`fasta` formatted
   files is through pysam_.

Reference
---------

"""
import os
import subprocess


class FastaRecord:
    """a :term:`fasta` record.

    Attributes
    ----------
    title: string
       the title of the sequence

    sequence: string
       the sequence

    fold : int
       the number of bases per line when writing out
    """

    def __init__(self, title, sequence, fold=False):

        self.title = title
        self.sequence = sequence
        self.fold = fold

    def __str__(self):
        ''' str method for writing out'''

        if self.fold:
            seq = [self.sequence[i:i + self.fold]
                   for i in range(0, len(self.sequence), self.fold)]
        else:
            seq = (self.sequence,)

        return ">%s\n%s" % (self.title, "\n".join(seq))


class FastaIterator:
    '''a iterator of :term:`fasta` formatted files.

    Yields
    ------
    FastaRecord

    '''

    def __init__(self, f, *args, **kwargs):
        self.iterator = iterate(f)

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.iterator)

    def next(self):
        return next(self.iterator)


def iterate(infile, comment="#", fold=False):
    '''iterate over fasta data in infile

    Lines before the first fasta record are
    ignored (starting with ``>``) as well as
    lines starting with the comment character.

    Parameters
    ----------
    infile : File
        the input file
    comment : char
        comment character
    fold : int
        the number of bases before line split when writing out

    Yields
    ------
    FastaRecord
    '''
    h = infile.readline()[:-1]

    if not h:
        return

    # skip everything until first fasta entry starts
    while h[0] != ">":
        h = infile.readline()[:-1]
        if not h:
            return
        continue

    h = h[1:]
    seq = []

    for line in infile:

        if line.startswith(comment):
            continue

        if line.startswith('>'):
            yield FastaRecord(h, ''.join(seq), fold)

            h = line[1:-1]
            seq = []
            continue

        seq.append(line[:-1])

    yield FastaRecord(h, ''.join(seq), fold)


def iterate_together(*args):
    """iterate synchronously over one or more fasta files.

    The iteration finishes once any of the files is exhausted.

    Arguments
    ---------

    :term:`fasta`-formatted files to be iterated upon

    Yields
    ------
    tuple
       a tuple of :class:`FastaRecord` corresponding to
       the current record in each file.
    """

    iterators = [FastaIterator(x) for x in args]
    for records in zip(*iterators):
        yield list(records)


def count(filename):
    '''count number of sequences in fasta file.

    Uses ``grep -c`` (and ``gzip -dc`` for ``.gz`` files) so that
    counting stays fast on large FASTA files.

    Arguments
    ---------
    filename : string
        The filename

    Raises
    ------
    OSError
        If the file does not exist

    Returns
    -------
    int
        The number of sequences in the file.
    '''

    if not os.path.exists(filename):
        raise OSError("file '%s' does not exist" % filename)

    grep = ["grep", "-c", "^>"]
    try:
        if filename.endswith(".gz"):
            gzip_proc = subprocess.Popen(
                ["gzip", "-dc", filename],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL)
            result = subprocess.run(
                grep,
                stdin=gzip_proc.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False)
            if gzip_proc.stdout is not None:
                gzip_proc.stdout.close()
            gzip_proc.wait()
        else:
            result = subprocess.run(
                grep + [filename],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False)
    except FileNotFoundError:
        raise OSError(
            "grep/gzip is required to count sequences in '%s'" % filename)

    # grep -c exits 1 when there are no matches
    if result.returncode == 1:
        return 0
    if result.returncode != 0:
        stderr = result.stderr.decode() if result.stderr else ""
        raise OSError(
            "failed to count sequences in '%s': %s" % (filename, stderr))
    return int(result.stdout.strip() or 0)
