'''
split_fasta.py
======================================================

:Tags: Python

Purpose
-------

.. todo::

   describe purpose of the script.

Usage
-----

Example::

   python split_fasta.py --help

Type::

   python split_fasta.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import os
import cgat.FastaIterator as FastaIterator
import cgatcore.iotools as iotools
import cgatcore.experiment as E


class Files:

    mFiles = {}

    def __init__(self,
                 output_pattern=None,
                 skip_identifiers=False):

        self.mOutputPattern = output_pattern
        self.mSkipIdentifiers = skip_identifiers
        self.mCounts = {}

    def __del__(self):
        """close all open files."""
        for file in list(self.mFiles.values()):
            file.close()

    def GetFile(self, identifier):
        return identifier

    def GetFilename(self, identifier):
        """get filename for an identifier."""

        if self.mOutputPattern:
            return re.sub("%s", str(identifier), self.mOutputPattern)
        else:
            return identifier

    def OpenFile(self, filename, mode="w"):
        """open file.

        If file is in a new directory, create directories.
        """
        if mode in ("w", "a"):
            dirname = os.path.dirname(filename)
            if dirname and not os.path.exists(dirname):
                os.makedirs(dirname)

        returniotools.open_file(filename, mode)

    def Write(self, identifier, sequence):

        filename = self.GetFilename(identifier)

        if filename not in self.mFiles:

            if len(self.mFiles) > 1000:
                for f in list(self.mFiles.values()):
                    f.close()
                self.mFiles = {}

            self.mFiles[filename] = self.OpenFile(filename, "a")

        if self.mSkipIdentifiers:
            self.mFiles[filename].write("%s\n" % (sequence.sequence))
        else:
            self.mFiles[filename].write(
                ">%s\n%s\n" % (sequence.title, sequence.sequence))

        if filename not in self.mCounts:
            self.mCounts[filename] = 0
        self.mCounts[filename] += 1

    def DeleteFiles(self, min_size=0):
        """delete all files below a minimum size."""

        ndeleted = 0
        for filename, counts in list(self.mCounts.items()):
            if counts < min_size:
                os.remove(filename)
                ndeleted += 1

        return ndeleted


class FilesChunks(Files):

    def __init__(self,
                 chunk_size, **kwargs):

        Files.__init__(self, **kwargs)
        self.mChunkSize = chunk_size
        self.mFilename = 0

    def GetFilename(self, identifier):

        if not self.mFilename or self.mCounts[self.mFilename] % self.mChunkSize == 0:
            self.mFilename = re.sub(
                "%s", str(len(self.mCounts) + 1), self.mOutputPattern)

        return self.mFilename


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument("-f", "--file", dest="input_filename", type=str,
                        help="input filename. If not given, stdin is used.")

    parser.add_argument("-i", "--input-pattern", dest="input_pattern", type=str,
                        help="input pattern. Parses description line in order to extract id.")

    parser.add_argument("-o", "--output-filename-pattern", dest="output_pattern", type=str,
                        help="output pattern. Gives filename for a given sequence.")

    parser.add_argument("-n", "--num-sequences", dest="num_sequences", type=int,
                        help="split by number of sequences (not implemented yet).")

    parser.add_argument("-m", "--map", dest="map_filename", type=str,
                        help="map filename. Map identifiers to filenames")

    parser.add_argument("-s", "--skip-identifiers", dest="skip_identifiers", action="store_true",
                        help="do not write identifiers.")

    parser.add_argument("--min-size", dest="min_size", type=int,
                        help="minimum cluster size.")

    parser.set_defaults(
        input_filename=None,
        map_filename=None,
        skip_identifiers=False,
        input_pattern="^(\S+)",
        min_size=0,
        num_sequences=None,
        output_pattern="%s")

    (args) = E.start(parser)

    if args.input_filename:
        infile = iotools.open_file(args.input_filename, "r")
    else:
        infile = sys.stdin

    if args.map_filename:
        map_id2filename = iotools.ReadMap(open(args.map_filename, "r"))
    else:
        map_id2filename = {}

    if args.num_sequences:
        files = FilesChunks(chunk_size=args.num_sequences,
                            output_pattern=args.output_pattern,
                            skip_identifiers=args.skip_identifiers)

    else:
        files = Files(output_pattern=args.output_pattern,
                      skip_identifiers=args.skip_identifiers)

    if args.input_pattern:
        rx = re.compile(args.input_pattern)
    else:
        rx = None

    ninput = 0
    noutput = 0
    identifier = None
    chunk = 0

    for seq in FastaIterator.iterate(infile):

        ninput += 1

        if rx:
            try:
                identifier = rx.search(seq.title).groups()[0]
            except AttributeError:
                print("# parsing error in description line %s" % (seq.title))
        else:
            identifier = seq.title

        if map_id2filename:
            if identifier in map_id2filename:
                identifier = map_id2filename[identifier]
            else:
                continue

        files.Write(identifier, seq)
        noutput += 1

    if args.input_filename:
        infile.close()

    # delete all clusters below a minimum size
    # Note: this has to be done at the end, because
    # clusters sizes are only available once both the fasta
    # file and the map has been parsed.
    if args.min_size:
        ndeleted = files.DeleteFiles(min_size=args.min_size)
    else:
        ndeleted = 0

    if args.loglevel >= 1:
        print("# input=%i, output=%i, ndeleted=%i" % (ninput, noutput, ndeleted))

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
