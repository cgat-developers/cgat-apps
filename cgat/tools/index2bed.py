"""
index2bed.py - convert indexed fasta file to bed file
=====================================================

:Tags: Python

Purpose
-------

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Command line options
--------------------

"""

import sys
import re
import cgat.IndexedFasta as IndexedFasta
import cgatcore.experiment as E


def getFixedWidthWindows(map_contig2size, window_size):
    """return a list of fixed contig sizes."""

    for contig, size in list(map_contig2size.items()):
        E.info("processing %s" % contig)
        for x in range(0, size, window_increment):
            if x + window_size > size:
                continue
            gff = GTF.Entry()
            gff.feature = "window"
            gff.source = "window"
            gff.contig = contig
            gff.start = x
            gff.end = min(size, x + window_size)
            yield gff


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument(
        "-g", "--genome-file", dest="genome_file", type=str,
        help="filename with genome.")

    parser.add_argument(
        "--remove-regex", dest="remove_regex",
        type=str,
        help="regular expression of contigs to remove.")

    parser.add_argument(
        "-e", "--gff-file", dest="gff_file", type=str,
        help="gff file to use for getting contig sizes.")

    parser.add_argument(
        "-f", "--fixed-width-windows",
        dest="fixed_width_windows", type=str,
        help="fixed width windows. Supply the window size as a "
        "parameter. Optionally supply an offset.")

    parser.set_defaults(
        genome_file=None,
        remove_regex=None,
        fixed_windows=None,
    )

    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    if args.remove_regex:
        remove_regex = re.compile(args.remove_regex)
    else:
        remove_regex = None

    if args.fixed_width_windows:
        v = list(map(int, args.fixed_width_windows.split(",")))
        if len(v) == 2:
            window_size, window_increment = v
        elif len(v) == 1:
            window_size, window_increment = v[0], v[0]
        else:
            raise ValueError(
                "could not parse window size '%s': should be size[,increment]" % args.fixed_width_windows)

    if args.gff_file:
        infile = iotools.open_file(args.gff_file, "r")
        gff = GTF.readFromFile(infile)
        infile.close()
        for g in gff:
            try:
                map_contig2size[g.mName] = max(map_contig2size[g.mName], g.end)
            except ValueError:
                map_contig2size[g.mName] = g.end

    else:
        gff = None

    if args.genome_file:
        fasta = IndexedFasta.IndexedFasta(args.genome_file)
        map_contig2size = fasta.getContigSizes(with_synonyms=False)
    else:
        fasta = None

    if map_contig2size is None:
        raise ValueError("no source of contig sizes supplied")

    # do sth
    counter = E.Counter()

    for contig, size in list(map_contig2size.items()):
        size = int(size)
        counter.input += 1

        if remove_regex and remove_regex.search(contig):
            counter.skipped += 1
            continue

        if args.fixed_width_windows:
            for x in range(0, size, window_increment):
                if x + window_size > size:
                    continue
                args.stdout.write(
                    "%s\t%i\t%i\n" % (contig, x, min(size, x + window_size)))
                counter.windows += 1
        else:
            args.stdout.write("%s\t%i\t%i\n" % (contig, 0, size))
            counter.windows += 1

        counter.output += 1

    E.info(str(counter))

    # write footer and output benchmark information.
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
