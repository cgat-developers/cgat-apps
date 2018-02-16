'''fasta2fasta - convert fasta files
===================================================

This tool manipulates :term:`fasta` formatted input
files.

Command line options
--------------------

'''
import sys
import pysam

import CGATCore.Experiment as E
import CGATCore.IOTools as IOTools


def main(argv=None):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-f", "--fasta", dest="input_filename_fasta",
        type="string",
        help="filename with fasta sequences. ")

    parser.add_option(
        "-o", "--output-filename-sequences", dest="output_filename_sequences",
        type="string",
        help="output per sequence information to filename")

    parser.add_option(
        "-m", "--method", dest="methods", action="append",
        type="choice",
        choices=["map-identifier", "nop"],
        help="method to apply")

    parser.add_option(
        "--input-filename-map", dest="input_filename_map",
        type="string",
        help="input filename with map for identifiers. The first row is a header")

    parser.add_option(
        "--fold-width", dest="fold_width", type="int",
        help="fold width for sequence output. 0 is unfolded [%default]")

    parser.set_defaults(
        input_filename_fasta="-",
        methods=[],
        input_filename_map=None,
        fold_width=80
    )

    (options, args) = E.start(parser, argv=argv)

    if len(args) > 0:
        options.input_filename_fasta = args[0]

    if options.methods is None:
        raise ValueError("please specify at least one method")

    map_identifier = "map-identifier" in options.methods
    if map_identifier:
        if options.input_filename_map is None:
            raise ValueError("for method=map-identifier use --input-filename-map")
        with IOTools.open_file(options.input_filename_map) as infile:
            map_identifier = IOTools.read_map(infile, has_header=True)

    iterator = pysam.FastxFile(options.input_filename_fasta)

    c = E.Counter()

    fold_width = options.fold_width

    def fold(s, w):
        return "\n".join([s[x:x+w] for x in range(0, len(s), w)])

    method = options.method

    for record in iterator:
        c.input += 1

        if map_identifier:
            record.name = map_identifier.get(record.name, record.name)

        if fold_width >= 0:
            options.stdout.write(">{} {}\n{}\n".format(
                record.name,
                record.comment,
                fold(record.sequence, fold_width)))
        else:
            options.stdout.write(str(record) + "\n")

        c.output += 1

    E.info(c)
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
