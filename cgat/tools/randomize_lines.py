'''
randomize_lines.py - randomize lines from stdin
===============================================

:Tags: Python

Purpose
-------

This script reads lines from stdin and outputs them
in randomized order.

Usage
-----

Example::

   cgat randomize-lines < in.lines > out.lines

Command line options
--------------------

'''

import sys
import random
import cgatcore.experiment as E


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("-k", "--keep-header", dest="keep_header", type=int,
                        help="randomize, but keep header in place ")

    parser.set_defaults(keep_header=0)

    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    inf = args.stdin
    outf = args.stdout
    c = E.Counter()
    for x in range(args.keep_header):
        c.header += 1
        outf.write(inf.readline())

    lines = inf.readlines()
    c.lines_input = len(lines)
    random.shuffle(lines)
    for line in lines:
        outf.write(line)
    c.lines_output = len(lines)

    E.info(c)

    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
