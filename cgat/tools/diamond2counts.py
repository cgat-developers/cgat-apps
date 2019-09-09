'''
diamond2counts.py - count alignments to reference
===================================================

:Tags: Python

Purpose
-------

Count the number of alignments to each reference in outfmt6.

Counts are based on various options specified by --method.

best       This will take the best alignment as judged by the highest
           bitscore.




TODO::
Add additional options

Usage
-----

Example::

   python diamond2counts.py

Type::

   python diamond2counts.py --help

for command line help.

Command line options
--------------------

'''

import sys

import cgatcore.experiment as E
from cgat.Diamond import *
import collections
import cgatcore.iotools as iotools


def readCogMap(cog_map):
    '''
    return a dictionary mapping gene to cog
    '''
    gene2cog = {}
    for line in iotools.open_file(cog_map):
        data = line[:-1].split("\t")
        gene2cog[data[0]] = data[1]
    return gene2cog


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--version", action='version', version="1.0")

    parser.add_argument("-m", "--method", dest="method", type=str,
                        choices=("best", None),
                        help="method for determing what to count")

    parser.add_argument("--sum-cog", dest="sum_cog", action="store_true",
                        help="sum counts over functions (COGs) in --cog-map")

    parser.add_argument("--evaluate-cog",
                        dest="evaluate_cog",
                        action="store_true",
                        help="""output the percent of
                              alignments for each read = best hit""")

    parser.add_argument("--cog-map", dest="cog_map", type=str,
                        help="file with gene to cog map")

    parser.add_argument("-n", "--nsamples", dest="nsamples", type=int,
                        help="""number of queries to evaluate-
                              will take the first n in the file""")

    parser.set_defaults(method=None,
                        sum_cog=False,
                        evaluate_cog=False,
                        cog_map=None,
                        nsamples=10000)

    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    if args.evaluate_cog:
        assert args.cog_map, """must specify an annotation
                                   mapping gene to function (COG)"""
        assert not args.method, """evaluation performed
                                      in the absence of counting"""

        E.info("reading gene to function (COG) map %s" % args.cog_map)
        gene2cog = readCogMap(args.cog_map)
        E.info("loaded gene to function (COG) map")

        E.info("retrieving alignment data")
        args.stdout.write("query\tpbest\tnalignments\n")

        c = 0
        for alignments in query_iterator(alignment_iterator(args.stdin)):
            c += 1
            scores = []
            if c <= args.nsamples:
                for alignment in alignments:
                    scores.append(alignment.score)
                    best = max(scores)
                    best_alignments = [
                        x for x in alignments if x.score == best]
                    if len(best_alignments) > 1:
                        best_alignments = random.sample(best_alignments, 1)
                    best_alignment = best_alignments[0]
                    best_cog = gene2cog[best_alignment.ref]
                pbest = float(len(
                    [gene2cog[x.ref]
                     for x in alignments
                     if gene2cog[x.ref] == best_cog])) / len(alignments) * 100
                nalignments = len(alignments)
                args.stdout.write(
                    "\t".join(map(
                        str, [alignments[0].qid,
                              pbest,
                              nalignments])) + "\n"
                )
            else:
                break
        return

    # container for counts
    counts = collections.defaultdict(int)
    E.info("counting alignments")
    assert args.method, "required option --method"
    if args.method == "best":
        if args.sum_cog:
            E.warn("""summing over functions (COGS)
                      will remove genes with no annotations
                      and those with multiple COG assignments""")
            assert args.cog_map, """a mapping between gene and
                                       function (COG) is required"""

            E.info("""reading gene to function (COG) mapping from %s"""
                   % args.cog_map)
            gene2cog = readCogMap(args.cog_map)
            E.info("loaded gene to function (COG) mapping")

            E.info("summing functional assignments")
            query_it = query_iterator(alignment_iterator(args.stdin))
            for best in best_alignment_iterator(query_it):
                cog = gene2cog[best.ref]
                # removing uassigned or multiple assignments
                if cog == "unknown" or cog.find(";") != -1:
                    continue
                counts[cog] += 1
        else:
            E.info("counting best alignments")
            query_it = query_iterator(alignment_iterator(args.stdin))
            for best in best_alignment_iterator(query_it):
                counts[best.ref] += 1
        E.info("finished counting")

        E.info("writing results")
        args.stdout.write("ref\tcount\n")
        for ref, count in sorted(counts.items()):
            args.stdout.write("\t".join([ref, str(count)]) + "\n")

    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
