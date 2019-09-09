'''
transfac2transfac.py - filter transfac motif files
====================================================

:Tags: Python

Purpose
-------

Filter a transfac motif file.

Usage
-----

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''


import sys
import re
import cgatcore.experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-f", "--filter-prefix", dest="filter_prefix", default=None,
        help="ID prefix to filter on, eg. V for vertebrates")

    parser.add_argument(
        "-p", "--pattern-identifier", dest="filter_pattern", default=None,
        help="ID pattern to filter (filter is case insensitive) eg. pax6. "
        "Multiple patterns should be specified as a comma separated list")

    (args) = E.start(parser)

    if args.filter_pattern:
        patterns = [x.strip() for x in args.filter_pattern.split(",")]
        E.info("Supplied patterns %s" % ", ".join(patterns))
    else:
        patterns = False

    filtered_motifs = []
    n = 0

    inmotif, tid, filter_emit, pattern_emit = False, False, False, False

    for line in args.stdin:

        # pick up motif start and ends.
        if line.startswith("AC") and inmotif is False:
            # print "in align"
            inmotif = True
            motif = line
            continue
        elif line.startswith("ID") and inmotif is True:
            # print line
            tid = line.split("  ")[1]
            motif += line
            continue

        elif line.startswith("//") and inmotif is True:

            motif += line

            if tid is False:
                raise ValueError("matrix ID not determined")

            if args.filter_prefix:
                if tid.startswith(args.filter_prefix):
                    filter_emit = True
            else:
                filter_emit = True

            if patterns is not False:
                for pat in patterns:
                    match = re.search(pat, tid, re.IGNORECASE)
                    if match is not None:
                        pattern_emit = True
                        break
            else:
                pattern_emit = True

            if filter_emit is True and pattern_emit is True:
                filtered_motifs.append(motif)
                n += 1

            inmotif, tid, filter_emit, pattern_emit = (
                False, False, False, False)
            continue

        elif inmotif is True:
            motif += line

        elif inmotif is False:
            continue

        else:
            raise ValueError("unknown parsing state")

    args.stdout.write("".join(filtered_motifs))

    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
