'''
clean.py - clean up output files from aborted runs
==================================================

:Tags: Python

Purpose
-------

This script checks one or more output files have they
have completed successfully. It will remove output files
for those jobs that are incomplete.

The script checks for the "job finished" tag at the
end of the file.

Usage
-----

Example::

   python clean.py --help

Type::

   python clean.py --help

for command line help.

Command line options
--------------------
'''

import os
import sys
import re
import glob
import os.path
import cgatcore.experiment as E


def getLastLine(filename, read_size=1024):
    """return last line of a file.
    """
    f = iotools.open_file(
        filename, 'rU')    # U is to open it with Universal newline support
    offset = read_size
    f.seek(0, 2)
    file_size = f.tell()
    if file_size == 0:
        return ""
    while 1:
        if file_size < offset:
            offset = file_size
        f.seek(-1 * offset, 2)
        read_str = f.read(offset)
        # Remove newline at the end
        if read_str[offset - 1] == '\n':
            read_str = read_str[:-1]
        lines = read_str.split('\n')
        if len(lines) >= 2:
            return lines[-1]
        if offset == file_size:   # reached the beginning
            return read_str
        offset += read_size
    f.close()


def checkPythonRuns(filename):
    """returns true if a python run is complete."""
    last_line = getLastLine(filename)
    return re.match("# job finished", last_line)


def isNewer(a, b):
    """return true if file a is newer than file b."""

    # get times of most recent access
    at = os.stat(a)[7]
    bt = os.stat(b)[7]

    return at > bt


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("-g", "--glob", dest="glob_pattern", type=str,
                        help="glob pattern to use for collecting files .")

    parser.add_argument("-n", "--dry-run", dest="dry_run", action="store_true",
                        help="only print out actions, do not execute them .")

    parser.add_argument("-f", "--file-pattern", dest="file_pattern", type=str,
                        help="only check files matching this pattern .")

    parser.set_defaults(glob_pattern="data.dir",
                        file_pattern=".out",
                        check_completeness="python",
                        skip_dirs=[],
                        dry_run=False,
                        )

    args, unknowns = E.start(parser,
                             add_pipe_options=True, unknowns=True)

    if unknowns:
        starts = unknowns
    elif args.glob_pattern:
        starts = glob.glob(args.glob_pattern)
    else:
        starts = "."

    ndirs, nfiles, ndeleted = 0, 0, 0

    if args.check_completeness == "python":
        isComplete = checkPythonRuns

    rx = re.compile(args.file_pattern)

    for start in starts:
        for root, dirs, files in os.walk(start):

            ndirs += 1
            # exclude directories
            for dir in args.skip_dirs:
                if dir in dirs:
                    dirs.remove(dir)

            for filename in files:
                p = os.path.join(root, filename)
                if rx.search(filename) and not isComplete(p):
                    if args.dry_run:
                        args.stdlog.write("# removing file %s\n" % p)
                    else:
                        os.remove(p)
                    ndeleted += 1

    if args.loglevel >= 1:
        args.stdlog.write("# ndirs=%i, nfiles=%i, ndeleted=%i\n" %
                          (ndirs, nfiles, ndeleted))

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
