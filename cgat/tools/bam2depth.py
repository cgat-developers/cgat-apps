'''output depth statistics for a BAM file.
'''

import collections
import subprocess
import re
import os
import shlex

import cgatcore.experiment as E
import cgatcore.iotools as iotools


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--input-filename-fasta", dest="input_filename_fasta", type=str,
        help="filename with reference sequence in fasta format ")

    parser.add_argument(
        "--counting-mode", dest="counting_mode", type=str,
        choices=("all", "pileup_defaults"),
        help="counting mode. all=all reads/bases. pileup-defaults= "
        "use default pileup thresholds. Options will be added to "
        "--mpileup-options. .")

    parser.add_argument(
        "--mpileup-options", dest="mpileup_options", type=str,
        help="pileup options to use ")

    parser.set_defaults(
        mpileup_options="",
        counting_mode="all",
        input_filename_fasta=None,
        report_step=1000000,
    )

    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv, add_output_options=True)

    bamfile = args[0]

    mpileup_options = args.mpileup_options

    if args.counting_mode == "all":
        mpileup_options += " -Q 0 -B -A"

    read_depth_histogram = collections.defaultdict(int)
    base_depth_histogram = collections.defaultdict(int)

    # deletions are marked by something like -2AA at the first
    # position and a '*' for subsequent positions
    rx_deletions = re.compile("([-][0-9]+|[*])")
    report_step = args.report_step
    npositions = 0

    samtools = iotools.which("samtools")

    statement = (
        "{samtools} mpileup "
        "-f {reference_fasta} "
        "{mpileup_options} "
        "{bamfile} ".format(
            samtools=samtools,
            reference_fasta=args.input_filename_fasta,
            mpileup_options=mpileup_options,
            bamfile=os.path.abspath(bamfile)))

    E.info("running the following statement: {}".format(statement))

    cmd_args = shlex.split(statement)
    proc = subprocess.Popen(
        cmd_args,
        shell=False,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        cwd=os.path.abspath(os.curdir))

    for line in proc.stdout:
        line = line.decode("utf-8")
        contig, pos, base, read_depth, info, qualities = line[:-1].split("\t")
        read_depth = int(read_depth)
        pos = int(pos)

        if pos % report_step == 0:
            E.info("working on {}: {}".format(contig, pos))

        ndeletions = len(rx_deletions.findall(info))
        base_depth = read_depth - ndeletions

        read_depth_histogram[read_depth] += 1
        base_depth_histogram[base_depth] += 1

    for line in proc.stderr:
        E.warn(line)

    keys = sorted(set(read_depth_histogram.keys()).union(
        base_depth_histogram.keys()))

    args.stdout.write("depth\tread_depth_positions\tbase_depth_positions\n")
    for key in keys:
        args.stdout.write("{}\t{}\t{}\n".format(
                key,
                read_depth_histogram[key],
                base_depth_histogram[key]))

    E.info("positions tested: {}".format(sum(read_depth_histogram.values())))
    E.stop()
