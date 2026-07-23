"""
bcl2fastq.py - Wrapper that outputs fastq files from bcl files
==============================================================

:Tags: BCL FASTQ Conversion

Purpose
-------

Convert the raw data from an Illumina Sequencing Run to :term:`fastq` formatted files.

Usage
-----

Example::

   python bcl2fastq.py -p "--runfolder-dir <RunFolder>"

This command demultiplexes and converts BCL files in the given run folder directory.
All arguments for Illumina's bcl2fastq software must be given with the -p argument. 
Type::

   python bcl2fastq.py --help

for command line help.

Documentation
-------------

Converts BCL files in a given run folder to fastq.

Command line options
--------------------

``--arguments``
    Supply arguments to be passed to Illumia's bcl2fastq software.

``--output-dir``
    Required if using --fastqc.

``--fastqc``
    After converting BCL files, run all fastq files in FastQC.

``--fastqc-options``
    Supply arguments to be passed to FastQC.

``--bcl2fastq-help``
    Prints help for Illumina's bcl2fastq software.

"""

import cgatcore.experiment as E
import subprocess
import sys
import glob
import gzip
import shlex
import shutil


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    parser = E.ArgumentParser()

    parser.add_argument("-p", "--arguments", type=str, dest="arguments",   
                        default="",
                        help="Pass options and arguments to the executable. Please surround options in \"\"")

    parser.add_argument("-o", "--output-dir", type=str, dest="output",
                        default=".",
                        help="Output for the fastq files.")

    parser.add_argument("-f", "--fastqc", dest="fastqc",
                        action="store_true",
                        help="After demultiplexing open the fastq files in FastQC.")

    parser.add_argument("-F", "--fastqc-options", type=str, dest="fastqc_options",
                        default="",
                        help="Options for FastQC. Please surround options in \"\"")

    parser.add_argument("-H", "--bcl2fastq-help", dest="bcl2fastq_help", 
                        action="store_true",
                        help="Print help for Illumina's bcl2fastq conversion software")

    (args) = E.start(parser)

    if shutil.which("bcl2fastq") is None:
        raise ValueError("bcl2fastq cannot be found")

    if args.bcl2fastq_help:
        subprocess.run(["bcl2fastq", "--help"], check=False)
        return

    bcl2fastq_cmd = ["bcl2fastq", "-o", args.output]
    if args.arguments:
        bcl2fastq_cmd[1:1] = shlex.split(args.arguments)
    result = subprocess.run(bcl2fastq_cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"bcl2fastq failed with exit code {result.returncode}")

    for infile in glob.glob(f"{args.output}/**/*.fastq.gz", recursive=True):
        with gzip.GzipFile(f"{infile}", "r") as f:
            if sum(1 for char in f.read().decode('utf-8') if char == "\n") % 4 != 0:
                raise ValueError(f"{infile} is either corrupt or incomplete.")

    if args.fastqc:
        fastqc_options = shlex.split(args.fastqc_options) if args.fastqc_options else []
        for infile in glob.glob(f"{args.output}/**/*.fastq.gz", recursive=True):
            subprocess.run(["fastqc", infile] + fastqc_options, check=False)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
