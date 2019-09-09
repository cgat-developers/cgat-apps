"""plot-variant-stats
=====================

Implemeted methods
------------------

plot-mutational-profile
+++++++++++++++++++++++

Plot mutational profile from output of vcf-stats.

"""

import sys
import pandas
import numpy
import re
import matplotlib.pyplot as plt
import cgatcore.experiment as E
from cgat.Plots.VariantPlots import MutationProfileBarPlot, DepthProfilePlot, \
    ManhattanPlot


def plot_mutation_profile_bar_plot(dataframe, section, map_key2label={}, **kwargs):

    for key, dataframe in dataframe.groupby(by="sample"):
        if key == "unique":
            continue

        if dataframe.empty:
            E.warn("no data for {}".format(key))
            continue

        ax = MutationProfileBarPlot()(dataframe)

        label = map_key2label.get(key, key)
        plt.savefig(E.get_output_file("-".join((section, label))))
        plt.close()


def plot_depth_profile_plot(dataframe, section, map_key2label={}, **kwargs):

    ax = DepthProfilePlot()(dataframe, map_sample2label={})
    plt.savefig(E.get_output_file(section))
    plt.close()


def plot_manhattan_plot(dataframe,
                        section,
                        filename_fasta,
                        map_key2label={},
                        **kwargs):

    plotter = ManhattanPlot(genome_size_file=filename_fasta)
    ax = plotter(dataframe, **kwargs)
    plt.savefig(E.get_output_file(section))
    plt.close()


def compute_log_depth_ratio(dataframe, min_depth=10):

    table = dataframe.set_index(["CHROM", "POS"])

    if len(table.columns) != 2:
        raise NotImplementedError("expected 2 columns in table {}".format(fn))

    normed_table = table / table.sum()
    # assumption: column order is normal, tumour
    normal, tumour = normed_table.columns
    logratios = numpy.log2(normed_table[tumour] / normed_table[normal])

    normed_table["l2fold_DP"] = logratios
    result = normed_table.drop([normal, tumour], axis=1)

    result = result[table[normal].gt(min_depth) |
                    table[tumour].gt(min_depth)]
    return result.reset_index()


def main(argv=None):

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-m", "--method", dest="method", type=str,
        choices=["mutation-profile-bar-plot",
                 "depth-profile-line-plot",
                 "manhattan-plot"],
        help="methods to apply ")

    parser.add_argument(
        "-t", "--transformation", dest="transformations", type=str,
        action="append",
        choices=["log-depth-ratio"],
        help="dataframe transformation options ")

    parser.add_argument(
        "-r", "--regex-filename", dest="regex_filename", type=str,
        help="")

    parser.add_argument(
        "-f", "--reference-fasta-file", dest="reference_fasta_file",
        help="reference genomic sequence in fasta format. "
        )

    parser.add_argument(
        "--input-file-format", dest="input_file_format", type=str,
        choices=("tsv", "bcftools-query"),
        help="input file format "
        )

    parser.add_argument(
        "--plot-options", dest="plot_options", type=str,
        help="plot options to pass through to the plotter. The string is "
        "eval'ed, for example: --plot-options='window_size=3, ylabel=\"12\"' "
        )

    parser.set_defaults(
        method=None,
        reference_fasta=None,
        input_file_format="tsv",
        plot_options=None,
        transformations=[],
    )

    (args, unknown) = E.start(parser,
                              argv=argv,
                              add_output_options=True,
                              unknowns=True)

    filenames = unknown

    if len(filenames) == 0:
        E.info("reading from stdin")
        filenames = [args.stdin]

    if args.plot_options is not None:
        plot_options = eval("dict({})".format(args.plot_options))
    else:
        plot_options = {}

    for index, filename in enumerate(filenames):

        E.info("working on {}".format(filename))

        try:
            if args.input_file_format == "bcftools-query":
                # for bctools query, header starts with "#".

                dataframe = pandas.read_csv(filename,
                                            sep="\t",
                                            skip_blank_lines=False,
                                            header=0,
                                            dtype={"CHROM": str})
                # names are of format [1]sample1:DP, extract sample1
                dataframe.columns = (
                    [re.search("\[\d+\]([^:]+)", x).groups()[0]
                     for x in dataframe.columns])
            else:
                dataframe = pandas.read_csv(filename, sep="\t",
                                            dtype={"CHROM": str})
        except pandas.io.common.EmptyDataError:
            E.warn("no data in {}, skipped".format(filename))
            continue

        E.info("read data from {}".format(filename))

        if args.regex_filename:
            section = re.search(args.regex_filename, filename).groups()[0]
        else:
            section = "{}".format(index + 1)

        for method in args.transformations:
            if method == "log-depth-ratio":
                dataframe = compute_log_depth_ratio(dataframe)

        if dataframe.empty:
            E.warn("dataframe from {} is empty - skipped".format(filename))
            continue

        if args.method == "mutation-profile-bar-plot":
            plot_mutation_profile_bar_plot(dataframe, section, **plot_options)

        elif args.method == "depth-profile-line-plot":
            plot_depth_profile_plot(dataframe, section, **plot_options)

        elif args.method == "manhattan-plot":
            plot_manhattan_plot(dataframe,
                                section,
                                filename_fasta=args.reference_fasta_file,
                                **plot_options)

    E.stop()


if __name__ == "__main__":
    sys.exit(main())
