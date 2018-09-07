import matplotlib
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy
import os
import pandas
import seaborn

import cgatcore.experiment as E


class MutationProfileBarPlot(object):
    """plot a mutation profile bar-plot.

    dataframe requires the following columns:

    signature:
    context:
    percent_sample
    """

    def __call__(self, dataframe):

        assert "signature" in dataframe.columns
        assert "context" in dataframe.columns
        assert "percent_sample" in dataframe.columns

        current_palette = seaborn.color_palette()
        factors = dataframe.signature.unique()
        map_signature2factor = dict([(y, x) for x, y in
                                     enumerate(factors)])

        colors = [current_palette[map_signature2factor[x]]
                  for x in dataframe.signature]

        handles = []
        for x, factor in enumerate(factors):
            handles.append(
                matplotlib.patches.Patch(color=current_palette[x],
                                         label=factor))

        ax = dataframe.plot(kind="bar",
                            x="context",
                            y="percent_sample",
                            color=colors,
                            legend=False)
        ax.legend(handles=handles)

        ax.figure.set_size_inches(30, 10)
        ax.set_xticklabels([x.get_text() for x in ax.get_xticklabels()],
                           family="monospace")

        return ax


class DepthProfilePlot(object):
    """
    """

    def __call__(self, dataframe, map_sample2label={}):

        df = dataframe.pivot(
            index="gc_bin",
            columns="sample",
            values="mean").reset_index()

        # remove duplicate sample names
        # (in Cancer analysis: two blood samples)
        to_drop = [x for x in df.columns if x.startswith("2:")]
        df.drop(to_drop, axis=1, inplace=True)

        df.colums = [map_sample2label.get(x, x) for x in df.columns]

        if df.empty:
            E.warn("no data, no plot will be output")
            return

        ax = df.plot(kind="line", x="gc_bin")

        return ax


class ChromosomeDataFromFile(object):
    """create melted data frame

    Table format is:

    CHROM   POS   sample1   sample2    ...
    """

    def __call__(self, track, value_label):
        fn = track

        table = pandas.read_csv(fn,
                                comment="#",
                                sep="\t",
                                dtype={"CHROM": object},
                                )

        table = pandas.melt(table, id_vars=["CHROM", "POS"])
        table.colums = ["CHROM", "POS", "sample", value_label]

        try:
            skip = table["depth"].str.startswith(".").fillna(False)
            table = table[~skip]
        except AttributeError:
            pass
        return table.set_index(["CHROM", "POS"])


class CoverageRatioFromFile(object):

    min_depth = 10

    def __call__(self, track, column_is_norm=None, *args, **kwargs):

        fn = track

        table = pandas.read_csv(fn,
                                comment="#",
                                sep="\t",
                                dtype={"CHROM": object},
                                ).set_index(["CHROM", "POS"])

        columns = table.columns

        if len(columns) > 2:
            E.warn("too many columns {}".format(columns))

        if column_is_norm:
            if len(columns) == 2:
                num, den = columns
                if num in column_is_norm:
                    num, den = den, num
            else:
                den = [x for x in columns if x in column_is_norm][0]
                num = [x for x in columns if x not in column_is_norm][0]

            if den not in column_is_norm:
                raise ValueError("denominator is {}, but not norm, cols = {}".format(
                    den, columns))

        table = table[(table[num] >= self.min_depth) &
                      (table[den] >= self.min_depth)]
        for column in columns:
            table[column] = table[column] / table[column].median()

        table["sample"] = num
        table["copy number"] = 2.0 * table[num] / table[den]
        table.drop(columns, axis=1, inplace=True)
        return table


class ChromosomePlot(object):

    position = None
    column = None
    plot_type = None
    window_size = 100
    rolling = None
    subplots = False
    plot_gaps = True

    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, dataframe, chromosome):

        if self.position is None or self.column is None:
            raise NotImplementedError(
                "missing .position or .column attributes")

        plt.figure()
        ax = plt.gca()
        colors = seaborn.color_palette()
        data = dataframe.reset_index()
        labels, plots = [], []
        x = 0
        grouper = data.groupby("track")

        if self.plot_gaps:
            gaps_fn = PARAMS.get("gaps_file", None)
            if gaps_fn is None or gaps_fn == "None":
                gaps = None
            else:
                gaps = pandas.read_csv(gaps_fn,
                                       comment="#",
                                       sep="\t",
                                       header=None)
                gaps = gaps[gaps[0] == chromosome][[0, 1, 2]]

        for name, group in grouper:
            if self.subplots:
                plt.subplot(len(list(grouper.groups.keys())), 1, x + 1)

            if self.plot_type == "scatter":
                plots.append(plt.scatter(group[self.position],
                                         group[self.column],
                                         s=20,
                                         c=colors[x % len(colors)]))
            elif self.plot_type == "line":
                plots.append(plt.plot(group[self.position],
                                      group[self.column],
                                      c=colors[x % len(colors)])[0])
            labels.append(name)
            if self.rolling == "average":
                mean = pandas.rolling_mean(group[self.column],
                                           self.window_size)
                plt.plot(group[self.position],
                         mean,
                         c=colors[x % len(colors)])

            x += 1
            plt.xlim(0, max(data[self.position]))

            plt.tick_params(axis='both', which='major', labelsize=20)
            plt.tick_params(axis='both', which='minor', labelsize=20)

            if gaps is not None:
                ax = plt.gca()
                ymin, ymax = ax.get_ylim()
                for i, row in gaps.iterrows():
                    ax.add_patch(
                        patches.Rectangle(
                            (row[1], ymin),
                            row[2] - row[1],
                            ymax - ymin,
                            color="black",
                            alpha=0.2))

        if self.subplots:
            plt.subplot(len(list(grouper.groups.keys())), 1, 1)
            plt.legend(
                plots, labels,
                bbox_to_anchor=(0., 1.0, 1., .102),
                loc=3,
                ncol=3,
                mode="expand",
                borderaxespad=0.,
                prop={'size': 24}
            )
        else:
            plt.legend(plots, labels)

        ax.figure.set_size_inches(40, 10 * len(list(grouper.groups.keys())))
        return ax


@E.cached_function
def get_chromosome_sizes(filename, vcf_chromosomes=None):

    if filename.endswith(".fa") or filename.endswith(".fasta"):
        fn = filename + ".fai"
        if not os.path.exists(fn):
            raise OSError("file {} with chromosome sizes does not exist"
                          .format(fn))

        with open(fn) as inf:
            chromosome_size = pandas.read_csv(
                inf,
                sep="\t",
                header=None,
                names=("chromosome", "size", "offset",
                       "lw", "lw2"),
                dtype={"chromosome": str},
                index_col="chromosome").drop(["offset", "lw", "lw2"], 1)

    elif filename.endswith(".bed") or filename.endswith(".bed.gz"):

        if not os.path.exists(filename):
            raise OSError("file {} with exons does not exist"
                          .format(filename))
        chromosome_size = pandas.read_csv(
            filename,
            sep="\t",
            header=None,
            names=("chromosome", "start", "end"),
            dtype={"chromosome": str},
        )
        chromosome_size["size"] = (chromosome_size["end"] -
                                   chromosome_size["start"] + 1)

        chromosome_size = chromosome_size.groupby(
            "chromosome").sum().drop(["start", "end"], 1)

    # intersect with VCF
    if vcf_chromosomes:
        chromosome_size = chromosome_size.loc[vcf_chromosomes]

    return chromosome_size.sort_values("size", ascending=False)


class ManhattanPlot(object):

    def __init__(self, genome_size_file, *args, **kwargs):
        if genome_size_file is None:
            raise ValueError(
                "ManhattanPlot requires a file with contig sizes.")

        if not os.path.exists(genome_size_file):
            raise OSError(
                "filename to retrieve contig sizes from does not exist: {}".format(genome_size_file))

        self.genome_size_file = genome_size_file

    def __call__(self,
                 dataframe,
                 chromosome_column="CHROM",
                 position_column="POS",
                 value_column=None,
                 plot_type="scatter",
                 plot_gaps=None,
                 subplots=True,
                 ylabel="read depth",
                 xlabel="chromosome position",
                 rolling=None,
                 window_size=100,
                 yrange=None):

        if chromosome_column not in dataframe.columns:
            raise ValueError(
                "missing chromosome column {}".format(chromosome_column))

        if position_column not in dataframe.columns:
            raise ValueError(
                "missing position column {}".format(position_column))

        if "sample" not in dataframe.columns:
            data = pandas.melt(
                dataframe,
                id_vars=[chromosome_column, position_column]).reset_index().drop(
                    ["index"], axis=1)

        else:
            if value_column is None:
                raise ValueError("a melted table requires value_column to be provided")

            data = dataframe[[chromosome_column,
                              position_column,
                              "sample",
                              value_column]]

        data.columns = [chromosome_column, position_column,
                        "sample", "value"]

        try:
            skip = data["value"].str.startswith(".").fillna(False)
            data = data[~skip]
        except AttributeError:
            pass

        plt.figure()
        ax = plt.gca()
        colors = seaborn.color_palette()
        labels, plots = [], []

        # build increment for chromosomal positions
        chromosome_sizes = get_chromosome_sizes(self.genome_size_file)
        chromosome_offsets = chromosome_sizes.cumsum().shift(1).fillna(0)
        chromosome_offsets_dict = dict(
            [(str(x), int(y)) for x, y in chromosome_offsets.iterrows()])

        # remove variants on unknown chromosomes
        known_chromosomes = set(chromosome_sizes.index)
        data = data[data[chromosome_column].isin(known_chromosomes)]

        if len(data) == 0:
            raise ValueError("after chromosome filtering, no data remained")

        # sort by order in chromonome sizes
        data[chromosome_column] = pandas.Categorical(
            data[chromosome_column].astype(str),
            chromosome_sizes.index)

        data.sort_values(by=[chromosome_column, position_column],
                         inplace=True)
        data[position_column] = [y + chromosome_offsets_dict[str(x)]
                                 for x, y in zip(data[chromosome_column],
                                                 data[position_column])]

        if plot_gaps:
            if not os.path.exits(plot_gaps):
                raise ValueError("filename {} with gap data not found".format(
                    plot_gaps))
            gaps = pandas.read_csv(plot_gaps,
                                   comment="#",
                                   sep="\t",
                                   header=None)
            gaps.columns = ["contig", "start", "end"]
            gaps = gaps[gaps.contig.isin(chromosome_offsets_dict)]
            gaps.start = [y + chromosome_offsets_dict[str(x)]
                          for x, y in zip(gaps.contig,
                                          gaps.start)]
            gaps.end = [y + chromosome_offsets_dict[str(x)]
                        for x, y in zip(gaps.contig,
                                        gaps.end)]
        else:
            gaps = None

        chromosome_colors = seaborn.color_palette(
            "Paired",
            len(chromosome_offsets))

        column = "value"
        grouper = data.groupby(by="sample")
        subp = 0
        for name, group in grouper:
            if subplots:
                plt.subplot(len(list(grouper.groups.keys())), 1, subp + 1)

            if yrange is not None:
                if yrange.startswith("percentile-"):
                    percentile = float(yrange.split("-")[1])
                    ymin = 0
                    ymax = numpy.percentile(group[column],
                                            percentile)
                elif yrange.startswith("median-"):
                    factor = float(yrange.split("-")[1])
                    ymin = 0
                    ymax = numpy.percentile(group[column],
                                            50.0) * factor
                elif yrange.startswith("max-"):
                    factor = float(yrange.split("-")[1])
                    ymin = 0
                    ymax = numpy.max(group[column]) * factor
                elif yrange.startswith("mean-"):
                    factor = float(yrange.split("-")[1])
                    ymin = 0
                    ymax = numpy.mean(group[column]) * factor
                elif yrange.startswith("balanced-"):
                    ymin = 0
                    # strike balance: truncate at high peaks,
                    # but at least certain multiple of mean
                    factor = float(yrange.split("-")[1])
                    ymax = max(
                        numpy.mean(group[column] * factor),
                        numpy.max(group[column]) * 0.5)
                elif yrange.startswith("threshold-"):
                    ymin = 0
                    # apply threshold, but only if max is below
                    factor = float(yrange.split("-")[1])
                    ymax = min(
                        factor,
                        numpy.max(group[column]))

                # truncate all dots to max value
                group[column][group[column] > ymax] = ymax

            if plot_type == "scatter":
                plots.append(plt.scatter(group[position_column],
                                         group[column],
                                         s=5,
                                         c=colors[subp % len(colors)]))
            elif plot_type == "line":
                plots.append(plt.plot(group[position_column],
                                      group[column],
                                      c=colors[subp % len(colors)])[0])
            labels.append(name)
            if rolling == "average":
                mean = pandas.rolling_mean(group[column],
                                           window_size)
                plt.plot(group[position_column],
                         mean,
                         c=colors[subp % len(colors)])

            subp += 1
            plt.xlim(0, max(data[position_column]))

            plt.tick_params(axis='both', which='major', labelsize=20)
            plt.tick_params(axis='both', which='minor', labelsize=20)

            if yrange is not None:
                ax.set_ylim((ymin, ymax * 1.01))

            ymin, ymax = ax.get_ylim()
            ax = plt.gca()

            if gaps is not None:
                ax = plt.gca()
                ymin, ymax = ax.get_ylim()
                for i, row in gaps.iterrows():
                    ax.add_patch(
                        patches.Rectangle(
                            (row[1], ymin),
                            row[2] - row[1],
                            ymax - ymin,
                            color="black",
                            alpha=0.2))

            last_pos = 0

            for i, row in enumerate(chromosome_offsets.iterrows()):
                label, pos = row
                pos = int(pos)
                ax.add_patch(
                    patches.Rectangle(
                        (last_pos, ymin),
                        pos - last_pos,
                        ymax - ymin,
                        color=chromosome_colors[i],
                        alpha=0.4))
                last_pos = pos

            # set major ticks and turn of labels
            ax.set_xticks(list(chromosome_offsets["size"]))
            ax.xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())

            # label at midpoint
            midpoints = [(a + b) / 2 for a, b in zip(chromosome_offsets["size"][:-1],
                                                     chromosome_offsets["size"][1:])]
            midpoints.append(chromosome_offsets["size"][-1])

            ax.xaxis.set_minor_locator(
                matplotlib.ticker.FixedLocator(midpoints))
            ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(
                chromosome_offsets.index))
            ax.set_ylabel(ylabel)

        if subplots:
            plt.subplot(len(list(grouper.groups.keys())), 1, 1)
            plt.legend(
                plots, labels,
                bbox_to_anchor=(0., 1.0, 1., .102),
                loc=3,
                ncol=3,
                mode="expand",
                borderaxespad=0.,
                prop={'size': 24}
            )
        else:
            plt.legend(plots, labels)

        ax.figure.set_size_inches(40, 10 * len(list(grouper.groups.keys())))
        plt.xlabel(xlabel)

        return ax


class ManhattanPlotWithCopyNumberBars(ManhattanPlot):

    def __init__(self, *args, **kwargs):
        ManhattanPlot.__init__(self, *args, **kwargs)

    def __call__(self, filenames_bed=None, *args, **kwargs):

        c = ManhattanPlot()

        if not filenames_bed:
            return c

        # build increment for chromosomal positions
        chromosome_sizes = get_chromosome_sizes(library="genome")
        chromosome_offsets = chromosome_sizes.cumsum().shift(1).fillna(0)
        chromosome_offsets_dict = dict(
            [(str(x), int(y)) for x, y in chromosome_offsets.iterrows()])
        colors = seaborn.color_palette()
        ax = plt.gca()
        row_index = 1
        y_width = 0.1

        for filename_bed in filenames_bed:
            df = pandas.read_csv(filename_bed,
                                 sep="\t",
                                 header=None,
                                 names=["contig", "start", "end", "copynumber"])
            df["start"] = [y + chromosome_offsets_dict[str(x)]
                           for x, y in zip(df.contig, df.start)]
            df["end"] = [y + chromosome_offsets_dict[str(x)]
                         for x, y in zip(df.contig, df.end)]

            colors = colors[row_index % len(colors)]
            for i, row in df.iterrows():
                alpha = 1.0
                y = row.copynumber
                ax.add_patch(
                    patches.Rectancle(
                        (row.start, y),
                        row.end - row.start,
                        y_width,
                        color=colors,
                        alpha=alpha))
            row_index += 1
        ylimits = ax.get_ylim()
        plt.ylim(0, min(10, ylimits[1]))

        return c


class HeterozygosityFrequenciesPlot(ChromosomePlot):
    position = "position"
    column = "density"
    plot_type = "line"
    ylabel = "proportion of heterozygous variants in window"

    def __init__(self, *args, **kwargs):
        ChromosomePlot.__init__(self, *args, **kwargs)

    def __call__(self, dataframe, path):

        c = ChromosomePlot()
        c.position = self.position
        c.column = self.column
        c.plot_type = self.plot_type
        result = c(dataframe, path)

        df = ROHRegions()(path)

        y_width = 0.02
        y_offset = 1.0 - y_width

        colors = seaborn.color_palette()
        x = 0
        ax = plt.gca()
        for key, group in df.groupby("sample"):
            color = colors[x % len(colors)]
            for i, row in group.iterrows():
                if row.state == 0:
                    alpha = 0.1
                else:
                    alpha = 1.0

                ax.add_patch(
                    patches.Rectangle(
                        (row.start, y_offset),
                        row.end - row.start,
                        y_width,
                        color=color,
                        alpha=alpha))

            y_offset -= y_width
            x += 1
        plt.ylim(0, 1)

        return result


class HeterozygosityFrequenciesManhattanPlot(ManhattanPlot):
    position = "position"
    column = "density"
    plot_type = "line"
    ylabel = "proportion of heterozygous variants in window"

    def __init__(self, *args, **kwargs):
        ManhattanPlot.__init__(self, *args, **kwargs)

    def __call__(self, dataframe, path):

        c = ManhattanPlot()
        c.position = self.position
        c.column = self.column
        c.plot_type = self.plot_type
        result = c(dataframe, path)

        df = ROHRegions()()

        # build increment for chromosomal positions
        chromosome_sizes = get_chromosome_sizes(library="genome")
        chromosome_offsets = chromosome_sizes.cumsum().shift(1).fillna(0)
        chromosome_offsets_dict = dict(
            [(str(x), int(y)) for x, y in chromosome_offsets.iterrows()])

        df["start"] = [y + chromosome_offsets_dict[str(x)]
                       for x, y in zip(df.contig,
                                       df.start)]
        df["end"] = [y + chromosome_offsets_dict[str(x)]
                     for x, y in zip(df.contig,
                                     df.end)]

        y_width = 0.02
        y_offset = 1.0 - y_width

        colors = seaborn.color_palette()
        x = 0
        ax = plt.gca()

        for key, group in df.groupby("sample"):
            color = colors[x % len(colors)]
            for i, row in group.iterrows():
                if row.state == 0:
                    alpha = 0.1
                else:
                    alpha = 1.0

                ax.add_patch(
                    patches.Rectangle(
                        (row.start, y_offset),
                        row.end - row.start,
                        y_width,
                        color=color,
                        alpha=alpha))

            y_offset -= y_width
            x += 1
        plt.ylim(0, 1)

        return result
