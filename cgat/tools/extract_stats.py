'''
extract_stats.py - extract and process tables from csvDB
========================================================

Purpose
-------

Extract tables from sqlite databases and process

Usage
-----

.. Example use case

Example::

   python extract_stats.py

Type::

   python extract_stats.py --help

for command line help.

Command line options
--------------------

tasks
+++++

`extract_table` - extract a table from a database
                  containing relevant information

`get_coverage` - calculate gene/transcript model
                 coverage stats - only works on single cell data
                 with filename format <seqRun>_<plate>_<well>_<mappper>

`aggregate` - aggregate together multiple stats tables,
              and select relevant measures

.. todo::
    Needs to be refactored to use sqlalchemy or Database instead of
    explicitely wrapping database frontends.


'''

import sys
import pandas
import numpy
import re
import cgatcore.experiment as E
import cgatcore.database as database


def getTableFromDb(database_url, table):
    '''
    Get a table from a database with pandas
    '''

    dbhandle = database.connect(url=database_url)
    df = pandas.read_sql("SELECT * FROM {}".format(table), con=dbhandle)
    df.index = df["track"]
    df.drop(labels="track", inplace=True, axis=1)

    return df


def cleanStatsTable(stats_file):
    '''
    Take in a table containing aggregated stats
    and clean by removing duplicate columns
    '''
    # , mangle_dupe_cols=False)
    # AH: disabled, because "ValueError: Setting mangle_dupe_cols=False is not supported yet"
    df = pandas.read_table(stats_file, sep="\t", header=0,
                           index_col=None)

    # drop duplicates is case sensitive, convert all to
    # same case - SQL is not case sensitive so will throw
    # a hissy fit for same column names in different cases
    df.columns = [cx.lower() for cx in df.columns]
    df = df.T.drop_duplicates().T
    df.index = df["track"]
    return df


def extractTranscriptCounts(con, table):
    '''
    Extract transcript model counts for a
    given sample

    Arguments
    ---------
    con: sqlite.connection
      An SQLite connection

    table: string
      the table to extract the transcript counts
      from.

    Returns
    -------
    coverages: pandas.Core.Series
    '''

    statement = '''
    SELECT coverage_sense_pcovered
    FROM %(table)s
    WHERE coverage_sense_nval > 0;
    ''' % locals()

    coverages = pandas.read_sql(statement, con)
    coverages = coverages.loc[:, "coverage_sense_pcovered"]
    return coverages


def summariseOverBins(coverages, bins):
    '''
    Summarise model coverages over a set of bins

    Argumnets
    ---------
    coverages: pandas.Core.Series
      coverages over gene/transcripts

    bins: list
      values corresponding to percentage bins

    Returns
    -------
    freqs: numpy.array
      frequency array of coverages over percentiles
    '''

    freqs = numpy.zeros(shape=len(bins), dtype=numpy.float64)
    for i in range(len(bins)):
        if i == 0:
            hits = coverages <= bins[i]
        else:
            hits = (coverages <= bins[i]) & (coverages > bins[i-1])

        freqs[i] = len(coverages[hits])

    return freqs


def getModelCoverage(database_url, table_regex, model_type="transcript"):
    '''
    Compute transcript model coverage stats

    Arguments
    ---------
    database_url: string
      database containing transcript counts

    table_regex: string
      regular expression for transcript count table

    model_type: string
      calculate coverages over either transcripts or
      genes.  Default is gene models

    Returns
    -------
    coverage_df: Pandas.Core.DataFrame
      model coverage stats summarised for each cell
    '''

    # need to regex for all the tables, one for each sample
    # fetch_all returns a list of tuples
    dbhandle = Database.connect(database_url)
    cc = dbhandle.execute("SELECT name FROM sqlite_master WHERE type='table';")

    tab_reg = re.compile(table_regex)
    table_list = [tx[0] for tx in cc.fetchall() if re.search(tab_reg, tx[0])]

    # pull out counts for each cell and compute coverages
    bins = range(0, 101)
    cov_dict = {}
    for tab in table_list:
        covs = extractTranscriptCounts(dbhandle, tab)
        freq_array = summariseOverBins(covs, bins)
        cov_dict[tab] = freq_array

    coverage_df = pandas.DataFrame(cov_dict).T
    # create a regex group to remove superfluous characters
    # from the track names
    ix_re = re.compile("_(?P<run>\d+)_(?P<plate>\d+)_(?P<well>\d+)_(?P<mapper>\S+)_transcript_counts")
    re_matches = [re.match(ix_re, ix) for ix in coverage_df.index]
    indx = ["%s_%s-%s.%s" % rm.group(1, 2, 3, 4) for rm in re_matches]
    coverage_df.index = indx
    coverage_df.columns = ["Bin%i" % bx for bx in coverage_df.columns]
    return coverage_df


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("--task", dest="task", type=str,
                        choices=["extract_table", "get_coverage",
                                 "clean_table"],
                        help="task to perform")

    parser.add_argument("-t", "--table-name", dest="table", type=str,
                        help="table in SQLite DB to extract")

    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv, add_database_options=True)

    if args.task == "extract_table":
        out_df = getTableFromDb(args.database_url, args.table)

    elif args.task == "get_coverage":
        out_df = getModelCoverage(args.database_url,
                                  table_regex="(\S+)_transcript_counts")

    elif args.task == "clean_table":
        infile = argv[-1]
        out_df = cleanStatsTable(infile)

    out_df.to_csv(args.stdout,
                  sep="\t", index_label="track")

    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
