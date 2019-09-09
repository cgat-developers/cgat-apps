'''GO.py - compute GO enrichment from gene lists
=============================================

:Tags: Python

Usage
-----

The script ``GO.py`` will test for enrichment or depletion of GO
categories within a gene list.

The script uses a hypergeometric test to check if a particular GO
category is enriched in a foreground set with respect to a background
set. Multiple testing is controlled by computing a empirical false
discovery rate using a sampling procedure.

A GO analysis proceeds in three steps:

   1. building gene to GO assignments
   2. create one or more gene lists with foreground and background
   3. run one or more GO analyses for each of the foreground gene lists

This script analyses multiple gene lists in parallel when a matrix of
gene lists is provided. If multiple gene lists are provided, the FDR
is controlled per gene list and not overall. However, my intuition is
that if the number of tests is large the results should be comparable
as if the FDR was controlled globally, though I have no proof for
this.

Building gene to GO assignments
+++++++++++++++++++++++++++++++

The easiest way to obtain a map from gene identifiers to GO assignments
is to down download GO assignments from the ENSEMBL database. The command
below will download go assignments for the human gene set
and save it in the file :file:`gene2go.data`::

   python runGO.py
      --filename-dump=gene2go.data
      --database-host=ensembldb.ensembl.org
      --database-user=anonymous
      --database-name=homo_sapiens_core_54_36p
      --database-port=5306
   > gene2go.log

In order to use GOslim categories, an additional mapping step needs to
be performed.  The sequence of commands is::

    wget http://www.geneontology.org/GO_slims/goslim_goa.obo
    wget http://www.geneontology.org/ontology/gene_ontology.obo
    map2slim -outmap go2goslim.map goslim_goa.obo gene_ontology.obo
    python runGO.py
                --go2goslim
                --filename-ontology=gene_ontology.obo
                --slims=go2goslim.map
                --log=goslim.log
        < gene2go.data > gene2goslim.data

The first two commands obtain GOslim information.  `map2slim
<http://search.cpan.org/~cmungall/go-perl/scripts/map2slim>`_ is part
of Chris Mungall's `go-perl
<http://search.cpan.org/~cmungall/go-perl/>`_ module and the last
command converts the gene-to-GO assignment into gene-to-GOSlim
assignments.

The gene-to-GO mapping can be constructed any other way. It is simply
a table of tab-separated values::

   go_type gene_id go_id   description     evidence
   biol_process    ENSG00000151729 GO:0000002      mitochondrial genome maintenance        NA
   biol_process    ENSG00000025708 GO:0000002      mitochondrial genome maintenance        NA
   biol_process    ENSG00000115204 GO:0000002      mitochondrial genome maintenance        NA
   ...

Building gene lists
+++++++++++++++++++

GO requires a list of genes to test for enrichment. This list is simply
a table with one column of gene identifiers. For example::

   gene_id
   ENSG00000116586
   ENSG00000065809
   ENSG00000164048
   ENSG00000115137
   ENSG00000121210

Alternatively, the gene list can be a multi-column table such as::

   gene_id             dataset1    dataset2
   ENSG00000116586     1           0
   ENSG00000065809     0           0
   ENSG00000164048     1           0
   ENSG00000115137     1           1
   ENSG00000121210     0           1

In this case, enrichment is computed for multiple datasets at once. Make sure
to add the ``%(set)s`` place holder to ``--filename-output-pattern``.

If no background is given, all genes that have GO assignments will constitute
the background.

Statistics
++++++++++

Enrichment is computed using the hypergeometric test.

.. todo::
    * apply filtering
    * more stats
    * more FDR

Running the GO analysis
+++++++++++++++++++++++

The command below runs a GO analysis, computing an FDR using 10.000 samples::

    python runGO.py
        --filename-input=gene2go.data
        --genes-tsv-file=foreground
        --background-tsv-file=background
        --method=sample --sample-size=10000
        --fdr
        --filename-ontology=gene_ontology.obo
        --output-filename-pattern='result/%(set)s.%(go)s.%(section)s'
   > go.log

The output will be stored in the directory :file:`result` and output
files will be created according to the pattern
``<set>.<go>.<section>``. ``<set>`` is the gene set that is analysed,
``<go>`` is one of ``biol_process``, ``mol_function`` and
``cell_location``.  ``<section>`` denotes the file contents. Files
output are:

+------------+----------------------------------------------+
|``section`` | contents                                     |
+------------+----------------------------------------------+
|samples     |sampling statistics                           |
+------------+----------------------------------------------+
|overall     |table with full results                       |
+------------+----------------------------------------------+
|results     |table with only the significant results       |
+------------+----------------------------------------------+
|parameters  |input and sampling parameters                 |
+------------+----------------------------------------------+
|fg          |assigments for genes in the foreground set    |
+------------+----------------------------------------------+

Other options
+++++++++++++

The script can accept other ontologies than just GO ontologies.

Command line options
--------------------

'''
import sys
import collections
import cgatcore.database as database
import cgatcore.experiment as E
import cgatcore.iotools as iotools

import cgat.GO as GO


def main(argv=None):

    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--version", action='version', version='runGO v1.0')

    parser.add_argument(
        "-s", "--species", dest="species", type=str,
        help="species to use.")

    parser.add_argument(
        "-i", "--slims", dest="filename_slims", type=str,
        help="filename with GO SLIM categories ")

    parser.add_argument(
        "-g", "--genes-tsv-file", dest="filename_genes", type=str,
        help="filename with genes to analyse ")

    parser.add_argument(
        "-b", "--background-tsv-file", dest="filename_background",
        type=str,
        help="filename with background genes to analyse ")

    parser.add_argument(
        "-m", "--min-counts", dest="minimum_counts",
        type=int,
        help="minimum count - ignore all categories that have "
        "fewer than # number of genes")

    parser.add_argument(
        "-o", "--sort-order", dest="sort_order", type=str,
        choices=("fdr", "pvalue", "ratio"),
        help="output sort order.")

    parser.add_argument(
        "--ontology", dest="ontology", type=str,
        action="append",
        help="go ontologies to analyze. Ontologies are tested "
        "separately.")

    parser.add_argument(
        "-t", "--threshold", dest="threshold", type=float,
        help="significance threshold [>1.0 = all ]. If --fdr is set, this "
        "refers to the fdr, otherwise it is a cutoff for p-values.")

    parser.add_argument(
        "--filename-dump", dest="filename_dump", type=str,
        help="dump GO category assignments into a flatfile ")

    parser.add_argument(
        "--gene2name-map-tsv-file", dest="filename_gene2name", type=str,
        help="optional filename mapping gene identifiers to gene names ")

    parser.add_argument(
        "--filename-ontology", dest="filename_ontology", type=str,
        help="filename with ontology in OBO format.")

    parser.add_argument(
        "--filename-input", dest="filename_input", type=str,
        help="read GO category assignments from a flatfile ")

    parser.add_argument(
        "--sample-size", dest="sample", type=int,
        help="do sampling (with # samples).")

    parser.add_argument(
        "--filename-output-pattern", "--output-filename-pattern",
        dest="output_filename_pattern", type=str,
        help="pattern with output filename pattern "
        "(should contain: %(go)s and %(section)s )")

    parser.add_argument(
        "--fdr", dest="fdr", action="store_true",
        help="calculate and filter by FDR.")

    parser.add_argument(
        "--go2goslim", dest="go2goslim", action="store_true",
        help="convert go assignments in STDIN to goslim assignments and "
        "write to STDOUT.")

    parser.add_argument(
        "--gene-pattern", dest="gene_pattern", type=str,
        help="pattern to transform identifiers to GO gene names ")

    parser.add_argument(
        "--filename-map-slims", dest="filename_map_slims", type=str,
        help="write mapping between GO categories and GOSlims ")

    parser.add_argument(
        "--get-genes", dest="get_genes", type=str,
        help="list all genes in the with a certain GOID.")

    parser.add_argument(
        "--strict", dest="strict", action="store_true",
        help="require all genes in foreground to be part of background. "
        "If not set, genes in foreground will be added to the background ")

    parser.add_argument(
        "-q", "--fdr-method", dest="qvalue_method", type=str,
        choices=("empirical", "storey", "BH"),
        help="method to perform multiple testing correction by controlling "
        "the fdr .")

    parser.add_argument(
        "--pairwise", dest="compute_pairwise", action="store_true",
        help="compute pairwise enrichment for multiple gene lists. ")

    # parser.add_argument( "--fdr-lambda", dest="qvalue_lambda", type=float,
    #                   help="fdr computation: lambda [default=%default]."  )

    # parser.add_argument( "--qvalue-pi0-method", dest="qvalue_pi0_method", type=str,
    #                    choices = ("smoother", "bootstrap" ),
    # help="fdr computation: method for estimating pi0 [default=%default]."  )

    parser.set_defaults(species=None,
                        filename_genes="-",
                        filename_background=None,
                        filename_slims=None,
                        minimum_counts=0,
                        ontology=[],
                        filename_dump=None,
                        sample=0,
                        fdr=False,
                        output_filename_pattern=None,
                        threshold=0.05,
                        filename_map_slims=None,
                        gene_pattern=None,
                        sort_order="ratio",
                        get_genes=None,
                        strict=False,
                        qvalue_method="empirical",
                        pairs_min_observed_counts=3,
                        compute_pairwise=False,
                        filename_gene2name=None
                        )

    (args) = E.start(parser, add_database_options=True)

    if args.go2goslim:
        GO.convertGo2Goslim(args)
        E.stop()
        sys.exit(0)

    if args.fdr and args.sample == 0:
        E.warn("fdr will be computed without sampling")

    #############################################################
    # dump GO
    if args.filename_dump:
        # set default orthologies to GO
        if not args.ontology:
            args.ontology = [
                "biol_process", "mol_function", "cell_location"]

        E.info("dumping GO categories to %s" % (args.filename_dump))

        dbhandle = database.connect(url=args.database_url)

        outfile = iotools.open_file(args.filename_dump, "w", create_dir=True)
        GO.DumpGOFromDatabase(outfile,
                              dbhandle,
                              args)
        outfile.close()
        E.stop()
        sys.exit(0)

    #############################################################
    # read GO categories from file
    if args.filename_input:
        E.info("reading association of categories and genes from %s" %
               (args.filename_input))
        infile = iotools.open_file(args.filename_input)
        gene2gos, go2infos = GO.ReadGene2GOFromFile(infile)
        infile.close()

    if args.filename_gene2name:
        E.info("reading gene identifier to gene name mapping from %s" %
               args.filename_gene2name)
        infile = iotools.open_file(args.filename_gene2name)
        gene2name = iotools.read_map(infile, has_header=True)
        infile.close()
        E.info("read %i gene names for %i gene identifiers" %
               (len(set(gene2name.values())),
                len(gene2name)))
    else:
        # use identity mapping
        gene2name = dict([(x, x) for x in list(gene2gos.keys())])

    #############################################################
    # read GO ontology from file
    if args.filename_ontology:
        E.info("reading ontology from %s" % (args.filename_ontology))

        infile = iotools.open_file(args.filename_ontology)
        ontology = GO.readOntology(infile)
        infile.close()

        def _g():
            return collections.defaultdict(GO.GOInfo)
        go2infos = collections.defaultdict(_g)

        # substitute go2infos
        for go in list(ontology.values()):
            go2infos[go.mNameSpace][go.mId] = GO.GOInfo(
                go.mId,
                go_type=go.mNameSpace,
                description=go.mName)

    #############################################################
    # get foreground gene list
    input_foreground, genelists = GO.ReadGeneLists(
        args.filename_genes,
        gene_pattern=args.gene_pattern)

    E.info("read %i genes for forground in %i gene lists" %
           (len(input_foreground), len(genelists)))

    #############################################################
    # get background
    if args.filename_background:

        # nick - bug fix: background is the first tuple element from
        # ReadGeneLists
        input_background = GO.ReadGeneLists(
            args.filename_background,
            gene_pattern=args.gene_pattern)[0]
        E.info("read %i genes for background" % len(input_background))
    else:
        input_background = None

    #############################################################
    # sort out which ontologies to test
    if not args.ontology:
        if args.filename_input:
            args.ontology = list(gene2gos.keys())

    E.info("found %i ontologies: %s" %
           (len(args.ontology), args.ontology))

    summary = []
    summary.append("\t".join((
        "genelist",
        "ontology",
        "significant",
        "threshold",
        "ngenes",
        "ncategories",
        "nmaps",
        "nforegound",
        "nforeground_mapped",
        "nbackground",
        "nbackground_mapped",
        "nsample_counts",
        "nbackground_counts",
        "psample_assignments",
        "pbackground_assignments",
        "messages")) + "\n")

    #############################################################
    # get go categories for genes
    for test_ontology in sorted(args.ontology):

        # store results for aggregate output of multiple gene lists
        all_results = []
        all_significant_results = []
        all_genelists_with_results = []

        E.info("working on ontology %s" % test_ontology)
        #############################################################
        # get/read association of GO categories to genes
        if args.filename_input:
            gene2go, go2info = gene2gos[test_ontology], go2infos[test_ontology]
        else:
            E.info("reading data from database ...")

            dbhandle.Connect(args)
            gene2go, go2info = GO.ReadGene2GOFromDatabase(
                dbhandle,
                test_ontology,
                args.database, args.species)

            E.info("finished")

        if len(go2info) == 0:
            E.warn(
                "could not find information for terms - "
                "could be mismatch between ontologies")

        ngenes, ncategories, nmaps, counts_per_category = GO.CountGO(gene2go)
        E.info("assignments found: %i genes mapped to %i categories "
               "(%i maps)" %
               (ngenes, ncategories, nmaps))

        if args.minimum_counts > 0:
            to_remove = set(
                [x for x, y in counts_per_category.items()
                 if y < args.minimum_counts])
            E.info("removing %i categories with less than %i genes" %
                   (len(to_remove), args.minimum_counts))
            GO.removeCategories(gene2go, to_remove)

            ngenes, ncategories, nmaps, counts_per_category = \
                GO.CountGO(gene2go)
            E.info("assignments after filtering: %i genes mapped "
                   "to %i categories (%i maps)" % (
                       ngenes, ncategories, nmaps))

        for genelist_name, foreground in sorted(genelists.items()):

            msgs = []
            E.info("processing %s with %i genes" %
                   (genelist_name, len(foreground)))
            ##################################################################
            ##################################################################
            ##################################################################
            # build background - reconcile with foreground
            ##################################################################
            if input_background is None:
                background = list(gene2go.keys())
            else:
                background = list(input_background)

            # nick - bug-fix backgorund included the foreground in a tuple.
            # background is the first tuple element
            missing = foreground.difference(set(background))

            if args.strict:
                assert len(missing) == 0, \
                    "%i genes in foreground but not in background: %s" % (
                        len(missing), str(missing))
            else:
                if len(missing) != 0:
                    E.warn("%i genes in foreground that are not in "
                           "background - added to background of %i" %
                           (len(missing), len(background)))

                background.extend(missing)

            E.info("(unfiltered) foreground=%i, background=%i" %
                   (len(foreground), len(background)))

            # sort foreground and background, important for reproducibility
            # under random seed
            foreground = sorted(foreground)
            background = sorted(background)

            #############################################################
            # sanity checks:
            # are all of the foreground genes in the dataset
            # missing = set(genes).difference( set(gene2go.keys()) )
            # assert len(missing) == 0, "%i genes in foreground set without GO annotation: %s" % (len(missing), str(missing))

            #############################################################
            # read GO slims and map GO categories to GO slim categories
            if args.filename_slims:
                go_slims = GO.GetGOSlims(
                    iotools.open_file(args.filename_slims, "r"))

                if args.loglevel >= 1:
                    v = set()
                    for x in list(go_slims.values()):
                        for xx in x:
                            v.add(xx)
                    args.stdlog.write(
                        "# read go slims from %s: go=%i, slim=%i\n" %
                        (args.filename_slims,
                         len(go_slims),
                         len(v)))

                if args.filename_map_slims:
                    if args.filename_map_slims == "-":
                        outfile = args.stdout
                    else:
                        outfile = iotools.open_file(
                            args.filename_map_slims, "w")

                    outfile.write("GO\tGOSlim\n")
                    for go, go_slim in sorted(list(go_slims.items())):
                        outfile.write("%s\t%s\n" % (go, go_slim))

                    if outfile != args.stdout:
                        outfile.close()

                gene2go = GO.MapGO2Slims(gene2go, go_slims, ontology=ontology)

                if args.loglevel >= 1:
                    ngenes, ncategories, nmaps, counts_per_category = \
                        GO.CountGO(gene2go)
                    args.stdlog.write(
                        "# after go slim filtering: %i genes mapped to "
                        "%i categories (%i maps)\n" % (
                            ngenes, ncategories, nmaps))

            #############################################################
            # Just dump out the gene list
            if args.get_genes:
                fg, bg, ng = [], [], []

                for gene, vv in list(gene2go.items()):
                    for v in vv:
                        if v.mGOId == args.get_genes:
                            if gene in genes:
                                fg.append(gene)
                            elif gene in background:
                                bg.append(gene)
                            else:
                                ng.append(gene)

                # skip to next GO class
                if not (bg or ng):
                    continue

                args.stdout.write(
                    "# genes in GO category %s\n" % args.get_genes)
                args.stdout.write("gene\tset\n")
                for x in sorted(fg):
                    args.stdout.write("%s\t%s\n" % ("fg", x))
                for x in sorted(bg):
                    args.stdout.write("%s\t%s\n" % ("bg", x))
                for x in sorted(ng):
                    args.stdout.write("%s\t%s\n" % ("ng", x))

                E.info("nfg=%i, nbg=%i, nng=%i" % (len(fg), len(bg), len(ng)))

                E.stop()
                sys.exit(0)

            #############################################################
            outfile = GO.getFileName(args,
                                     go=test_ontology,
                                     section='foreground',
                                     set=genelist_name)

            outfile.write("gene_id\n%s\n" % ("\n".join(sorted(foreground))))
            if args.output_filename_pattern:
                outfile.close()

            outfile = GO.getFileName(args,
                                     go=test_ontology,
                                     section='background',
                                     set=genelist_name)

            # Jethro bug fix - see section 'build background' for assignment
            outfile.write("gene_id\n%s\n" % ("\n".join(sorted(background))))
            if args.output_filename_pattern:
                outfile.close()

            #############################################################
            # do the analysis
            go_results = GO.AnalyseGO(gene2go, foreground, background)

            if len(go_results.mSampleGenes) == 0:
                E.warn("%s: no genes with GO categories - analysis aborted" %
                       genelist_name)
                continue

            pairs = list(go_results.mResults.items())

            #############################################################
            # calculate fdr for each hypothesis
            if args.fdr:
                fdrs, samples, method = GO.computeFDRs(go_results,
                                                       foreground,
                                                       background,
                                                       args,
                                                       test_ontology,
                                                       gene2go,
                                                       go2info)
                for x, v in enumerate(pairs):
                    v[1].mQValue = fdrs[v[0]][0]
            else:
                fdrs, samples, method = {}, {}, None

            msgs.append("fdr=%s" % method)

            if args.sort_order == "fdr":
                pairs.sort(key=lambda x: x[1].mQValue)
            elif args.sort_order == "ratio":
                pairs.sort(key=lambda x: x[1].mRatio)
            elif args.sort_order == "pvalue":
                pairs.sort(key=lambda x: x[1].mPValue)

            #############################################################
            #############################################################
            #############################################################
            # output the full result
            outfile = GO.getFileName(args,
                                     go=test_ontology,
                                     section='overall',
                                     set=genelist_name)

            GO.outputResults(
                outfile, pairs, go2info, args, fdrs=fdrs, samples=samples)

            if args.output_filename_pattern:
                outfile.close()

            #############################################################
            #############################################################
            #############################################################
            # filter significant results and output
            filtered_pairs = GO.selectSignificantResults(pairs, fdrs, args)

            nselected = len(filtered_pairs)
            nselected_up = len([x for x in filtered_pairs if x[1].mRatio > 1])
            nselected_down = len(
                [x for x in filtered_pairs if x[1].mRatio < 1])

            assert nselected_up + nselected_down == nselected

            outfile = GO.getFileName(args,
                                     go=test_ontology,
                                     section='results',
                                     set=genelist_name)

            GO.outputResults(outfile,
                             filtered_pairs,
                             go2info,
                             args,
                             fdrs=fdrs,
                             samples=samples)

            if args.output_filename_pattern:
                outfile.close()

            #############################################################
            #############################################################
            #############################################################
            # save results for multi-gene-list analysis
            all_results.append(pairs)
            all_significant_results.append(filtered_pairs)
            all_genelists_with_results.append(genelist_name)

            #############################################################
            #############################################################
            #############################################################
            # output parameters
            ngenes, ncategories, nmaps, counts_per_category = \
                GO.CountGO(gene2go)

            outfile = GO.getFileName(args,
                                     go=test_ontology,
                                     section='parameters',
                                     set=genelist_name)

            nbackground = len(background)
            if nbackground == 0:
                nbackground = len(go_results.mBackgroundGenes)

            outfile.write(
                "# input go mappings for gene list '%s' and category '%s'\n" %
                (genelist_name, test_ontology))
            outfile.write("parameter\tvalue\tdescription\n")
            outfile.write("mapped_genes\t%i\tmapped genes\n" % ngenes)
            outfile.write(
                "mapped_categories\t%i\tmapped categories\n" % ncategories)
            outfile.write("mappings\t%i\tmappings\n" % nmaps)
            outfile.write("genes_in_fg\t%i\tgenes in foreground\n" %
                          len(foreground))
            outfile.write(
                "genes_in_fg_with_assignment\t%i\tgenes in foreground with GO assignments\n" %
                (len(go_results.mSampleGenes)))
            outfile.write(
                "genes_in_bg\t%i\tinput background\n" % nbackground)
            outfile.write(
                "genes_in_bg_with_assignment\t%i\tgenes in background with GO assignments\n" % (
                    len(go_results.mBackgroundGenes)))
            outfile.write(
                "associations_in_fg\t%i\tassociations in sample\n" %
                go_results.mSampleCountsTotal)
            outfile.write(
                "associations_in_bg\t%i\tassociations in background\n" %
                go_results.mBackgroundCountsTotal)
            outfile.write(
                "percent_genes_in_fg_with_association\t%s\tpercent genes in sample with GO assignments\n" % (
                    iotools.pretty_percent(len(go_results.mSampleGenes),
                                           len(foreground), "%5.2f")))
            outfile.write(
                "percent_genes_in_bg_with_associations\t%s\tpercent genes background with GO assignments\n" % (
                    iotools.pretty_percent(len(go_results.mBackgroundGenes),
                                           nbackground, "%5.2f")))
            outfile.write(
                "significant\t%i\tsignificant results reported\n" % nselected)
            outfile.write(
                "significant_up\t%i\tsignificant up-regulated results reported\n" % nselected_up)
            outfile.write(
                "significant_down\t%i\tsignificant up-regulated results reported\n" % nselected_down)
            outfile.write(
                "threshold\t%6.4f\tsignificance threshold\n" % args.threshold)

            if args.output_filename_pattern:
                outfile.close()

            summary.append("\t".join(map(str, (
                genelist_name,
                test_ontology,
                nselected,
                args.threshold,
                ngenes,
                ncategories,
                nmaps,
                len(foreground),
                len(go_results.mSampleGenes),
                nbackground,
                len(go_results.mBackgroundGenes),
                go_results.mSampleCountsTotal,
                go_results.mBackgroundCountsTotal,
                iotools.pretty_percent(
                    len(go_results.mSampleGenes), len(foreground), "%5.2f"),
                iotools.pretty_percent(
                    len(go_results.mBackgroundGenes), nbackground, "%5.2f"),
                ",".join(msgs)))) + "\n")

            #############################################################
            #############################################################
            #############################################################
            # output the fg patterns
            outfile = GO.getFileName(args,
                                     go=test_ontology,
                                     section='withgenes',
                                     set=genelist_name)

            GO.outputResults(outfile, pairs, go2info, args,
                             fdrs=fdrs,
                             samples=samples,
                             gene2go=gene2go,
                             foreground=foreground,
                             gene2name=gene2name)

            if args.output_filename_pattern:
                outfile.close()

        if len(genelists) > 1:

            ###################################################################
            # output various summary files
            # significant results
            GO.outputMultipleGeneListResults(all_significant_results,
                                             all_genelists_with_results,
                                             test_ontology,
                                             go2info,
                                             args,
                                             section='significant')

            # all results
            GO.outputMultipleGeneListResults(all_results,
                                             all_genelists_with_results,
                                             test_ontology,
                                             go2info,
                                             args,
                                             section='all')

            if args.compute_pairwise:
                GO.pairwiseGOEnrichment(all_results,
                                        all_genelists_with_results,
                                        test_ontology,
                                        go2info,
                                        args)

    outfile_summary = args.stdout
    outfile_summary.write("".join(summary))

    E.stop()


if __name__ == "__main__":
    sys.exit(main())
