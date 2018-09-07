#' filtering single cell data based on QC metrics
#'
#' Example usage:
#' 
#' cgat sc-counts2counts --counts-filename=featurecounts.tsv --phenotypes-filename=phenodata.tsv --factor=group,mouse_id,collection_date,slice_depth,slice_number,pipette_visual,timepoint > filtered_counts.tsv
#'
#' `feature_counts.tsv` is a table (tab-separated) of ngenes x ncells,
#' that is the genes are in rows and the columns are cells.
#'
#' `phenodata.tsv` is a table (tab-separated) of ncells x nfeatures,
#' that is rows are cells and features are in columns. The table should contain
#' a column called `sample_id` that will match the columns in the table
#' `feature_counts.tsv`.
#'
#' Features can then be selected in the `--factor` option to be
#' plotted.
#' 
#' -> todo: parameterize detection of ERCC (pattern?)
#' -> todo: parameterize definition of mitochondrial genes - currently hardcoded for mouse.

## conda dependencies: bioconductor-scater r-cairo

suppressMessages(library(futile.logger))
suppressMessages(library(getopt))
suppressMessages(library(Cairo))
suppressMessages(library(scater))

source(file.path(Sys.getenv("R_ROOT"), "experiment.R"))


start_plot <- function(section, height = 6, width = 10, type = "png") {
    file = get_output_filename(paste0(section, ".", type))
    Cairo(file = file,
          type = type,
          width = width,
          height = height,
          units="in",
          dpi = 300,
          bg = "white",
          pointsize = 12)
    opar <- par(lwd=0.5)
}

end_plot <- function() {
    dev.off()
}


plot_qc <- function(data_sceset, section, factors, opt) {

    flog.info("histogram of total counts")
    start_plot(paste("qc_total_counts_histogram", section, sep = "-"))
    hist(data_sceset$total_counts / 1e3,
         xlab="Library size (thousands)",
         main="",
         breaks=20,
         col="dodgerblue3",
         border="white",
         cex.lab=1.25,
         ylab="Number of cells")
    end_plot()

    flog.info("histogram of total features")
    start_plot(paste("qc_total_features_histogram.png", section, sep = "-"))
    hist(data_sceset$total_features,
         xlab="Number of expressed genes",
         main="",
         breaks=20,
         col="dodgerblue3",
         border="white",
         cex.lab=1.25,
         ylab="Number of cells")
    end_plot()

    flog.info("histogram of mitochondrial proportions")
    start_plot(paste("percent_mitochondria_histogram", section, sep = "-"))
    hist(data_sceset$pct_counts_Mt,
         xlab="Proportion of reads that map to mitochondrial genes (%)",
         main="", breaks=20, col="dodgerblue3", border="white", cex.lab=1.25, ylab="Number of cells")
    end_plot()

    flog.info("scatter plot of mitochondrial reads with +/- mad to show cutoff used (red line) for cells")
    start_plot(paste("percent_mitochondria_histogram_mad", section, sep = "-"))
    plot(data_sceset$pct_counts_Mt, ylim=c(-10, 30), xlab="Number of cells", main="",
         pch=20, col="dodgerblue3",
         ylab="Proportion of reads that map to mitochondrial genes (%)",
         cex.lab=1.25, axes=FALSE)
    axis(2, cex.axis=1)
    axis(1, cex.axis=1)
    x1 = median(data_sceset$pct_counts_Mt)
    y1 = mad(data_sceset$pct_counts_Mt)
    z1 = x1 - (opt$percent_counts_mito_nmads * y1)
    z2 = x1 + (opt$percent_counts_mito_nmads * y1)
    abline(z2,0, col = "red", lwd=0.75)
    text(data_sceset$pct_counts_Mt, colnames(data_sceset), cex=0.4, pos=4, col="orange")
    abline(z1,0, col = "snow4", lwd=0.75)
    end_plot()
    
    flog.info("histogram of ERCC proportions")
    start_plot(paste("percent_ERCC_histogram", section, sep = "-"))
    hist(data_sceset$pct_counts_ERCC,
         xlab="Proportion of reads that map to ERCCs (%)",
         main="",
         breaks=20,
         col="dodgerblue3",
         border="white",
         cex.lab=1.25,
         ylab="Number of cells")
    end_plot()
    
    flog.info("scatter plot of ERCC reads with +/- mad=3 to show cutoff used (red line) for cells")
    start_plot(paste("percent_ERCC_scatter_mad", section, sep = "-"))
    x1 = median(data_sceset$pct_counts_ERCC)
    y1 = mad(data_sceset$pct_counts_ERCC)
    z1 = x1 - (opt$percent_counts_ercc_nmads * y1)
    z2 = x1 + (opt$percent_counts_ercc_nmads * y1)
    plot(data_sceset$pct_counts_ERCC,
         ylim=c(-1, 30),
         xlab="Number of cells",
         main="",
         pch=20, col="dodgerblue3",
         ylab="Proportion of reads that map to ERCCs (%)",
         cex.lab=1.25,
         axes=FALSE)
    axis(2, cex.axis=1);axis(1, cex.axis=1);abline(z1,0, col = "snow4", lwd=0.75)
    abline(z2,0, col = "red", lwd=0.75)
    text(data_sceset$pct_counts_ERCC,
         colnames(data_sceset),
         cex=0.4,
         pos=4,
         col="orange")
    end_plot()

    flog.info(paste("plotting highest expressed genes"))
    start_plot(paste("qc_highest_expression", section, sep = "-"))
    print(plotQC(data_sceset, type = "highest-expression"))
    end_plot()

    flog.info("scatter plot of log10 of total features with +/- mad to show cutoff used (red line) for cells")
    start_plot(paste("qc_total_features_scatter_mad", section, sep = "-"))
    plot(log10(data_sceset$total_features),
         ylim=c(3.25, 4.5),
         xlab="Number of cells",
         main="",
         pch=20,
         col="dodgerblue3",
         ylab="Number of expressed genes (log10 scale)",
         cex.lab=1.25,
         axes=FALSE);
    axis(2, cex.axis=1)
    axis(1, cex.axis=1)
    x1 = median(log10(data_sceset$total_features))
    y1 = mad(log10(data_sceset$total_features))
    z1 = x1 - (opt$total_features_nmads * y1)
    z2 = x1 + (opt$total_features_nmads * y1)
    abline(z1, 0, col = "red", lwd=0.75)
    abline(z2, 0, col = "snow4", lwd=0.75)
    end_plot()

    for (factor in factors) {
        flog.info(paste("creating plots for factor", factor))

        for (feature in c("total_counts", "total_features",
                          "pct_counts_Mt", "pct_counts_ERCC")) {
            
            df <- data.frame(data_sceset[[feature]],
                             data_sceset[[factor]])
            colnames(df) <- c(feature, factor)

            flog.info(paste("histogram of", feature))
            start_plot(paste("qc", feature, "histogram", section, factor, sep = "-"))
            print(ggplot(df,
                         aes_string(x = feature)) + geom_histogram(bins=50) +
                  facet_wrap(factor) +
                  labs(x=feature, y="Number of cells"))
            end_plot()

            flog.info(paste("kde of", feature))
            start_plot(paste("qc", feature, "kde", section, factor, sep = "-"))
            print(ggplot(df,
                         aes_string(x = feature, colour = factor)) + geom_density() +
                  labs(x=feature, y="Number of cells"))
            end_plot()
        }
        
        flog.info(paste("plotting percentage counts in mitochondrial genes into"))
        start_plot(paste("qc_percentage_counts_mitochondrial_genes", section, factor, sep = "-"))
        print(scater::plotPhenoData(
                          data_sceset,
                          aes_string(
                              x = "total_features",
                              y = "pct_counts_Mt",
                              colour = factor
                              )
                      ))
        end_plot()

        flog.info(paste("plotting percentage counts in ERCC genes"))
        start_plot(paste("qc_percentage_counts_ERCC_genes", section, factor, sep = "-"))
        print(scater::plotPhenoData(
                          data_sceset,
                          aes_string(
                          x = "total_features",
                          y = "pct_counts_ERCC",
                          colour = factor
                      )
                      ))
        end_plot()

        flog.info(paste("plotting PCA2"))
        start_plot(paste("qc_pca2_counts", section, factor, sep = "-"))
        print(plotPCA(data_sceset,
                      ntop = 500,
                      exprs_values = "counts",
                      ncomponents = 2,
                      colour_by = factor,
                      size_by = "total_features"))
        end_plot()

        flog.info(paste("plotting PCA4"))
        start_plot(paste("qc_pca4_counts", section, factor, sep = "-"))
        print(plotPCA(data_sceset,
                      ntop = 500,
                      exprs_values = "counts",
                      ncomponents = 4,
                      colour_by = factor))
        end_plot()
        
        ## endog_genes <- !rowData(reads.qc)$is_feature_control
        ## note: do not use logcounts_raw, but logcounts for downstream
        ## start_plot(paste("qc_find_pcs.png"))
        ## print(plotQC(
        ##     data_sceset,
        ##     type = "find-pcs",
        ##     exprs_values = "logcounts_raw",
        ##     variable = "total_features"
        ## ))
        ## end_plot()
        
    }
}


plot_average_counts <- function(ave.counts, opt) {
    
    start_plot("counts_frequency_cutoff_raw")
    hist(ave.counts,
         breaks = 40000,
         xlim=c(0, 200),
         main="",
         col="dodgerblue2",
         cex.lab=1.25,
         xlab="Average count")
    end_plot()
    
    start_plot("counts_frequency_cutoff")
    hist(log10(ave.counts),
         breaks=250,
         main="",
         col="dodgerblue2",
         cex.lab=1.25,
         xlab=expression(Log[10]~"average count"))
    abline(v=log10(opt$average_counts_threshold), col="red", lwd=0.5)
    end_plot()
}

run <- function(opt) {

    options(stringsAsFactors = FALSE)
    set.seed(1234567)

    flog.info(paste("reading counts data from", normalizePath(opt$counts_filename)))

    counts_data <-read.table(opt$counts_filename, header = TRUE, row.names = 1, sep = "\t")
    flog.info(paste("read counts data", paste(dim(counts_data), collapse = ",")))

    flog.info(paste("reading phenotype data from", normalizePath(opt$phenotypes_filename)))
    annotation_data <- read.table(opt$phenotypes_filename, sep = "\t", header = TRUE)
    flog.info(paste("read phenotype data", paste(dim(annotation_data), collapse = ",")))
    
    flog.info(paste("read annotation data", paste(dim(annotation_data), collapse = "-")))
    annotation_data$sample_id <- gsub("-", ".", annotation_data$sample_id)
    row.names(annotation_data) <- annotation_data$sample_id

    flog.info("building SingleCellExperiment data set")
    all_sceset <- SingleCellExperiment(assays = list(counts = as.matrix(counts_data)), colData = annotation_data)
    flog.info(paste("built SingleCellExperiment data set", paste(dim(all_sceset), collapse=",")))

    flog.info("removing genes not expressed in any cell")
    keep_feature <- rowSums(counts(all_sceset) > 0) > 0
    flog.info(paste("keeping", sum(keep_feature), "genes"))
    all_sceset <- all_sceset[keep_feature, ]
    
    ercc <- rownames(all_sceset)[grep("ERCC", rownames(all_sceset))]
    mt <- c("ENSMUSG00000064336","ENSMUSG00000064337","ENSMUSG00000064338",
            "ENSMUSG00000064339","ENSMUSG00000064340","ENSMUSG00000064341",
            "ENSMUSG00000064342","ENSMUSG00000064343","ENSMUSG00000064344",
            "ENSMUSG00000064345","ENSMUSG00000064346","ENSMUSG00000064347",
            "ENSMUSG00000064348","ENSMUSG00000064349","ENSMUSG00000064350",
            "ENSMUSG00000064351","ENSMUSG00000064352","ENSMUSG00000064353",
            "ENSMUSG00000064354","ENSMUSG00000064355","ENSMUSG00000064356",
            "ENSMUSG00000064357","ENSMUSG00000064358","ENSMUSG00000064359",
            "ENSMUSG00000064360","ENSMUSG00000064361","ENSMUSG00000064363",
            "ENSMUSG00000064364","ENSMUSG00000064365","ENSMUSG00000064366",
            "ENSMUSG00000064367","ENSMUSG00000064368","ENSMUSG00000064369",
            "ENSMUSG00000064370","ENSMUSG00000064371","ENSMUSG00000064372",
            "ENSMUSG00000065947")

    is.spike <- (rownames(all_sceset) %in% ercc)
    isSpike(all_sceset, type="ERCC") <- is.spike

    is.mito <- (rownames(all_sceset) %in% mt)
    isSpike(all_sceset, type="Mt") <- is.mito
    
    flog.info(paste("marking", sum(is.spike), "spike-in genes"))
    flog.info(paste("marking", sum(is.mito), "mitochondrial genes"))

    flog.info("calculating SCRAN QC metrics using spike-in and mitochondrial controls")
    all_sceset <- scater::calculateQCMetrics(all_sceset,
                                              feature_controls = list(ERCC = is.spike, Mt = is.mito))
    all_sceset <- scater::arrange(all_sceset, group)

    ## print(colnames(colData(all_sceset)))
    ## print(colnames(rowData(all_sceset)))

    flog.info("QC plotting of unfiltered data")
    plot_qc(all_sceset, section = "unfiltered", factors = opt$factors, opt = opt)

    ## CELL FILTERING#####
    flog.info("stage 1 filtering: filter by cell")
    flog.info(paste("filter on library size smaller than nmads of", opt$library_size_nmads))
    libsize.drop <- isOutlier(all_sceset$total_counts,
                              nmads =  opt$library_size_nmads,
                              type = "lower",
                              log = TRUE)

    flog.info(paste("filter on total features smaller than nmads of", opt$total_features_nmads))
    feature.drop <- isOutlier(all_sceset$total_features,
                              nmads = opt$total_features_nmads,
                              type = "lower",
                              log = TRUE)

    flog.info(paste("filter on the percent of counts that map to mitochondrial genes higher than nmads of",
              opt$percent_counts_mito_nmads))
    mito.drop <- isOutlier(all_sceset$pct_counts_Mt,
                           nmads = opt$percent_counts_mito_nmads,
                           type="higher")

    flog.info(paste("filter on the percent of counts that map to ERCCS higher than nmads of",
                    opt$percent_counts_ercc_nmads))
    spike.drop <- isOutlier(all_sceset$pct_counts_ERCC,
                            nmads = opt$percent_counts_ercc_nmads,
                            type = "higher")

    filtered_sceset <- all_sceset[,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
    flog.info(paste("removed by library size", sum(libsize.drop)))
    flog.info(paste("removed by total features", sum(feature.drop)))
    flog.info(paste("removed by mitochondrial expression", sum(mito.drop)))
    flog.info(paste("removed by spike-in", sum(spike.drop)))
    flog.info(paste("after stage 1 filtering", paste(dim(filtered_sceset), collapse = ",")))

    flog.info("stage 2 filtering: filtering average counts")
    flog.info("filter low abundance genes")
    ave.counts <- calcAverage(filtered_sceset)
    plot_average_counts(ave.counts, opt = opt)
    keep <- ave.counts >= opt$average_counts_threshold
    filtered_sceset<- filtered_sceset[keep,]
    flog.info(paste("after stage 2 filtering", paste(dim(filtered_sceset), collapse = ",")))
    
    flog.info(paste("stage 3 filtering: filter out genes expressed in fewer than",
                    opt$num_cells_threshold, "cells"))
    numcells <- nexprs(filtered_sceset, byrow=TRUE)
    alt.keep <- numcells >= opt$num_cells_threshold
    filtered_sceset <- filtered_sceset[alt.keep,]
    flog.info(paste("after stage 3 filtering", paste(dim(filtered_sceset), collapse = ",")))

    flog.info("QC plotting of filtered data")
    plot_qc(filtered_sceset, section = "filtered", factors = opt$factors, opt = opt)
    
    ## alt method: hard cut-off of at least a total of x number of reads in at least x number of cells
    ## keep_feature <- rowSums(counts(reads) > 100) > 5
    ## reads <- reads[keep_feature, ]
    ## head(reads)

    ## alt method: outlier detection based on PCA using scater
    ## new_sceset<- scater::plotPCA(
    ##     old_sceset,
    ##     size_by = "total_features",
    ##     pca_data_input = "pdata",
    ##     detect_outliers = TRUE,
    ##     return_SCE = TRUE
    ##     )
    
    ## Set up gene lengths for RPKM
    flog.info("outputting filtered counts data")
    write.table(counts(filtered_sceset),
                file = "",
                sep = "\t",
                quote = FALSE,
                row.names = TRUE,
                col.names = TRUE)
}

main <- function() {

    option_list <- list(
        make_option(
            "--counts-filename",
            dest = "counts_filename",
            type = "character",
            default = "featurecounts.tsv",
            help = paste("filename with input data of counts")
        ),
        make_option(
            "--phenotypes-filename",
            dest = "phenotypes_filename",
            type = "character",
            default = "phenodata.tsv",
            help = paste("filename with phenotype data")
        ),
        make_option(
            "--factor",
            dest = "factors",
            type = "character",
            default = "group,collection_date",
            help = paste("factors to colour QC plots by.")
        ),
        make_option(
            "--num-cells-threshold",
            dest = "num_cells_threshold",
            type = "numeric",
            default = 5,
            help = paste("remove low genes in at least # cells")
        ),
        make_option(
            "--average-counts-threshold",
            dest = "average_counts_threshold",
            type = "numeric",
            default = 1,
            help = paste("remove genes with average count below #")
        ),
        make_option(
            "--library-size-nmads",
            dest = "library_size_nmads",
            type = "numeric",
            default = 3,
            help = paste("remove cells with library size outside # number of median-",
                         "absolute-deviations away from median")
        ),
        make_option(
            "--total-features-nmads",
            dest = "total_features_nmads",
            type = "numeric",
            default = 3,
            help = paste("remove cells with total features outside # number of median-",
                         "absolute-deviations away from median")
        ),
        make_option(
            "--percent-counts-mito-nmads",
            dest = "percent_counts_mito_nmads",
            type = "numeric",
            default = 5,
            help = paste(
                "remove cells with percent counts in mitochondrial genes outside # number of median-",
                "absolute-deviations away from median")
        ),
        make_option(
            "--percent-counts-ercc-nmads",
            dest = "percent_counts_ercc_nmads",
            type = "numeric",
            default = 4,
            help = paste(
                "remove cells with percent counts in ERCC genes outside # number of median-",
                "absolute-deviations away from median")
        )
    )
    opt <- experiment_start(option_list = option_list,
                            description = description)

    if (!is.null(opt$factors)) {
        opt$factors = unlist(strsplit(opt$factors, ","))
    }
    run(opt)
    
    experiment_stop()
}

main()
