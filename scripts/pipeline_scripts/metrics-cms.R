library(optparse)

#------------------------------------------------------------------------------#
# COMMAND LINE ARGUMENTS
#------------------------------------------------------------------------------#

option_list <- list(
    make_option(c("-a", "--original_h5ad"),
                            default = NA,
                            type = "character",
                            help = "File location of the original h5ad input file"),
    make_option(c("-e", "--embedding_tsv"),
                            default = NA,
                            type = "character",
                            help = "File location of the embedding TSV file (cell barcodes as row index)"),
    make_option(c("-o", "--outfile"),
                            default = NA,
                            type = "character",
                            help = "File location of the output file where the metrics"),
    make_option(c("-s", "--sample"),
                            default = NA,
                            type = "character",
                            help = "sample name"),
    make_option(c("-m", "--integration_method"),
                            default = NA,
                            type = "character",
                            help = "integration method. One of ['harmony','scvi','scanvi']"),
    make_option(c("-g", "--gene_selection"),
                            default = NA,
                            type = "character",
                            help = "gene selection method. One of ['hvg', 'angl', 'full', 'rand', 'topbtvr', 'bottombtvr']"),
    make_option(c("-b", "--batch_key"),
                            default = "batch",
                            type = "character",
                            help = "batch key specifying a column in the metadata that divides the batches")
)
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

suppressPackageStartupMessages({
    library(anndataR)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(CellMixS)
    library(scater)
    library(scuttle)
    library(checkmate)
})

#------------------------------------------------------------------------------#
# FUNCTION
#------------------------------------------------------------------------------#
#' Calculate CMS (Cell Mixing Score)
#' 
#' @description
#' This function calculates the Cell Mixing Score (CMS) of an integrated dataset
#' based on the CellMixS package. The CMS is a measure of cell mixing, which
#' reflects how well individual cells represent their respective batches.
#' 
#' @param integrated_sce SingleCellExperiment object containing the integrated
#' dataset
#' @param batch_key Key for batch information in metadata SCE object
#'
#' @return The CMS score, which is a value between 0 and 1, where 1 indicates
#' complete mixing and 0 indicates no mixing.
calculate_cms <- function(
    integrated_sce,
    batch_key
) {
    set.seed(1)
    k <- min(200, ncol(integrated_sce) - 1L)
    integrated_sce <- suppressMessages(CellMixS::cms(
        integrated_sce,
        k = k,
        group = batch_key,
        dim_red = "emb",
        n_dim = ncol(SingleCellExperiment::reducedDim(integrated_sce, "emb"))
    ))

    cms_scores <- SummarizedExperiment::colData(integrated_sce)$cms
    score <- 1 - mean(cms_scores < 0.1, na.rm = TRUE)
    return(score)
}

#------------------------------------------------------------------------------#
# MAIN
#------------------------------------------------------------------------------#
# TEMPORARY--------------
# args <- list()
# args$original_h5ad <- "/path/to/original.h5ad"
# args$embedding_tsv <- "/path/to/embedding.tsv"
# args$batch_key <- "Batch"
# args$sample <- "test"
# args$gene_selection <- "test"
# args$integration_method <- "test"
# args$outfile <- sprintf("%s%s%s_%s_cms.tsv", args$sample, args$integration_method, args$sample, args$gene_selection)
# ----------------------
original_sce <- anndataR::read_h5ad(
    args$original_h5ad,
    to = "SingleCellExperiment"
)
emb_df <- read.table(args$embedding_tsv, sep = "\t", header = TRUE, row.names = 1)
integrated_sce <- original_sce[, rownames(emb_df)]
SingleCellExperiment::reducedDim(integrated_sce, "emb") <- as.matrix(emb_df)
# CellMixS requires the batch column to be a factor; anndataR reads it as character
integrated_sce[[args$batch_key]] <- factor(integrated_sce[[args$batch_key]])

cms <- calculate_cms(
    integrated_sce = integrated_sce,
    batch_key = args$batch_key
)

metrics <- data.frame(
    cms = cms,
    sample = args$sample,
    integration_method = args$integration_method,
    gene_selection = args$gene_selection
)
write.table(
    metrics,
    file = args$outfile,
    row.names = FALSE,
    quote = FALSE
)
