library(optparse)

#------------------------------------------------------------------------------#
# COMMAND LINE ARGUMENTS
#------------------------------------------------------------------------------#
option_list <- list(

    make_option(c("-i", "--original_h5ad"),
                            default = NA,
                            type = "character",
                            help = "File location of the original h5ad input file"),
    make_option(c("-d", "--integrated_h5ad"),
                            default = NA,
                            type = "character",
                            help = "File location of the integrated h5ad file"),
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
                            help = "gene selection method. One of ['hvg', 'angl', 'full', 'rand']"),
    make_option(c("-b", "--batch_key"),
                            default = "batch",
                            type = "character",
                            help = "batch key specifying a column in the metadata that divides the batches"),
    make_option(c("-l", "--label_key"),
                            default = "CellType",
                            type = "character",
                            help = "label key specifying a column in the metadata that specifies the cell type information")
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
#' Calculate the Difference in Local Density Factor metric for an integrated
#' dataset
#'
#' @param original_h5ad SingleCellExperiment object containing the expression
#' matrix
#' @param integrated_h5ad SingleCellExperiment object containing the integrated
#' dataset
#' @param batch_key Key for batch information in metadata SCE object
#' @param label_key Key for cell type/grouping information in metadata of
#' the SCE object
#'
#' @returns The ldfDiff metric score
calculate_ldfDiff <- function(
    original_sce,
    integrated_sce,
    batch_key,
    label_key
) {
    set.seed(1)
    n_dim <- ncol(SingleCellExperiment::reducedDim(integrated_sce, "emb"))

    message("Calculating batch PCAs...")
    batches <- sort(unique(original_sce[[batch_key]]))
    batch_objects <- lapply(batches, function(.batch) {
        message("Batch '", .batch, "'...")
        batch_sce <- original_sce[, original_sce[[batch_key]] == .batch]
        batch_sce <- scuttle::logNormCounts(batch_sce)
        batch_sce <- scater::runPCA(
            batch_sce,
            exprs_values = "logcounts",
            ncomponents = n_dim,
            name = "pca"
        )
    })
    names(batch_objects) <- batches

    # Use k=75 unless one of the batches is smaller than that
    # (should only happen for the test dataset)
    k <- min(sapply(batch_objects, ncol)) - 1
    if (k < 75) {
        warning(
            "k was set to ",
            k,
            " because one of the batches has fewer than 76 cells"
        )
    } else {
        k <- 75
    }
    message("Calculating cell ldfDiff scores...")
    integrated_sce <- CellMixS::ldfDiff(
        sce_pre_list = batch_objects,
        sce_combined = integrated_sce,
        group = batch_key,
        k = k,
        dim_red = "pca",
        dim_combined = "emb",
        n_dim = n_dim
    )

    message("Calculating final ldfDiff score...")
    # Scores closer to 0 are better so use absolute values
    scores <- abs(SummarizedExperiment::colData(integrated_sce)$diff_ldf)
    if (any(is.na(scores))) {
        message(
            "Warning: Ignoring ",
            sum(is.na(scores)),
            " cells with NA ldfDiff scores"
        )
        scores <- scores[!is.na(scores)]
    }
    # Score are unbounded so we set any scores greater than 1 to 1
    if (any(scores > 1)) {
        warning(
            sum(scores > 1),
            " absolute cell scores (",
            round(sum(scores > 1) / ncol(integrated_sce) * 100, 2),
            "%) are greater than 1. Setting these to 1."
        )
        scores[scores > 1] <- 1
    }
    # Take the mean of the scores and subtract it from 1 so higher is better
    score <- 1 - mean(scores)

    return(score)
}


#------------------------------------------------------------------------------#
# MAIN
#------------------------------------------------------------------------------#
original_sce <- anndataR::read_h5ad(
    args$original_h5ad,
    to = "SingleCellExperiment"
)
integrated_sce <- anndataR::read_h5ad(
    args$integrated_h5ad,
    to = "SingleCellExperiment"
)
ldfdiff_score <- calculate_ldfDiff(
    original_sce = original_sce,
    integrated_sce = integrated_sce,
    batch_key = args$batch_key,
    label_key = args$label_key
)


metrics <- data.frame(
    ldfDiff = ldfdiff_score,
    sample = args$sample,
    integration_method = args$integration_method,
    gene_selection = args$gene_selection
)
message("Writing metrics to ", args$outfile)
write.table(
    metrics,
    file = args$outfile,
    row.names = FALSE,
    quote = FALSE
)