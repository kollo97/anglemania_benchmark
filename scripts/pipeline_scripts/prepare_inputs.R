.libPaths(.Library)

library(optparse)

#------------------------------------------------------------------------------#
# COMMAND LINE ARGUMENTS
#------------------------------------------------------------------------------#
option_list <- list(

    make_option(c("-i", "--infile"), 
                            default=NA,
                            type = "character",
                            help="File location of the h5ad input file"),
    make_option(c("-o", "--outfile"), 
                            default=NA,
                            type = "character",
                            help="File location of the output dir where the integrated data/embeddings should be stored"),
    make_option(c("-b", "--batch_key"), 
                            default="Batch",
                            type = "character",
                            help="batch key specifying a column in the metadata that divides the batches"),
    make_option(c("-d", "--dataset_key"),
                            default=NA,
                            type = "character",
                            help="dataset key specifying a column in the metadata that divides the datasets"),
    make_option(c("-g", "--gene_selection"),
                            default=NA,
                            type = "character",
                            help="gene selection method. One of ['hvg', 'angl', 'full', 'rand']"),
    make_option(c("-m", "--anglemania_mode"),
                            default="cosine",
                            type = "character",
                            help="mode for creating the similarity matrix in anglemania. One of ['cosine', 'spearman', 'diem']"),
    make_option(c("-s", "--zscore_snr_threshold"),
                            default=2.5,
                            type = "double",
                            help="threshold for the signal-to-noise ratio. Default is 2.5"),
    make_option(c("-z", "--zscore_mean_threshold"),
                            default = 2.5,
                            type = "double",
                            help = "threshold for the mean zscore. Default is 2.5")
)
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

suppressPackageStartupMessages({
    library(Seurat)
    library(anndataR)
    library(anglemania)
})

#------------------------------------------------------------------------------#
# CHECK INPUTS
#------------------------------------------------------------------------------#
# COMMENT: TOY
# args$infile <- "/home/akollot/projects/CRC1588/dataset_integration/anglemania_benchmark/output/tcells_hao_pbmc.h5ad"
# args$gene_selection <- "angl"
# args$outfile <- "/home/akollot/projects/CRC1588/dataset_integration/anglemania_benchmark/output/tcells_hao_pbmc.tsv"
# args$batch_key <- "orig.ident"
# check if infile, gene_selection or outfile are NULL
if (is.na(args$infile) | is.na(args$gene_selection) | is.na(args$outfile)) {
    stop(paste0("PARAMETERS: ", paste(names(args[sapply(args, is.na)]), collapse = ", "), " ARE MISSING"))
}

if (!file.exists(args$infile)) {
    stop("No input file found. Please provide a valid h5ad file")
}


if (!(args$gene_selection %in% c("hvg", "angl", "full", "rand"))) {
    stop("Invalid gene selection method. Please provide one of ['hvg', 'angl', 'full', 'rand']")
}

#------------------------------------------------------------------------------#
# MAIN
#------------------------------------------------------------------------------#


# Load data
message(paste0("Processing ", args$infile, " with batch key ", args$batch_key))
se <- anndataR::read_h5ad(args$infile)
se <- se$to_Seurat()
se

if (!(args$batch_key %in% colnames(se[[]]))) {
    stop("Batch key not found in metadata")
}


if (args$gene_selection == "hvg") {
# Create seurat object with only HVGs
    batch_list <- SplitObject(se, split.by = args$batch_key)
    batch_list <- lapply(batch_list, function(tmp_se){
        tmp_se <- NormalizeData(tmp_se, normalization.method = "LogNormalize", scale.factor = 10000)
        return(tmp_se)
        })
    hvgs <- SelectIntegrationFeatures(batch_list, nfeatures = 2000)
    hvgs <- data.frame(hgnc_symbol = hvgs, row.names = names(hvgs))
    write.table(hvgs, args$outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
} else if (args$gene_selection == "angl") {

    # Create seurat object with only anglemania genes
    angl <- create_anglemania_object(se,
                          batch_key = args$batch_key,
                          dataset_key = args$dataset_key,
                          min_cells_per_gene = 1)

    angl <- anglemania(angl,
                            zscore_mean_threshold = args$zscore_mean_threshold,
                            zscore_sn_threshold  = args$zscore_snr_threshold,
                            max_n_genes = 2000,
                            method = args$anglemania_mode) # "cosine", "spearman"

    angl_genes <- get_anglemania_genes(angl)
    angl_genes <- data.frame(hgnc_symbol = angl_genes, row.names = names(angl_genes))
    write.table(angl_genes, args$outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

} else if (args$gene_selection == "full") {

    # Create seurat object with all genes
    full <- data.frame(hgnc_symbol = rownames(se), row.names = rownames(se))
    write.table(full, args$outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

} else if (args$gene_selection == "rand") {

    # Create seurat object with random genes
    set.seed(42)
    rand_genes <- sample(rownames(se), 2000)
    rand_genes <- data.frame(hgnc_symbol = rand_genes, row.names = names(rand_genes))
    write.table(rand_genes, args$outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}
