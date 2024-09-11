.libPaths(.Library)
suppressPackageStartupMessages({
    library(Seurat)
    library(schard)
    library(optparse)
    library(reticulate)
    use_condaenv(condaenv = "/home/akollot/miniforge3/envs/angl_BM/bin/python", conda = "/home/akollot/miniforge3/envs/angl_BM/bin/python")
})
devtools::load_all("/home/akollot/projects/CRC1588/dataset_integration/anglemania/")

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
                            help="batch key specifying a column in the metadat that divides the batches"),
    make_option(c("g", "--gene_selection"),
                            default=NA,
                            type = "character",
                            help="gene selection method. One of ['hvg', 'angl', 'full', 'rand']")

)
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)


#------------------------------------------------------------------------------#
# CHECK INPUTS
#------------------------------------------------------------------------------#
# TOY
# args$infile <- "/local/Projects/AAkalin_Neuroblastoma/Results/simulate_single_cell/h5ad/batch.facLoc0.1_de.facLoc0.1_ngroup2.h5ad"
# args$gene_selection <- "hvg"
# args$outfile <- "test_hvg.tsv"
# args$batch_key <- "Batch"
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
se <- anndataR::read_h5ad(args$infile, to = "Seurat", mode = "r")


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
    angl <- create_anglem(se,
                          batch_key = args$batch_key,
                          min_cells_per_gene = 1)

    angl <- big_anglemanise(angl,
                            zscore_mean_threshold = 2.5,
                            zscore_sn_threshold  = 2.5,
                            max_n_genes = 2000)

    angl_genes <- extract_integration_genes(angl)
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
