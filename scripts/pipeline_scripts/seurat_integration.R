suppressPackageStartupMessages({
    library(Seurat)
    library(optparse)
    library(dplyr)
    library(anndataR)
})
#------------------------------------------------------------------------------#
# COMMAND LINE ARGUMENTS
#------------------------------------------------------------------------------#
option_list <- list(

    make_option(c("-i", "--infile"), 
                            default="/local/Projects/AAkalin_Neuroblastoma/Results/simulate_single_cell/h5ad/batch.facLoc0.1_de.facLoc0.1_ngroup2.h5ad",
                            type = "character",
                            help="File location of the h5ad input file"),
    make_option(c("f", "--feature_subset"),
                            default="hvg.tsv",
                            type = "character",
                            help="File location of the tsv file with the subset of genes to use for integration"),
    make_option(c("-o", "--outfile"), 
                            default="/data/akalin/akollot/anglemania_benchmark/results/pbmcsca/hvg/test.h5ad",
                            type = "character",
                            help="File location of the output file where the integrated data/embeddings should be stored"),
    make_option(c("-b", "--batch_key"), 
                            default="Batch",
                            type = "character",
                            help="batch key specifying a column in the metadat that divides the batches"),
    make_option(c("-l", "--label_key"), 
                            default="Group", 
                            type = "character",
                            help="label key specifying a column in the metadat that specifies the cell type information")

)
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

#------------------------------------------------------------------------------#
# CHECK INPUTS
#------------------------------------------------------------------------------#

if (!file.exists(args$infile)) {
    stop("No input file found. Please provide a valid h5ad file")
}


#------------------------------------------------------------------------------#
# MAIN
#------------------------------------------------------------------------------#
message(paste0("Processing ", args$infile, " with batch key ", args$batch_key, " and label key ", args$label_key))

se <- anndataR::read_h5ad(args$infile, to = "Seurat", mode = "r")
se <- NormalizeData(se, normalization.method = "LogNormalize", scale.factor = 10000)

features <- read.csv(args$feature_subset, sep = "\t", header = TRUE)
features <- as.character(features$hgnc_symbol)
se <- subset(se, features = features)

if (!(args$batch_key %in% colnames(se[[]]))) {
    stop("Batch key not found in metadata")
}

batch_list <- SplitObject(se, split.by = args$batch_key)

anchors <- FindIntegrationAnchors(batch_list, anchor.features = rownames(se)) # the object put into the scripts are already filtered down to only contain HVGs or only the anglemania genes, so we take all of those

integrated_se <- IntegrateData(anchors)

Seurat::DefaultAssay(integrated_se) <- "integrated" 
integrated_se[["integrated"]]$counts <- integrated_se[["RNA"]]$counts
# for writing the h5ad I need to set the default assay to integrated
# then when writing the h5ad, it will return an anndata object where X = raw_counts and data = integrated_counts

message(paste0("Saving ", args$outfile))
anndataR::write_h5ad(se, args$outfile)


# Rscript pipeline_scripts/seurat_integration.R --infile /data/akalin/akollot/anglemania_benchmark/data/pbmcsca/adata_hvg.h5ad --outdir /data/akalin/akollot/anglemania_benchmark/results/pbmcsca/hvg/ --batch_key batch --label_key CellType