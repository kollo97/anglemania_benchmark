library(Seurat)
library(schard)
library(optparse)
library(SeuratDisk)
#------------------------------------------------------------------------------#
# COMMAND LINE ARGUMENTS
#------------------------------------------------------------------------------#
option_list <- list(

    make_option(c("-i", "--infile"), 
                            default="/data/akalin/akollot/anglemania_benchmark/data/pbmcsca/adata_hvg.h5ad",
                            type = "character",
                            help="File location of the h5ad input file"),
    make_option(c("-o", "--outdir"), 
                            default="/data/akalin/akollot/anglemania_benchmark/results/pbmcsca/hvg/",
                            type = "character",
                            help="File location of the output dir where the integrated data/embeddings should be stored"),
    make_option(c("-b", "--batch_key"), 
                            default="batch",
                            type = "character",
                            help="batch key specifying a column in the metadat that divides the batches"),
    make_option(c("-l", "--label_key"), 
                            default="CellType", 
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

if (!dir.exists(args$outdir)) {
    stop("No output directory found. Please provide a valid output directory")
}


#------------------------------------------------------------------------------#
# MAIN
#------------------------------------------------------------------------------#
message(paste0("Processing ", args$infile, " with batch key ", args$batch_key, " and label key ", args$label_key))

# Load data
se <- schard::h5ad2seurat(args$infile)
test <- LayerData(se, assay = "RNA", layer = "counts")

if (!(args$batch_key %in% colnames(se[[]]))) {
    stop("Batch key not found in metadata")
}

batch_list <- SplitObject(se, split.by = args$batch_key)

anchors <- FindIntegrationAnchors(batch_list, anchor.features = rownames(se)) # the object put into the scripts are already filtered down to only contain HVGs or only the anglemania genes, so we take all of those

integrated_se <- IntegrateData(anchors)

Seurat::DefaultAssay(integrated_se) <- "integrated"
integrated_se <- integrated_se %>% NormalizeData() %>% ScaleData()
integrated_se <- Seurat::RunPCA(integrated_se, reduction.name = "X_emb")

message(paste0("Saving ", args$outdir, "/seurat.H5Seurat and converting to", args$outdir, "/seurat.h5ad"))
SeuratDisk::SaveH5Seurat(integrated_se, paste0(args$outdir, "seurat.H5Seurat"), overwrite=TRUE)
SeuratDisk::Convert(source = paste0(args$outdir, "seurat.H5Seurat"),
                    dest = paste0(args$outdir, "seurat.h5ad"),
                    assay="integrated",
                    overwrite=TRUE)

# Rscript seurat_integration.R --infile /data/akalin/akollot/anglemania_benchmark/data/pbmcsca/adata_hvg.h5ad --outdir /data/akalin/akollot/anglemania_benchmark/results/pbmcsca/hvg/ --batch_key batch --label_key CellType