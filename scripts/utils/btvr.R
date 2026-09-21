## Ground-truth BTVR (Between-to-Total Variance Ratio) straight from a
## splatter-simulated .h5ad.
##
## Same definition used by scripts/pipeline_scripts/prepare_inputs.py
## (select_btvr, for the topbtvr/bottombtvr gene selections) and by
## scripts/notebooks/btvr_sanity_check_visualization.Rmd, factored out here so
## the stability scripts can reuse it:
##   BTVR = rescale( log( var(DEFacGroup*) / var(BatchFacBatch*) ) )
## i.e. per gene, the log ratio of its across-cell-type DE factor spread
## (biological signal) to its across-batch factor spread (technical signal),
## min-max rescaled to [0, 1]. High BTVR = gene whose variation is mostly
## biological, which is what a gene selection ought to prefer.

suppressPackageStartupMessages({
    library(rhdf5)
})

## var(.) of an h5ad, as a data frame with gene names as rownames.
get_var_from_h5ad <- function(h5ad_path) {
    stopifnot(file.exists(h5ad_path))
    v <- rhdf5::h5read(h5ad_path, "/var")
    rhdf5::h5closeAll()
    index_col <- intersect(c("_index", "X_index", "index"), names(v))
    stopifnot(length(index_col) >= 1)
    genes <- as.character(v[[index_col[1]]])
    v[index_col] <- NULL
    df <- as.data.frame(lapply(v, as.vector), stringsAsFactors = FALSE)
    rownames(df) <- genes
    df
}

## Per-gene BTVR from an already-loaded var data frame.
get_btvr <- function(var_df) {
    bio_cols <- grep("^DEFacGroup", colnames(var_df), value = TRUE)
    tech_cols <- grep("^BatchFacBatch", colnames(var_df), value = TRUE)
    if (!length(bio_cols) || !length(tech_cols)) {
        stop("BTVR needs DEFacGroup*/BatchFacBatch* columns in var ",
             "(splatter-simulated data only). Found: ",
             paste(colnames(var_df), collapse = ", "))
    }
    raw_bio_var <- apply(as.matrix(var_df[, bio_cols, drop = FALSE]), 1, var)
    raw_tech_var <- apply(as.matrix(var_df[, tech_cols, drop = FALSE]), 1, var)
    btvr_raw <- log((raw_bio_var + 1e-15) / (raw_tech_var + 1e-15))
    rng <- range(btvr_raw[is.finite(btvr_raw)])
    data.frame(
        gene = rownames(var_df),
        btvr = (btvr_raw - rng[1]) / (rng[2] - rng[1]),
        btvr_raw = btvr_raw,
        raw_bio_var = raw_bio_var,
        raw_tech_var = raw_tech_var,
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}

get_btvr_from_h5ad <- function(h5ad_path) get_btvr(get_var_from_h5ad(h5ad_path))

## BTVR indexed by the integer gene id used in pair_key (Gene1234 -> 1234),
## for O(1) lookup of ~10M gene pairs. Returns a numeric vector where
## position i holds the BTVR of "Gene<i>" (NA where absent).
btvr_by_gene_id <- function(btvr_df, column = "btvr") {
    ids <- suppressWarnings(as.integer(sub("^Gene", "", btvr_df$gene)))
    stopifnot(!anyNA(ids))
    out <- rep(NA_real_, max(ids))
    out[ids] <- btvr_df[[column]]
    out
}
