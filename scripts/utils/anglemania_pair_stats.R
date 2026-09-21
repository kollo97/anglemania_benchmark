## Shared helpers for reading the per-gene-pair anglemania stats tables
## written by scripts/python/stability_anglemania_stats.py, used by
## scripts/R/stability_lobo_zscore_shift.R and
## scripts/R/stability_lobo_ranking_decomposition.R.
##
## Gene pairs are keyed by a single integer-valued key so 10-17M-row tables
## from different runs can be joined cheaply. The stats tables hold the
## upper triangle only (_select.py's triu_indices), so each unordered pair
## appears once, but the key canonicalises (min/max) anyway in case two runs
## ordered their intersected gene lists differently.

suppressPackageStartupMessages(library(data.table))

## gene names in the simulated data are Gene1..Gene10000
GENE_ID_MULT <- 1e5

read_pair_stats <- function(stats_dir, sample,
                            cols = c("geneA", "geneB", "mean_zscore", "sd_zscore", "rank")) {
    f <- file.path(stats_dir, paste0(sample, "_angl_stats.tsv"))
    stopifnot(file.exists(f))
    dt <- fread(f, select = cols, showProgress = FALSE)
    a <- as.integer(substring(dt$geneA, 5L))     # "Gene1234" -> 1234
    b <- as.integer(substring(dt$geneB, 5L))
    stopifnot(!anyNA(a), !anyNA(b), max(a, b) < GENE_ID_MULT)
    dt[, c("geneA", "geneB") := NULL]
    dt[, pair_key := pmin(a, b) * GENE_ID_MULT + pmax(a, b)]
    stopifnot(!anyDuplicated(dt$pair_key))
    setkey(dt, pair_key)
    dt[]
}

pair_gene_ids <- function(pair_key) {
    list(gene_a = as.integer(pair_key %/% GENE_ID_MULT),
         gene_b = as.integer(pair_key %% GENE_ID_MULT))
}

pair_genes <- function(pair_key) {
    ids <- pair_gene_ids(pair_key)
    list(geneA = paste0("Gene", ids$gene_a), geneB = paste0("Gene", ids$gene_b))
}
