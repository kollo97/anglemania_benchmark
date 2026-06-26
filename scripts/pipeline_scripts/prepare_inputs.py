import argparse
import sys

import anndata as ad
import numpy as np
import pandas as pd

#------------------------------------------------------------------------------#
# COMMAND LINE ARGUMENTS
#------------------------------------------------------------------------------#
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", required=True,
                     help="File location of the h5ad input file")
parser.add_argument("-o", "--outfile", required=True,
                     help="File location of the output TSV with the selected genes")
parser.add_argument("-b", "--batch_key", default="Batch",
                     help="batch key specifying a column in the metadata that divides the batches")
parser.add_argument("-d", "--dataset_key", default=None,
                     help="dataset key specifying a column in the metadata that divides the datasets")
parser.add_argument("-g", "--gene_selection", required=True,
                     choices=["hvg", "angl", "full", "rand", "topbtvr", "bottombtvr"],
                     help="gene selection method")
parser.add_argument("-m", "--anglemania_mode", default="cosine",
                     choices=["cosine", "spearman", "phi_s"],
                     help="mode for creating the similarity matrix in anglemania")
parser.add_argument("-n", "--n_genes", default=2000, type=int,
                     help="number of genes to select")
parser.add_argument("-p", "--permutation_function", default="sample",
                     choices=["sample", "permute_nonzero"],
                     help="function to permute the counts for the null distribution")
parser.add_argument("--normalization_mode", default="classical",
                     choices=["classical", "pflog1ppf"],
                     help="normalization applied inside anglemania before computing angles: "
                          "'classical' (CP10K + log1p) or 'pflog1ppf' (shifted-CLR)")

NORMALIZATION_METHODS = {
    "classical": "divide_by_total_counts",
    "pflog1ppf": "pflog1ppf",
}


def select_hvg(adata, batch_key, n_genes):
    import scanpy as sc

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat",
        n_top_genes=n_genes,
        batch_key=batch_key,
    )
    return adata.var_names[adata.var["highly_variable"]].tolist()


def select_angl(adata, batch_key, dataset_key, n_genes, method, permutation_function, normalization_mode):
    from pyanglemania import preprocessing as pp

    _move_to_gpu(adata)
    pp.anglemania(
        adata,
        batch_key=batch_key,
        dataset_key=dataset_key,
        max_n_genes=n_genes,
        min_cells_per_gene=1,
        method=method,
        permutation_function=permutation_function,
        do_normalize=True,
        normalization_method=NORMALIZATION_METHODS[normalization_mode],
    )
    return list(adata.uns["anglemania"]["anglemania_genes"])


def _move_to_gpu(adata):
    import cupy as cp
    from cupyx.scipy import sparse as csp
    from scipy import sparse as sp

    try:
        has_gpu = cp.cuda.runtime.getDeviceCount() > 0
    except Exception:
        has_gpu = False
    if not has_gpu:
        raise RuntimeError(
            "No CUDA device reachable; the 'angl' gene selection step requires a GPU."
        )
    adata.X = csp.csr_matrix(adata.X) if sp.issparse(adata.X) else cp.asarray(adata.X)


def select_full(adata):
    return adata.var_names.tolist()


def select_rand(adata, n_genes):
    rng = np.random.default_rng(42)
    n_select = min(n_genes, adata.n_vars)
    return rng.choice(adata.var_names, size=n_select, replace=False).tolist()


def select_btvr(adata, gene_selection, n_genes):
    # Select genes by ground-truth BTVR (Between-to-Total Variance Ratio).
    # Requires splatter-simulated h5ad files that store DEFacGroup* (biological
    # variance factors) and BatchFacBatch* (technical variance factors) in var.
    var_df = adata.var
    bio_cols = [c for c in var_df.columns if c.startswith("DEFacGroup")]
    tech_cols = [c for c in var_df.columns if c.startswith("BatchFacBatch")]
    if not bio_cols or not tech_cols:
        raise ValueError(
            "topbtvr/bottombtvr requires DEFacGroup* and BatchFacBatch* columns "
            f"in var (splatter-simulated data only). Found: {list(var_df.columns)}"
        )
    bio_var = var_df[bio_cols].var(axis=1)
    tech_var = var_df[tech_cols].var(axis=1)
    btvr_raw = np.log((bio_var + 1e-15) / (tech_var + 1e-15))
    btvr = (btvr_raw - btvr_raw.min()) / (btvr_raw.max() - btvr_raw.min())
    ascending = gene_selection == "bottombtvr"
    n_select = min(n_genes, var_df.shape[0])
    return btvr.sort_values(ascending=ascending).index[:n_select].tolist()


def main(args):
    print(f"Processing {args.infile} with batch key {args.batch_key}", file=sys.stderr)
    adata = ad.read_h5ad(args.infile)

    if args.gene_selection == "hvg":
        genes = select_hvg(adata, args.batch_key, args.n_genes)
    elif args.gene_selection == "angl":
        genes = select_angl(
            adata, args.batch_key, args.dataset_key, args.n_genes,
            args.anglemania_mode, args.permutation_function, args.normalization_mode,
        )
    elif args.gene_selection == "full":
        genes = select_full(adata)
    elif args.gene_selection == "rand":
        genes = select_rand(adata, args.n_genes)
    else:
        genes = select_btvr(adata, args.gene_selection, args.n_genes)

    print(f"Selected {len(genes)} genes for integration.", file=sys.stderr)
    pd.DataFrame({"hgnc_symbol": genes}).to_csv(args.outfile, sep="\t", index=False)


if __name__ == "__main__":
    main(parser.parse_args())
