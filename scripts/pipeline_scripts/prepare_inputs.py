import argparse
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd

# gene_level_scores lives next to this script; make it importable regardless of
# the cwd Snakemake happens to invoke us from.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gene_level_scores import VARIANT_SPECS as ANGLGENE_VARIANTS  # noqa: E402

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
                     choices=["hvg", "angl", "hvg_rfilter", "hvg_rweight", "hvg_intersect",
                              "full", "rand", "topbtvr", "bottombtvr",
                              *ANGLGENE_VARIANTS],
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
parser.add_argument("--n_bins", default=20, type=int,
                     help="anglgene only: equal-frequency mean-expression bins the "
                          "per-gene signal/R z-scores are computed within")
parser.add_argument("--weight_signal", default=0.5, type=float,
                     help="anglgene only: weight of the binned signal z-score in the score")
parser.add_argument("--weight_R", default=0.5, type=float,
                     help="anglgene only: weight of the binned R z-score in the score")
parser.add_argument("--r_keep_frac", default=0.5, type=float,
                     help="hvg_rfilter only: fraction of the scored genes to keep "
                          "(the best by binned R z-score) before running HVG on them")
parser.add_argument("--r_hvg_weight", default=1.0, type=float,
                     help="hvg_rweight only: weight of the binned R z-score added to "
                          "scanpy's dispersions_norm (both are within-expression-bin "
                          "z-scores, so 1.0 weights them equally)")
parser.add_argument("--loess_span", default=0.5, type=float,
                     help="anglgene_* only: LOESS span for the two residual axes")
parser.add_argument("--moderate_n_bins", default=10, type=int,
                     help="anglgene_*_mod only: mean-expression bins per gene; the "
                          "moderated-S prior trend is fitted on the unordered pair "
                          "of the two genes' bins")
parser.add_argument("--score_table", default=None,
                     help="anglgene_* only: reuse a per-gene score table written by "
                          "make_genelevel_variant_lists.py instead of recomputing the "
                          "cross-batch M/S (which is identical across these arms)")

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


def select_anglgene(adata, tag, batch_key, dataset_key, n_genes, method,
                    permutation_function, normalization_mode, n_bins, weight_signal,
                    weight_R, loess_span, moderate_n_bins, score_table=None):
    """Per-gene selection -- see pipeline_scripts/gene_level_scores.py.

    Same cross-batch z-score statistics anglemania itself computes, reduced to a
    per-gene score (row sums of M^2 and S^2) instead of ranking gene pairs and
    taking genes on first appearance. Which column the arm ranks on, in which
    direction, and whether S is moderated first, all come from ``VARIANT_SPECS``.

    Every ``anglgene_*`` arm derives from the *same* M/S, so recomputing it per
    arm is pure waste; ``--score_table`` lets the caller hand in one already
    computed by ``make_genelevel_variant_lists.py``.
    """
    from gene_level_scores import VARIANT_SPECS, score_genes, top_genes

    column, ascending, needs_moderation = VARIANT_SPECS[tag]

    if score_table is not None:
        df = pd.read_csv(score_table, sep="\t")
        print(f"Reusing score table {score_table} ({len(df)} genes)", file=sys.stderr)
    else:
        df = score_genes(
            adata,
            batch_key=batch_key,
            dataset_key=dataset_key,
            method=method,
            normalization_mode=normalization_mode,
            permutation_function=permutation_function,
            min_cells_per_gene=1,
            gpu=True,
            n_bins=n_bins,
            weight_signal=weight_signal,
            weight_R=weight_R,
            loess_span=loess_span,
            moderate=needs_moderation,
            moderate_n_bins=moderate_n_bins,
        )
    return top_genes(df, n_genes, column=column, ascending=ascending)


def select_hvg_r(adata, batch_key, dataset_key, n_genes, method, permutation_function,
                 normalization_mode, n_bins, mode, r_keep_frac, r_hvg_weight):
    """Use anglemania's per-gene repeatability R to steer HVG selection.

    Three modes, all of which end up ranking genes inside anglemania's all-batch
    intersection (``R`` is only defined there):

    ``hvg_rfilter``
        Hard prefilter: drop the ``1 - r_keep_frac`` worst genes by binned ``R``
        z-score, then run the *same* HVG selection the ``hvg`` arm uses on what
        survives. "Don't let HVG pick a gene whose relational structure doesn't
        reproduce across batches."

    ``hvg_rweight``
        Soft version: rank on ``dispersions_norm + r_hvg_weight * R_z_binned``.
        Both terms are within-expression-bin z-scores (see ``hvg_dispersion``),
        so they are on the same scale and the weight is interpretable.

    ``hvg_intersect``
        The control: plain HVG, restricted to the intersection but with no ``R``
        involvement at all. Needed because both variants above also silently
        drop every gene missing from some batch, and on real data that is a
        large fraction -- without this arm, an ``R`` effect can't be told apart
        from an intersection effect.
    """
    from gene_level_scores import (
        hvg_dispersion,
        intersect_gene_names,
        score_genes,
    )

    if mode == "hvg_intersect":
        genes = intersect_gene_names(
            adata, batch_key=batch_key, dataset_key=dataset_key, min_cells_per_gene=1
        )
        print(f"Restricting to {len(genes)} intersected genes (of {adata.n_vars})",
              file=sys.stderr)
        return select_hvg(adata[:, genes].copy(), batch_key, n_genes)

    df = score_genes(
        adata,
        batch_key=batch_key,
        dataset_key=dataset_key,
        method=method,
        normalization_mode=normalization_mode,
        permutation_function=permutation_function,
        min_cells_per_gene=1,
        gpu=True,
        n_bins=n_bins,
    )

    if mode == "hvg_rfilter":
        n_keep = max(int(round(r_keep_frac * len(df))), n_genes)
        keep = df.nlargest(n_keep, "R_z_binned")["gene"].to_numpy()
        print(f"Keeping the top {len(keep)} of {len(df)} scored genes by R, "
              f"then selecting {n_genes} HVGs among them", file=sys.stderr)
        # select_hvg normalizes in place, so hand it a copy of the subset
        return select_hvg(adata[:, keep].copy(), batch_key, n_genes)

    # hvg_rweight
    disp = hvg_dispersion(adata, batch_key, n_top_genes=n_genes,
                          genes=df["gene"].to_numpy())
    merged = df.merge(disp, on="gene", how="inner", validate="1:1")
    assert len(merged) == len(df)
    merged["score"] = merged["dispersions_norm"] + r_hvg_weight * merged["R_z_binned"]
    return merged.nlargest(min(n_genes, len(merged)), "score")["gene"].tolist()


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
    elif args.gene_selection in ANGLGENE_VARIANTS:
        genes = select_anglgene(
            adata, args.gene_selection, args.batch_key, args.dataset_key, args.n_genes,
            args.anglemania_mode, args.permutation_function, args.normalization_mode,
            args.n_bins, args.weight_signal, args.weight_R,
            args.loess_span, args.moderate_n_bins, args.score_table,
        )
    elif args.gene_selection in ("hvg_rfilter", "hvg_rweight", "hvg_intersect"):
        genes = select_hvg_r(
            adata, args.batch_key, args.dataset_key, args.n_genes,
            args.anglemania_mode, args.permutation_function, args.normalization_mode,
            args.n_bins, args.gene_selection, args.r_keep_frac, args.r_hvg_weight,
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
