import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from balanced_clustering import balanced_v_measure
from scib.metrics import cluster_optimal_resolution
import click
import os
from os.path import join

#------------------------------------------------------------------------------#
# COMMAND LINE ARGUMENTS
#------------------------------------------------------------------------------#
@click.command()
@click.option("--original_h5ad",
              default=None,
              help="File location of the original h5ad input file")
@click.option("--embedding_tsv",
              default=None,
              help="File location of the embedding TSV file (cell barcodes as row index)")
@click.option("--outfile",
              default="/data/akalin/akollot/anglemania_benchmark/results/pbmcsca/hvg/",
              help="File location of the output file where the metrics")
@click.option("--sample",
              prompt="The sample name",
              help="Sample name")
@click.option("--integration_method",
              prompt="The implementation of the integration method to use from the scib package",
              help="Specify integration method. One of ['harmony','scvi','scanvi']",
              type=click.Choice(['harmony', 'scvi', 'scanvi', 'seurat', 'scanorama']))
@click.option("--gene_selection",
              prompt="The implementation of the integration method to use from the scib package",
              help="Specify integration method. One of ['harmony','scvi','scanvi']",
              type=click.Choice(['hvg', 'angl', 'anglmean', 'anglmod05', 'anglmod1', 'anglsiggene', 'anglgene', 'anglgene_nobin', 'anglgene_anchor', 'anglgene_biology', 'anglgene_noiselo', 'anglgene_noisehi', 'anglgene_mod', 'anglgene_anchor_mod', 'anglgene_biology_mod', 'hvg_rfilter', 'hvg_rweight', 'hvg_intersect', 'full', 'rand', 'topbtvr', 'bottombtvr']))
@click.option("--label_key",
              default="CellType",
              help="Label key specifying a column in the metadata that specifies the cell type information")

def calc_bNMI(
    original_h5ad,
    embedding_tsv,
    outfile,
    sample,
    gene_selection,
    integration_method,
    label_key
    ):
    """
    Calculate balanced NMI (bNMI) for an integrated AnnData object.
    """
    emb_df = pd.read_csv(embedding_tsv, sep="\t", index_col=0)
    adata_orig = sc.read_h5ad(original_h5ad)
    adata_integrated = adata_orig[emb_df.index, :].copy()
    adata_integrated.obsm["X_emb"] = emb_df.values
    adata_integrated.obs[label_key] = adata_integrated.obs[label_key].astype("category").cat.remove_unused_categories()

    # Build nearest neighbor graph on the integrated embedding
    sc.pp.neighbors(adata_integrated, use_rep = "X_emb")
    
    click.secho("[bNMI] Finding optimal clustering resolution...", fg='bright_yellow', err=True)

    # Find optimal Leiden resolution for balanced NMI
    cluster_optimal_resolution(
        adata_integrated,
        label_key=label_key,
        cluster_key="Cluster",
        metric=bNMI_metric
    )
    click.secho("[bNMI] Calculating bNMI...", fg='bright_yellow', err=True)
    # Compute final score
    bNMI = balanced_v_measure(
        labels_true=adata_integrated.obs[label_key],
        labels_pred=adata_integrated.obs["Cluster"]
    )

    metrics = pd.DataFrame(
        data = {
            "bNMI" : [bNMI],
            "sample" : [sample],
            "integration_method" : [integration_method],
            "gene_selection" : [gene_selection]
            }
    )
    #put method as the first 
    metrics.to_csv(outfile, sep = "\t")


def bNMI_metric(adata_integrated, label_key, cluster_key):
    """Metric wrapper for cluster_optimal_resolution."""
    clusters = adata_integrated.obs[cluster_key].to_numpy()
    labels = adata_integrated.obs[label_key].to_numpy()
    if len(np.unique(clusters)) == 1:
        return 0
    return balanced_v_measure(labels_true=labels, labels_pred=clusters)


if __name__ == "__main__":
    calc_bNMI()

