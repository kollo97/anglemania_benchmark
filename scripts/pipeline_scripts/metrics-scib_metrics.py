import click
# Python packages
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
# import scib
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import os
from os.path import join
from scib.metrics import cluster_optimal_resolution

#------------------------------------------------------------------------------#
# COMMAND LINE ARGUMENTS
#------------------------------------------------------------------------------#
@click.command()
@click.option("--original_h5ad",
              default="/data/akalin/akollot/anglemania_benchmark/data/pbmcsca/adata_hvg.h5ad",
              help="File location of the original h5ad input file")
@click.option("--integrated_h5ad",
              default="/data/akalin/akollot/anglemania_benchmark/results/pbmcsca/hvg/harmony.h5ad", 
              help="File location of the integrated h5ad file")
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
              type=click.Choice(['hvg', 'angl', 'full', 'rand']))
@click.option("--batch_key",
              default="batch",
              help="Batch key specifying a column in the metadata that divides the batches")
@click.option("--label_key",
              default="CellType",
              help="Label key specifying a column in the metadata that specifies the cell type information")

#------------------------------------------------------------------------------#
# MAIN FUNCTION
#------------------------------------------------------------------------------#

def calc_metrics(original_h5ad, 
                 integrated_h5ad, 
                 outfile,
                 sample,
                 gene_selection,
                 integration_method, 
                 batch_key, 
                 label_key):
    outfile=outfile
    if os.path.exists(outfile):
        # stop and throw message
        raise FileExistsError(f"{outfile} already exists. Please delete it and try again.")
        
    click.secho(f"Input files are: {original_h5ad}, {integrated_h5ad}", fg='bright_yellow', err=True)

    adata = sc.read_h5ad(original_h5ad)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["logcounts"] = adata.X.copy()
    
    adata_integrated = sc.read_h5ad(integrated_h5ad) # is already log normalized => located in .X 
    
    if integration_method == 'seurat':
        adata_integrated.layers["raw_counts"] = adata.X
        # make category type, so that the metrics can be calculated
        adata_integrated.obs[batch_key] = adata_integrated.obs[batch_key].astype("category")
        adata_integrated.obs[label_key] = adata_integrated.obs[label_key].astype("category")
        sc.tl.pca(adata_integrated, layer="data") # anndataR::write_h5ad() writes the default assay of the seurat object (which was "integrated") to the "data" layer. This is the log-normalized data.
        
        adata_integrated.obsm["X_emb"] = adata_integrated.obsm['X_pca']
    
    bm = Benchmarker(
        adata_integrated,
        batch_key = batch_key,
        label_key = label_key,
        bio_conservation_metrics = BioConservation(
            silhouette_label=True,
            isolated_labels=True,
            nmi_ari_cluster_labels_leiden=False,
            nmi_ari_cluster_labels_kmeans=False,
            clisi_knn=True
        ),
        batch_correction_metrics = BatchCorrection(
            bras=True,
            ilisi_knn=True,
            kbet_per_label=False, # also exclude this according to the feature selection benchmarking paper
            graph_connectivity=False,
            pcr_comparison=True
        ),
        embedding_obsm_keys=['X_emb'],
        n_jobs = 10
    )
    bm.benchmark()
    metrics = bm.get_results(min_max_scale=False).reset_index(drop=True)
    metrics = metrics.iloc[0:1,:]

    click.secho(f"Saving metrics to {outfile}", fg='bright_yellow', err=True)
    
    hvg_angl = os.path.basename(original_h5ad)
    if "anglemania" in hvg_angl:
        hvg_angl = "anglemania"
    elif "hvg" in hvg_angl:
        hvg_angl = "hvg"
    # metrics = metrics.T
    # IF SNAKEMAKE IS BEING USED:
    metrics["sample"] = sample
    metrics["integration_method"] = integration_method
    metrics["gene_selection"] = gene_selection
    #put method as the first 
    metrics.to_csv(outfile, sep = "\t")
    
if __name__ == "__main__":
    calc_metrics()