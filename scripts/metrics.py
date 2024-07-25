import click
# Python packages
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import scib
import os
from os.path import join

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
@click.option("--integration_method", 
              prompt="The implementation of the integration method to use from the scib package", 
              help="Specify integration method. One of ['harmony','scvi','scanvi']", 
              type=click.Choice(['harmony','scvi','scanvi', 'seurat']))
@click.option("--batch_key", 
              default="batch", 
              help="Batch key specifying a column in the metadata that divides the batches")
@click.option("--label_key", 
              default="CellType", 
              help="Label key specifying a column in the metadata that specifies the cell type information")

#------------------------------------------------------------------------------#
# MAIN FUNCTION
#------------------------------------------------------------------------------#

def calc_metrics(original_h5ad, integrated_h5ad, outfile, integration_method, batch_key, label_key):
    outfile=outfile
    if os.path.exists(outfile):
        # stop and throw message
        raise FileExistsError(f"{outfile} already exists. Please delete it and try again.")
        
    click.secho(f"Input files are: {original_h5ad}, {integrated_h5ad}", fg='bright_yellow', err=True)

    adata = sc.read_h5ad(original_h5ad)
    adata_integrated = sc.read_h5ad(integrated_h5ad)
    if integration_method == 'seurat':
        # apparently the SeuratDisk Convert function produces an old format of the h5ad format, so just reconstruct it
        counts = adata_integrated.X.copy()
        meta = adata_integrated.obs.copy()
        cell_names = adata_integrated.obs_names
        gene_names = adata_integrated.var_names
        adata_integrated = ad.AnnData(counts)
        adata_integrated.obs_names = cell_names
        adata_integrated.var_names = gene_names
        adata_integrated.obs = meta
        # make category type, so that the metrics can be calculated
        adata_integrated.obs[batch_key] = adata_integrated.obs[batch_key].astype("category")
        adata_integrated.obs[label_key] = adata_integrated.obs[label_key].astype("category")

        adata_integrated
        sc.tl.pca(adata_integrated)
        adata_integrated.obsm["X_emb"] = adata_integrated.obsm['X_pca']
    
    selected_metrics = {
        'ari_' : True,
        'nmi_' : True,
        'silhouette_' : True,
        'pcr_' : True,
        'hvg_score_' : True,
        'isolated_labels_' : True,
        'isolated_labels_f1_' : True,
        'isolated_labels_asw_' : True,
        'graph_conn_' : True,
        'kBET_' : True,
        'lisi_graph_' : True,
        'ilisi_' : True,
        'clisi_' : True,
        'n_cores' : 20
    }
    click.secho(f"Calculating metrics: {selected_metrics}", fg='bright_yellow', err=True)

    metrics = scib.metrics.metrics(adata, adata_integrated, batch_key = batch_key, label_key = label_key, embed = 'X_emb', **selected_metrics)

    click.secho(f"Saving metrics to {outfile}", fg='bright_yellow', err=True)
    
    hvg_angl = os.path.basename(original_h5ad)
    if "anglemania" in hvg_angl:
        hvg_angl = "anglemania"
    elif "hvg" in hvg_angl:
        hvg_angl = "hvg"
    method = {"method" : f"{integration_method}_{hvg_angl}"}
    method = pd.DataFrame.from_dict(method, orient="index")
    metrics = pd.concat([method, metrics], axis = "rows")
    
    #put method as the first 
    metrics.to_csv(outfile, sep = "\t")
    
if __name__ == "__main__":
    calc_metrics()