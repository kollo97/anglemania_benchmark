import click
# Python packages
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import harmony
import scvi
import os
from os.path import join

#------------------------------------------------------------------------------#
# COMMAND LINE ARGUMENTS
#------------------------------------------------------------------------------#
@click.command()
@click.option("--infile", 
              default="/data/akalin/akollot/anglemania_benchmark/data/pbmcsca/adata.h5ad", 
              help="File location of the original h5ad input file")
@click.option("--feature_subset", 
              default="/data/akalin/akollot/anglemania_benchmark/data/pbmcsca/adata_hvg.tsv", 
              help="File location of the h5ad input file")
@click.option("--outfile", 
              default="/data/akalin/akollot/anglemania_benchmark/results/pbmcsca/integrated/test.h5ad", 
              help="File location of the output file where the integrated data/embeddings should be stored")
@click.option("--integration_method", 
              prompt="The implementation of the integration method to use from the scib package", 
              help="Specify integration method. One of ['harmony','scvi','scanvi']", 
              type=click.Choice(["harmony","scvi","scanvi"]))
@click.option("--batch_key", 
              default="batch", 
              help="Batch key specifying a column in the metadata that divides the batches")
@click.option("--label_key", 
              default="CellType", 
              help="Label key specifying a column in the metadata that specifies the cell type information")

#------------------------------------------------------------------------------#
# MAIN FUNCTION
#------------------------------------------------------------------------------#

def integrate(infile, feature_subset, outfile, integration_method, batch_key, label_key):
    click.secho(f"Input file is: {infile}", fg="bright_yellow", err=True)
    click.secho("Reading input file...", fg="bright_yellow", err=True)

    selected_features  = pd.read_csv(feature_subset, sep = "\t")
    selected_features = selected_features.hgnc_symbol.to_list()
    adata = ad.read_h5ad(infile)
    
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["logcounts"] = adata.X.copy()

    adata = adata[:, selected_features].copy()
    #and layers["data"] is the normalized data
    
    if integration_method == "harmony":
        click.secho("Performing PCA...", fg="bright_yellow", err=True)
        if "X_pca" not in adata.obsm:
            sc.tl.pca(adata) # pca on adata.X i.e. normalized counts
                    # integrated = join(out_integration, "{sample}/{sample}_{integration_method}.h5ad")
        click.secho("Integrating using Harmony...", fg="bright_yellow", err=True)
        adata.obsm["X_emb"] = harmony.harmonize(adata.obsm["X_pca"], adata.obs, batch_key=batch_key)
        click.secho(f"Saving Harmony integrated data to {outfile}", fg="bright_yellow", err=True)
        adata.write(outfile)

    elif integration_method == "scvi":
        click.secho("Setting up SCVI model...", fg="bright_yellow", err=True)
        scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
        max_epochs_scvi = np.min([round((20000 / adata.n_obs) * 400), 400])
        
        click.secho("Training SCVI model...", fg="bright_yellow", err=True)
        model = scvi.model.SCVI(adata)
        model.train()
        adata.obsm["X_emb"] = model.get_latent_representation()
        
        click.secho(f"Saving SCVI integrated data to {outfile}", fg="bright_yellow", err=True)
        adata.write_h5ad(outfile)
    
    elif integration_method == "scanvi":
        click.secho("Setting up SCVI model...", fg="bright_yellow", err=True)
        scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
        max_epochs_scvi = np.min([round((20000 / adata.n_obs) * 400), 400])
        
        click.secho("Training SCVI model...", fg="bright_yellow", err=True)
        model = scvi.model.SCVI(adata)
        model.train()
        adata.obsm["X_emb"] = model.get_latent_representation()
        max_epochs_scanvi = int(np.min([10, np.max([2, round(max_epochs_scvi / 3.0)])]))
        
        click.secho("Setting up SCANVI model...", fg="bright_yellow", err=True)
        model = scvi.model.SCANVI.from_scvi_model(model, labels_key=label_key, unlabeled_category="unlabelled")
        
        click.secho("Training SCANVI model...", fg="bright_yellow", err=True)
        model.train(max_epochs=max_epochs_scanvi)
        adata.obsm["X_emb"] = model.get_latent_representation()
        
        click.secho(f"Saving SCANVI integrated data to {outfile}", fg="bright_yellow", err=True)
        adata.write_h5ad(outfile)

    click.secho("Integration process completed successfully.", fg="bright_yellow", err=True)




if __name__ == "__main__":
    integrate()
