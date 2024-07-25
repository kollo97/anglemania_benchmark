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
              default="/data/akalin/akollot/anglemania_benchmark/data/pbmcsca/adata_hvg.h5ad", 
              help="File location of the h5ad input file")
@click.option("--outdir", 
              default="/data/akalin/akollot/anglemania_benchmark/results/hvg/", 
              help="File location of the output dir where the integrated data/embeddings should be stored")
@click.option("--integration_method", 
              prompt="The implementation of the integration method to use from the scib package", 
              help="Specify integration method. One of ['harmony','scvi','scanvi']", 
              type=click.Choice(['harmony','scvi','scanvi']))
@click.option("--batch_key", 
              default="batch", 
              help="Batch key specifying a column in the metadata that divides the batches")
@click.option("--label_key", 
              default="CellType", 
              help="Label key specifying a column in the metadata that specifies the cell type information")

#------------------------------------------------------------------------------#
# MAIN FUNCTION
#------------------------------------------------------------------------------#

def integrate(infile, outdir, integration_method, batch_key, label_key):
    click.secho(f"Input file is: {infile}", fg='bright_yellow', err=True)
    click.secho("Reading input file...", fg='bright_yellow', err=True)
    adata = ad.read_h5ad(infile)
    
    if integration_method == 'harmony':
        click.secho("Performing PCA...", fg='bright_yellow', err=True)
        if 'X_pca' not in adata.obsm:
            sc.tl.pca(adata)
            
        click.secho("Integrating using Harmony...", fg='bright_yellow', err=True)
        adata.obsm["X_emb"] = harmony.harmonize(adata.obsm['X_pca'], adata.obs, batch_key=batch_key)
        output_file = join(outdir, 'harmony.h5ad')
        click.secho(f"Saving Harmony integrated data to {output_file}", fg='bright_yellow', err=True)
        adata.write(output_file)

    elif integration_method in ['scvi', 'scanvi']:
        click.secho("Setting up SCVI model...", fg='bright_yellow', err=True)
        scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
        max_epochs_scvi = np.min([round((20000 / adata.n_obs) * 400), 400])
        
        click.secho("Training SCVI model...", fg='bright_yellow', err=True)
        model = scvi.model.SCVI(adata)
        model.train()
        adata.obsm["X_emb"] = model.get_latent_representation()
        
        scvi_outfile = join(outdir, 'scvi.h5ad')
        click.secho(f"Saving SCVI integrated data to {scvi_outfile}", fg='bright_yellow', err=True)
        adata.write_h5ad(scvi_outfile)
    
        max_epochs_scanvi = int(np.min([10, np.max([2, round(max_epochs_scvi / 3.0)])]))
        click.secho("Setting up SCANVI model...", fg='bright_yellow', err=True)
        model = scvi.model.SCANVI.from_scvi_model(model, labels_key=label_key, unlabeled_category="unlabelled")
        
        click.secho("Training SCANVI model...", fg='bright_yellow', err=True)
        model.train(max_epochs=max_epochs_scanvi)
        adata.obsm["X_emb"] = model.get_latent_representation()
        
        scanvi_outfile = join(outdir, 'scanvi.h5ad')
        click.secho(f"Saving SCANVI integrated data to {scanvi_outfile}", fg='bright_yellow', err=True)
        adata.write_h5ad(scanvi_outfile)

    click.secho("Integration process completed successfully.", fg='bright_yellow', err=True)




if __name__ == "__main__":
    integrate()
