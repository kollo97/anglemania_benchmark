# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a benchmarking pipeline for [anglemania](https://github.com/omnideconv/anglemania), an R package for gene selection in single-cell RNA-seq batch integration. The benchmark compares anglemania gene selection (`angl`) against highly variable gene selection (`hvg`), full gene sets (`full`), and random gene selection (`rand`) across multiple integration methods and datasets.

## Environment Setup

All commands for this project are run inside the **Guix environment** using the `bm_guix` shell alias (defined in `~/.bashrc`):

```bash
bm_guix
# expands to:
# guix time-machine -C guix/channels.scm -- shell -m guix/manifest.scm --no-grafts python-numpyro
```

The Guix environment spec is in `guix/channels.scm` and `guix/manifest.scm`.

A conda environment (`angl_BM2`) is also defined in `envs/env.yml` as an alternative, with additional R packages that must be installed manually:
- `kBET`: `devtools::install_github("theislab/kBET")`
- `schard`: `devtools::install_github("cellgeni/schard")`
- `anglemania` R package (the package being benchmarked)
- `CellMixS`, `balanced_clustering` Python package

## Running the Pipeline

The main pipeline is a Snakemake workflow. All commands run from `scripts/`:

```bash
cd scripts/

# Dry run to check pipeline
snakemake -s Snakefile.smk --configfile config_files/<config>.yml -n

# Run pipeline (local)
snakemake -s Snakefile.smk --configfile config_files/<config>.yml --cores <N>

# Run with SLURM
snakemake -s Snakefile.smk --configfile config_files/<config>.yml --cluster "sbatch ..."
```

Config files for different experiments are in `scripts/config_files/`. Key config parameters:
- `samplesheet`: TSV/CSV/JSON file with `sample_name` and `file_path` columns (pointing to `.h5ad` files)
- `output_dir`: root output directory
- `integration_methods`: currently only `["scvi", "scanorama"]` are used
- `gene_selection`: currently only `["hvg", "angl"]` are compared (anglemania vs HVG)
- `batch_key` / `label_key`: column names in `.obs` metadata
- `anglemania_mode`: `"cosine"` | `"spearman"` | `"diem"`
- `permutation_function`: `"sample"` | `"permute_nonzero"`
- `n_genes`: integer (default 2000)

## Pipeline Architecture

The Snakemake workflow (`scripts/Snakefile.smk`) chains three rule files:

```
rules/preprocess.smk   → scripts/pipeline_scripts/prepare_inputs.R
rules/integrate.smk    → scripts/pipeline_scripts/integration.py  (Python methods)
                       → scripts/seurat_integration.R             (Seurat method)
rules/metrics.smk      → scripts/pipeline_scripts/metrics-scib_metrics.py  (scib_metrics)
                       → scripts/pipeline_scripts/metrics-bNMI.py           (balanced NMI)
                       → scripts/pipeline_scripts/metrics-cms.R             (CellMixS CMS)
                       → scripts/pipeline_scripts/metrics-ldfDiff.R         (ldfDiff)
                       → scripts/pipeline_scripts/combine_metrics.py        (merge all)
```

**Data flow:**
1. **Preprocess** (`prepare_inputs.R`): reads `.h5ad`, selects genes using the specified method, writes a TSV with column `hgnc_symbol`.
2. **Integrate** (`integration.py` or `seurat_integration.R`): reads original `.h5ad` + gene list TSV, runs integration, writes integrated `.h5ad` with embedding in `obsm["X_emb"]`.
3. **Metrics** (4 parallel scripts): each reads original + integrated `.h5ad`, computes one metric group, writes a TSV with columns `[metric, sample, integration_method, gene_selection]`.
4. **Combine** (`combine_metrics.py`): merges all metric TSVs, recomputes Bio conservation / Batch correction / Total scores (60/40 weighted mean).

**Metrics computed:**
- Bio conservation: Silhouette label, Isolated labels, bNMI, cLISI, ldfDiff
- Batch correction: PCR comparison, iLISI, CMS, BRAS
- Total = 0.6 × Bio conservation + 0.4 × Batch correction

## Running Individual Scripts Manually

```bash
# Gene selection
Rscript scripts/pipeline_scripts/prepare_inputs.R \
    --infile data/sample.h5ad --outfile out/genes.tsv \
    --batch_key Batch --gene_selection angl

# Integration (Python methods)
python3 scripts/pipeline_scripts/integration.py \
    --infile=data/sample.h5ad --feature_subset=out/genes.tsv \
    --outfile=out/integrated.h5ad --integration_method=harmony \
    --batch_key=Batch --label_key=Group

# scib metrics
python3 scripts/pipeline_scripts/metrics-scib_metrics.py \
    --original_h5ad=data/sample.h5ad --integrated_h5ad=out/integrated.h5ad \
    --outfile=out/metrics.tsv --sample=mysample \
    --gene_selection=angl --integration_method=harmony \
    --batch_key=Batch --label_key=Group
```

## Output Structure

```
output_dir/
  preprocessed/           # {sample}_{gene_selection}.tsv  (gene lists)
  integration/
    {sample}/{method}/    # {sample}_{gene_selection}.h5ad (integrated data)
  metrics/
    {sample}/{method}/    # *_scibmetrics.tsv, *_bnmi.tsv, *_cms.tsv, *_ldfDiff.tsv
    {name}_combined_metrics.tsv
```

## Key Notes

- All integrated objects must store their embedding in `obsm["X_emb"]`
- scanorama requires cells sorted by batch (done automatically in `integration.py`)
- The `JAX_PLATFORMS=cpu` environment variable is set before Python integrations to prevent JAX GPU conflicts
- GPU resources are assigned dynamically for `scvi`/`scanvi` via Snakemake `resources: gpu`
- Notebooks for visualization and simulation are in `scripts/notebooks/`
- Utility ggplot functions shared across notebooks are in `scripts/utils/ggplot_utils.R` and `scripts/ggplot_utils.R`
