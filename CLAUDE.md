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
- `anglemania` R package (the legacy R implementation, superseded by `pyanglemania` for gene selection)
- `CellMixS`, `balanced_clustering` Python package

Gene selection (`scripts/pipeline_scripts/prepare_inputs.py`) runs in a separate **`pyanglemania` conda environment**, not in the Guix environment, because it depends on the GPU-accelerated [`pyanglemania`](https://github.com/omnideconv/pyanglemania) package (`~/projects/pyanglemania`, env spec at `~/projects/pyanglemania/envs/pyanglemania.yml`) and on `scikit-misc` (only needed if `hvg` is switched to scanpy's `seurat_v3` flavor; not used by the current `seurat` flavor, kept installed for that option). Sync it with `mamba env update -f ~/projects/pyanglemania/envs/pyanglemania.yml`. The Snakemake rules invoke it via `conda run -n pyanglemania`, so `snakemake` itself still runs from `bm_guix`.

## Running the Pipeline

The main pipeline is a Snakemake workflow. All commands run from `scripts/`:

```bash
cd scripts/

# Dry run to check pipeline
snakemake -s Snakefile.smk --configfile config_files/<config>.yml -n

# Run pipeline (local)
snakemake -s Snakefile.smk --configfile config_files/<config>.yml --cores <N>

# Run with SLURM (preferred — each rule becomes a separate SLURM job)
snakemake -s Snakefile.smk --configfile config_files/<config>.yml --profile slurm_profile/
```

### SLURM execution details

The SLURM profile is at `scripts/slurm_profile/config.yaml`. It translates each
rule's `resources` block into `sbatch` flags. Resource defaults:

| Resource      | Default | Notes                                          |
|---------------|---------|------------------------------------------------|
| `cpus_per_task` | 4     | `--cpus-per-task`                              |
| `mem_mb`      | 16 000  | `--mem` in MB                                  |
| `runtime`     | 60      | `--time` in minutes                            |
| `gres`        | —       | `--gres`; GPU rules use `gres="gpu:1"`         |
| `gpu_slots`   | —       | not an sbatch flag; local Snakemake resource pool (see below) |

`preprocess_gpu` is the only rule that sets `gres="gpu:1"` — `integrate` runs
scvi/scanvi on CPU (Guix's python-pytorch has no CUDA). If your cluster routes
GPU jobs through a dedicated partition, also add `partition="gpu"` (or
equivalent) to that rule.

**Important:** Do NOT add `gpu=N` to a rule's `resources` block. The native
SLURM executor plugin converts any resource named `gpu` into `--gpu=N` as an
sbatch flag, which this SLURM version rejects as ambiguous (exit 255) and
causes the job to fail silently. Do NOT use `slurm_extra="--gres=gpu:1"` —
when multiple GPU jobs submit simultaneously, the `slurm_extra` code path in
the executor deadlocks (all submission threads wait on the same mutex). Use
`gres="gpu:1"` instead; SLURM's scheduler handles concurrent GPU allocation.

`preprocess_gpu` also sets `gpu_slots=1`, and `slurm_profile/config.yaml` caps the
global `gpu_slots` pool at 1 (`resources: [gpu_slots=1]`). Since `gres` isn't
a resource Snakemake itself schedules on (it's only translated into an sbatch
flag), without this Snakemake will submit as many `gres="gpu:1"` jobs at once
as `jobs:` allows, and SLURM's own scheduler ends up handing out however many
GPUs are actually free on the cluster concurrently. `gpu_slots` is a plain,
arbitrarily-named Snakemake-local resource (not recognized by the SLURM
executor, so it's never turned into an sbatch flag) that throttles how many
`gres="gpu:1"` jobs Snakemake will have in flight pipeline-wide at once. Raise
the `gpu_slots` cap in `slurm_profile/config.yaml` to allow more concurrent
GPU jobs.

SLURM logs go to `scripts/logs/slurm/` (created automatically on pipeline start).

**No setup needed.** `snakemake-executor-plugin-slurm` and
`-slurm-jobstep` are defined in `guix/manifest.scm`, so they are on the path
inside `bm_guix`. Nothing has to be pip-installed.

**Important:** run Snakemake from inside `bm_guix` so that the Guix python/R
binaries and the executor plugins are on PATH — the profile passes
`--export=ALL` to sbatch, so each submitted job inherits the full Guix
environment.

Config files for different experiments are in `scripts/config_files/`. Key config parameters:
- `samplesheet`: TSV/CSV/JSON file with `sample_name` and `file_path` columns (pointing to `.h5ad` files)
- `output_dir`: root output directory
- `integration_methods`: currently only `["scvi", "scanorama"]` are used
- `gene_selection`: currently only `["hvg", "angl"]` are compared (anglemania vs HVG)
- `batch_key` / `label_key`: column names in `.obs` metadata
- `anglemania_mode`: `"cosine"` | `"spearman"` | `"phi_s"`
- `permutation_function`: `"sample"` | `"permute_nonzero"`
- `normalization_mode`: `"classical"` (CP10K + log1p) | `"pflog1ppf"` (shifted-CLR). Only affects `angl` (passed through to pyanglemania's `normalization_method`). `hvg` always normalizes with plain CP10K + log1p before scanpy's dispersion-based `seurat` flavor — `pflog1ppf` is incompatible with it (its per-cell mean-centering breaks `seurat`'s dispersion calculation; verified ~40% of genes get NaN dispersion), so it's hardcoded rather than exposed as a choice there.
- `n_genes`: integer (default 2000)

## Repository Layout

`scripts/` separates pipeline code from one-off analysis code:

```
scripts/
  Snakefile.smk           # workflow entrypoint
  rules/                  # preprocess.smk, integrate.smk, metrics.smk
  pipeline_scripts/       # everything the Snakemake rules invoke (.py and .R)
  config_files/           # *.yml configs + samplesheets
  slurm_profile/          # SLURM executor profile
  R/                      # standalone R analysis/plotting scripts (not called by the pipeline)
  python/                 # standalone Python analysis scripts (not called by the pipeline)
  utils/                  # R helpers sourced by scripts/R/* and the notebooks
  notebooks/              # .Rmd visualization and simulation notebooks
  *.sh                    # shell wrappers (run_genelevel_variant_lists.sh)
```

Only `pipeline_scripts/` is referenced from `rules/*.smk`. The `R/` and `python/`
scripts are run by hand — their header comments give the exact invocation, and
paths in them are relative to the repo root, so run them from there.

## Pipeline Architecture

The Snakemake workflow (`scripts/Snakefile.smk`) chains three rule files:

```
rules/preprocess.smk   → scripts/pipeline_scripts/prepare_inputs.py
                          rule preprocess_cpu (hvg/full/rand/topbtvr/bottombtvr, CPU)
                          rule preprocess_gpu (angl, resources: gres="gpu:1", gpu_slots=1)
rules/integrate.smk    → scripts/pipeline_scripts/integration.py  (Python methods)
                       → scripts/pipeline_scripts/seurat_integration.R (Seurat method)
rules/metrics.smk      → scripts/pipeline_scripts/metrics-scib_metrics.py  (scib_metrics)
                       → scripts/pipeline_scripts/metrics-bNMI.py           (balanced NMI)
                       → scripts/pipeline_scripts/metrics-cms.R             (CellMixS CMS)
                       → scripts/pipeline_scripts/metrics-ldfDiff.R         (ldfDiff)
                       → scripts/pipeline_scripts/combine_metrics.py        (merge all)
```

`preprocess_cpu` and `preprocess_gpu` both produce `preprocessed/{sample}_{gene_selection}.tsv`; each rule restricts which `gene_selection` values it can produce via `wildcard_constraints`, so Snakemake resolves the right rule per file with no input function needed — `rules/integrate.smk` just depends on the shared `PREPROCESS_OUTPUT` path pattern (defined in `Snakefile.smk`).

**Data flow:**
1. **Preprocess** (`prepare_inputs.py`): reads `.h5ad`, selects genes using the specified method (`hvg` → CP10K + log1p then `scanpy.pp.highly_variable_genes(flavor="seurat", batch_key=...)`, `angl` → `pyanglemania.preprocessing.anglemania` on GPU with `normalization_method` set from `normalization_mode`, plus `full`/`rand`/`topbtvr`/`bottombtvr`), writes a TSV with column `hgnc_symbol`.
2. **Integrate** (`integration.py` or `seurat_integration.R`): reads original `.h5ad` + gene list TSV, runs integration, writes the cell embedding as a TSV to `embedding/{sample}/{method}/{sample}_{gene_selection}.tsv`.
3. **Metrics** (4 parallel scripts): each takes `--original_h5ad` + `--embedding_tsv`, computes one metric group, writes a TSV with columns `[metric, sample, integration_method, gene_selection]`.
4. **Combine** (`combine_metrics.py`): merges all metric TSVs, recomputes Bio conservation / Batch correction / Total scores (60/40 weighted mean).

**Metrics computed:**
- Bio conservation: Silhouette label, Isolated labels, bNMI, cLISI, ldfDiff
- Batch correction: PCR comparison, iLISI, CMS, BRAS
- Total = 0.6 × Bio conservation + 0.4 × Batch correction

## Running Individual Scripts Manually

```bash
# Gene selection (runs in the `pyanglemania` conda env: conda run -n pyanglemania ...)
python3 scripts/pipeline_scripts/prepare_inputs.py \
    --infile data/sample.h5ad --outfile out/genes.tsv \
    --batch_key Batch --gene_selection angl

# Integration (Python methods)
python3 scripts/pipeline_scripts/integration.py \
    --infile=data/sample.h5ad --feature_subset=out/genes.tsv \
    --outfile=out/embedding.tsv --integration_method=harmony \
    --batch_key=Batch --label_key=Group

# scib metrics
python3 scripts/pipeline_scripts/metrics-scib_metrics.py \
    --original_h5ad=data/sample.h5ad --embedding_tsv=out/embedding.tsv \
    --outfile=out/metrics.tsv --sample=mysample \
    --gene_selection=angl --integration_method=harmony \
    --batch_key=Batch --label_key=Group
```

## Output Structure

```
output_dir/
  preprocessed/           # {sample}_{gene_selection}.tsv  (gene lists)
  embedding/
    {sample}/{method}/    # {sample}_{gene_selection}.tsv (cell embeddings)
  metrics/
    {sample}/{method}/    # *_scibmetrics.tsv, *_bnmi.tsv, *_cms.tsv, *_ldfDiff.tsv
    {name}_combined_metrics.tsv
```

## Key Notes

- `integration.py` stores the embedding in `obsm["X_emb"]` internally, then writes it out as the embedding TSV the metrics rules consume
- scanorama requires cells sorted by batch (done automatically in `integration.py`)
- The `JAX_PLATFORMS=cpu` environment variable is set before Python integrations to prevent JAX GPU conflicts
- GPU use is confined to gene selection and regulated in three places: both preprocess rules shell out to `conda run -n pyanglemania ...` (that env, not Guix, has CUDA-enabled `cupy` + pyanglemania); `preprocess_gpu` sets `gres="gpu:1"` to get a device from SLURM; and `gpu_slots=1` caps pipeline-wide GPU concurrency (see SLURM execution details above). `prepare_inputs.py` raises if no CUDA device is reachable — `angl`/`anglgene_*` have no CPU fallback
- The `gpu:` key in some config files is vestigial — nothing in `Snakefile.smk`, `rules/` or `pipeline_scripts/` reads it
- Notebooks for visualization and simulation are in `scripts/notebooks/`
- Utility ggplot functions shared across notebooks are in `scripts/R/ggplot_utils.R` (the full set: `draw_heatmap*`, `theme_big_text`, `today`) and `scripts/utils/ggplot_utils.R` (only `big_text_theme`) — two separate files, not copies
