# anglemania Benchmark

Benchmarking pipeline for [anglemania](https://github.com/omnideconv/anglemania), a gene selection method for single-cell RNA-seq batch integration. The pipeline systematically compares **anglemania** (`angl`) against **highly variable gene** (`hvg`) selection across multiple integration methods and datasets, using a comprehensive set of batch correction and biological conservation metrics.

Gene selection itself runs through [pyanglemania](https://github.com/omnideconv/pyanglemania), the GPU-accelerated Python port; the original R package is still used by the standalone analysis scripts.

## Experimental Design

| Dimension | Options |
|-----------|---------|
| Gene selection | `angl` (anglemania), `hvg` (scanpy `seurat` HVGs) — plus the per-gene and HVG-hybrid arms below |
| Integration methods | `scvi`, `scanorama` (`harmony`, `scanvi`, `seurat` also implemented) |
| Datasets | Splatter-simulated (varying batch/bio effects, mixed effects, dropout), NeurIPS 2021 CITE-seq, Luecken & Theis NatMethods 2021 |
| Bio conservation | Silhouette label, Isolated labels, bNMI, cLISI, ldfDiff |
| Batch correction | PCR comparison, iLISI, CMS, BRAS |
| Final score | **Total = 0.6 × Bio conservation + 0.4 × Batch correction** |

## Environment Setup

Snakemake and everything except gene selection run inside the **Guix environment**:

```bash
# Via shell alias (defined in ~/.bashrc):
bm_guix
# Equivalent to:
guix time-machine -C guix/channels.scm -- shell -m guix/manifest.scm --no-grafts python-numpyro
```

Gene selection (`prepare_inputs.py`) instead runs in a dedicated `pyanglemania` conda environment (GPU-accelerated, spec at `~/projects/pyanglemania/envs/pyanglemania.yml`), which the Snakemake rules invoke via `conda run -n pyanglemania`. Snakemake itself is still started from `bm_guix`.

Keep the conda env in sync with:

```bash
mamba env update -f ~/projects/pyanglemania/envs/pyanglemania.yml
```

Nothing else needs installing for cluster runs: the SLURM executor plugins
(`snakemake-executor-plugin-slurm` and `-slurm-jobstep`) are part of
`guix/manifest.scm`, so they are already on the path inside `bm_guix`. Start
Snakemake from that shell and the `--profile slurm_profile/` run works as-is.
The profile passes `--export=ALL` to sbatch, so each submitted job inherits the
same Guix environment.

### Which stages use the GPU

Only gene selection does, and only through the conda env — it is regulated in
three places:

| Layer | Mechanism |
|-------|-----------|
| Environment | Both preprocess rules shell out to `conda run -n pyanglemania python3 pipeline_scripts/prepare_inputs.py`. That env, not Guix, provides CUDA-enabled `cupy` and pyanglemania. |
| Allocation | `preprocess_gpu` sets `resources: gres="gpu:1"`, which the SLURM profile turns into `sbatch --gres=gpu:1`. It is the **only** rule that requests a GPU. |
| Concurrency | The same rule sets `gpu_slots=1` against the profile's global `gpu_slots` pool, capping how many GPU jobs Snakemake has in flight pipeline-wide. |

`prepare_inputs.py` moves the matrix onto the device with `cupy` before calling
pyanglemania and **raises** if no CUDA device is reachable — the `angl` and
`anglgene_*` arms have no CPU fallback. `preprocess_cpu` uses the same conda
env but never touches the GPU, so its arms (`hvg`, `full`, `rand`, `topbtvr`,
`bottombtvr`, `hvg_intersect`) run anywhere.

Integration is the opposite case: it runs under Guix, whose `python-pytorch` is
built without CUDA, so scVI/scanVI always run on CPU. The `integrate` rule
accordingly requests no `gres` — don't set `accelerator="gpu"` in
`integration.py` expecting it to work.

> The `gpu:` key present in some config files is vestigial — nothing in
> `Snakefile.smk`, `rules/` or `pipeline_scripts/` reads it. GPU use is decided
> by the rule a `gene_selection` value routes to, not by config.

## Usage

All pipeline commands run from the `scripts/` directory, inside `bm_guix`:

```bash
cd scripts/

# Dry run — check what will be executed
snakemake -s Snakefile.smk --configfile config_files/<config>.yml -n

# Cluster run (preferred — one SLURM job per rule)
snakemake -s Snakefile.smk --configfile config_files/<config>.yml --profile slurm_profile/

# Local run
snakemake -s Snakefile.smk --configfile config_files/<config>.yml --cores 8
```

### SLURM profile

`scripts/slurm_profile/config.yaml` translates each rule's `resources` block into
`sbatch` flags. Defaults: `cpus_per_task=4`, `mem_mb=16000`, `runtime=60`.
`preprocess_gpu` additionally requests `gres="gpu:1"`, and pipeline-wide GPU
concurrency is capped separately by a Snakemake-local `gpu_slots` pool (default
1 — raise it in the profile to allow more concurrent GPU jobs); see
[Which stages use the GPU](#which-stages-use-the-gpu). Logs go to
`scripts/logs/slurm/`, created automatically on pipeline start.

Two gotchas baked into the profile, both of which cause silent job failures if
undone:

- Do **not** name a resource `gpu` — the SLURM executor turns it into `--gpu=N`, which this SLURM version rejects as ambiguous (exit 255).
- Do **not** use `slurm_extra="--gres=gpu:1"` — concurrent GPU submissions deadlock in that code path. Use `gres="gpu:1"`.

### Available Configs

See `docs/DATASETS.md` for the datasets each config runs on, with sizes.

| Config file | Dataset |
|-------------|---------|
| `LueckenTheis_config.yml` | Luecken, Theis NatMethods 2021 (all 6 samples) |
| `NeurIPS_CITE-seq.yml` | NeurIPS 2021 CITE-seq (full) |
| `NeurIPS_CITE-seq_Tcells.yml` | NeurIPS 2021 CITE-seq (T cells) |
| `NeurIPS_CITE-seq_filtered_no_ribo.yml` | NeurIPS 2021 CITE-seq without ribosomal genes |
| `config_varying_batch_biovar_facScale0.1.yml` | Splatter `de.facLoc 0.1`, mild `de.facScale` (54 samples) |
| `config_varying_batch_biovar_facScale0.2.yml` | Splatter `de.facLoc 0.1`, stronger `de.facScale` (54 samples) |
| `config_varying_batch_biovar_dropout_sample.yml` | Splatter-simulated with dropout (54 samples) |
| `Theis_dropout_simulations.yml` | Simulations based on Luecken, Theis Splatter parameters (20 samples) |
| `config_mixed_effects.yml` | Mixed-effects stability set (baseline + 11 perturbations) |
| `config_mixed_effects_meanonly.yml` | `anglmean` ranking counterfactual (baseline + 4 LOBO) |
| `config_mixed_effects_altranking.yml` | Alternative pair-ranking criteria (baseline + 4 LOBO) |
| `config_genelevel_pancreas.yml` | Per-gene signal/R score arms on the pancreas dataset |
| `config_genelevel_pancreas_cosine_diag.yml` | As above, `cosine` + `classical` diagnostic contrast |
| `config_genelevel_sim.yml` | Per-gene score arms on the mixed-effects baseline |
| `config_genelevel_batchsweep.yml` | Per-gene score across `batch.facLoc` {0.1, 0.3, 0.6} |
| `config_pyanglemania_test.yml` | Single-sample smoke test for the pyanglemania rules |

### Key Config Parameters

```yaml
name: my_experiment
samplesheet: path/to/samplesheet.tsv   # columns: sample_name, file_path (.h5ad)
output_dir: path/to/output/
out_preprocessed: preprocessed/         # sub-dirs of output_dir, one per stage
out_integration: integration/
out_embedding: embedding/
out_metrics: metrics/

batch_key: Batch                        # obs column for batch labels
label_key: Group                        # obs column for cell type labels
integration_methods: [scvi, scanorama]
gene_selection: [hvg, angl]
anglemania_mode: cosine                 # cosine | spearman | phi_s
permutation_function: sample            # sample | permute_nonzero
normalization_mode: classical           # classical (CP10K + log1p) | pflog1ppf (shifted-CLR, angl only)
n_genes: 2000
```

`normalization_mode` only affects `angl`: it is passed through to pyanglemania's
`normalization_method`. `hvg` always normalizes with plain CP10K + log1p, because
`pflog1ppf`'s per-cell mean-centering breaks scanpy's dispersion calculation
(~40% of genes end up with NaN dispersion).

Optional parameters, only read by the per-gene and HVG-hybrid arms:

```yaml
anglgene_n_bins: 20                     # equal-frequency mean-expression bins
anglgene_weight_signal: 0.5             # weight of the binned signal z-score
anglgene_weight_R: 0.5                  # weight of the binned R z-score
anglgene_loess_span: 0.5                # LOESS span for the two residual axes
anglgene_moderate_n_bins: 10            # bins for the moderated-S prior trend
anglgene_score_table_dir: path/to/dir   # reuse a precomputed M/S score table
r_keep_frac: 0.5                        # hvg_rfilter: fraction of genes kept
r_hvg_weight: 1.0                       # hvg_rweight: weight added to dispersions_norm
```

## Gene Selection Methods

`gene_selection` accepts any of the following. Which Snakemake rule handles a
value is decided by the `wildcard_constraints` in `rules/preprocess.smk`:
CPU-only arms go to `preprocess_cpu`, anything needing pyanglemania's
cross-batch statistics goes to `preprocess_gpu`.

**Baselines** (`preprocess_cpu`)

| Value | Method |
|-------|--------|
| `hvg` | CP10K + log1p, then `scanpy.pp.highly_variable_genes(flavor="seurat", batch_key=...)` |
| `full` | all genes |
| `rand` | random subset |
| `topbtvr` / `bottombtvr` | top / bottom genes by between-to-total variance ratio |

**anglemania** (`preprocess_gpu`)

| Value | Method |
|-------|--------|
| `angl` | `pyanglemania.preprocessing.anglemania` — ranks *gene pairs* by a 0.4/0.6 blend of mean and sd cross-batch z-score, takes genes on first appearance |

**Per-gene scores** (`preprocess_gpu`, `pipeline_scripts/gene_level_scores.py`)

These reduce the same cross-batch mean `M` and sd `S` that anglemania computes
to a *per-gene* score (row sums of `M²` and `S²`) instead of ranking pairs, and
z-score it within equal-frequency mean-expression bins:

| Value | Ranked on |
|-------|-----------|
| `anglgene` | binned signal/R blend |
| `anglgene_nobin` | same, without expression binning |
| `anglgene_anchor` / `anglgene_biology` | the two diagonals of the signal/R plane |
| `anglgene_noiselo` / `anglgene_noisehi` | the noise-residual axis, low / high end |
| `anglgene_mod`, `anglgene_anchor_mod`, `anglgene_biology_mod` | as above, with `S` empirical-Bayes moderated first |

Every `anglgene_*` arm is a different ranking of the *same* `M`/`S`, so
recomputing it per arm is wasted GPU time. Setting `anglgene_score_table_dir`
makes `preprocess_gpu` reuse one precomputed score table per sample instead.

**HVG × anglemania hybrids**

| Value | Method |
|-------|--------|
| `hvg_rfilter` | drop the worst `1 - r_keep_frac` of genes by binned R, then run plain HVG on the rest |
| `hvg_rweight` | rank on `dispersions_norm + r_hvg_weight × R_z_binned` |
| `hvg_intersect` | **control:** plain HVG restricted to anglemania's all-batch gene intersection, with no R involvement |

`hvg_intersect` is not optional when interpreting the other two: both of them
also silently drop every gene missing from some batch, which on real data is a
large fraction. Without this arm an R effect cannot be told apart from an
intersection effect.

> The `anglmean`, `anglmod05`, `anglmod1` and `anglsiggene` arms used by
> `config_mixed_effects_meanonly.yml` and `config_mixed_effects_altranking.yml`
> are **not** produced by `prepare_inputs.py`. Their gene lists are written
> directly into `preprocessed/` by standalone R scripts before the pipeline runs;
> the configs' header comments say which.

## Pipeline Overview

```
Input .h5ad files
       │
       ▼
┌─────────────────────┐
│  1. Preprocess      │  prepare_inputs.py (preprocess_cpu / preprocess_gpu)
│  Gene selection     │  + gene_level_scores.py for the anglgene_*/hvg_r* arms
│                     │  → preprocessed/{sample}_{gene_selection}.tsv
└─────────────────────┘
       │
       ▼
┌─────────────────────┐
│  2. Integrate       │  integration.py / seurat_integration.R
│  Batch correction   │  → embedding/{sample}/{method}/{sample}_{sel}.tsv
└─────────────────────┘
       │
       ▼
┌──────────────────────────────────────────────────┐
│  3. Metrics (four rules, run in parallel)        │
│   • metrics-scib_metrics.py → scib metrics       │
│   • metrics-bNMI.py         → balanced NMI       │
│   • metrics-cms.R           → CMS, BRAS          │
│   • metrics-ldfDiff.R       → ldfDiff            │
└──────────────────────────────────────────────────┘
       │
       ▼
┌─────────────────────┐
│  4. Combine         │  combine_metrics.py
│  Final scores       │  → metrics/{name}_combined_metrics.tsv
└─────────────────────┘
```

Stages 2–4 pass cell embeddings around as **TSV**, not `.h5ad`: `integration.py`
writes the embedding matrix, and each metrics script takes `--original_h5ad`
plus `--embedding_tsv`.

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
  utils/                  # R helpers sourced by scripts/R/*
  notebooks/              # .Rmd visualization and simulation notebooks
```

Only `pipeline_scripts/` is referenced from `rules/*.smk`. The `R/` and
`python/` scripts are run by hand — their header comments give the exact
invocation, and paths in them are relative to the repo root, so run them from
there.

## Output Structure

```
output_dir/
├── preprocessed/                          # {sample}_{gene_selection}.tsv  (gene lists)
├── embedding/
│   └── {sample}/{method}/                 # {sample}_{gene_selection}.tsv  (cell embeddings)
└── metrics/
    ├── {sample}/{method}/                 # *_scibmetrics.tsv, *_bnmi.tsv, *_cms.tsv, *_ldfDiff.tsv
    └── {name}_combined_metrics.tsv        # final benchmark results
```

The combined metrics TSV has columns: `sample`, `integration_method`,
`gene_selection`, all individual metrics, and summary scores (`Bio conservation`,
`Batch correction`, `Total`). `combine_metrics.py` recomputes the summary scores
itself rather than using scib's, so that the non-scib metrics (bNMI, ldfDiff,
CMS, BRAS) are included.

## Running Individual Scripts Manually

```bash
# Gene selection (in the pyanglemania conda env)
conda run -n pyanglemania python3 scripts/pipeline_scripts/prepare_inputs.py \
    --infile data/sample.h5ad --outfile out/genes.tsv \
    --batch_key Batch --gene_selection angl

# Integration (Python methods)
python3 scripts/pipeline_scripts/integration.py \
    --infile=data/sample.h5ad --feature_subset=out/genes.tsv \
    --outfile=out/embedding.tsv --integration_method=scvi \
    --batch_key=Batch --label_key=Group

# scib metrics
python3 scripts/pipeline_scripts/metrics-scib_metrics.py \
    --original_h5ad=data/sample.h5ad --embedding_tsv=out/embedding.tsv \
    --outfile=out/metrics.tsv --sample=mysample \
    --gene_selection=angl --integration_method=scvi \
    --batch_key=Batch --label_key=Group
```

## Analysis Notebooks

Notebooks are in `scripts/notebooks/`.

| Notebook | Purpose |
|----------|---------|
| `bm_simulated_visualization.Rmd` | Main benchmark results on simulated data |
| `bm_simulated_dropout_visualization.Rmd` | Benchmark results with dropout effects |
| `bm_LueckenTheisNatMethods2021_visualization.Rmd` | Benchmark on Luecken & Theis real datasets |
| `compare_gene_selections.Rmd` | Overlap and comparison of selected genes |
| `gene_selection_stability_mixed_effects.Rmd` | LOGO/LOBO/sparsity/composition stability of the gene sets |
| `group_composition_effect.Rmd` | Why `angl` is sensitive to group composition and `hvg` is not |
| `investigate_sim2_subbatch_anglemania.Rmd` | Sub-batch structure in the `sim2` simulation |
| `btvr_sanity_check_visualization.Rmd` | BTVR metric validation |
| `explore_anglemania_genes.Rmd` | Exploratory analysis of anglemania gene selections |
| `simulate_mixed_effects.Rmd` | Generate the mixed-effects simulation set |
| `simulate_single_cell.rmd` | Generate splatter-simulated datasets |
| `simulate_single_cell_with_dropout.rmd` | Generate dropout-simulated datasets |

## Key Dependencies

- [pyanglemania](https://github.com/omnideconv/pyanglemania) — GPU-accelerated Python port of anglemania, the gene selection method being benchmarked
- [anglemania](https://github.com/omnideconv/anglemania) — the original R implementation, used by the standalone analysis scripts
- [scanpy](https://scanpy.readthedocs.io/) — `seurat` (dispersion-based) batch-aware HVG selection
- [scib-metrics](https://github.com/YosefLab/scib-metrics) — batch integration benchmarking metrics
- [scvi-tools](https://github.com/scverse/scvi-tools) — scVI/scANVI integration
- [CellMixS](https://bioconductor.org/packages/CellMixS/) — CMS and ldfDiff metrics
- [balanced-clustering](https://github.com/scverse/balanced-clustering) — balanced NMI
- [Seurat v5](https://satijalab.org/seurat/) — Seurat-based integration method
- [Snakemake](https://snakemake.readthedocs.io/) — workflow management
