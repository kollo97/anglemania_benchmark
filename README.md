# anglemania Benchmark

Benchmarking pipeline for [anglemania](https://github.com/omnideconv/anglemania), an R package for gene selection in single-cell RNA-seq batch integration. The pipeline systematically compares **anglemania** (`angl`) against **highly variable gene** (`hvg`) selection across multiple integration methods and datasets, using a comprehensive set of batch correction and biological conservation metrics.

## Experimental Design

| Dimension | Options |
|-----------|---------|
| Gene selection | `angl` (anglemania), `hvg` (Seurat HVGs) |
| Integration methods | `scvi`, `scanorama` |
| Datasets | Splatter-simulated (varying batch/bio effects, dropout), NeurIPS 2021 CITE-seq, Luecken & Theis NatMethods 2021 |
| Metrics | Silhouette label, Isolated labels, bNMI, cLISI, ldfDiff, Graph connectivity, PCR comparison, iLISI, CMS, BRAS |
| Final score | **Total = 0.6 × Bio conservation + 0.4 × Batch correction** |

## Environment Setup

All commands must be run inside the **Guix environment** (recommended for reproducibility):

```bash
# Via shell alias (defined in ~/.bashrc):
bm_guix
# Equivalent to:
guix time-machine -C guix/channels.scm -- shell -m guix/manifest.scm --no-grafts python-numpyro
```

## Usage

All pipeline commands run from the `scripts/` directory:

```bash
cd scripts/

# Dry run — check what will be executed
snakemake -s Snakefile.smk --configfile config_files/<config>.yml -n

# Local run
snakemake -s Snakefile.smk --configfile config_files/<config>.yml --cores 8
```

### Available Configs

| Config file | Dataset |
|-------------|---------|
| `config_varying_batch_biovar_facScale0.1.yml` | Splatter-simulated, moderate biological difference |
| `config_varying_batch_biovar_facScale0.2.yml` | Splatter-simulated, strong biological difference in all Sim |
| `config_varying_batch_biovar_dropout_sample.yml` | Splatter-simulated with dropout |
| `config_btvr_sanity_check.yml` | BTVR sanity check - is high BTVR good? (simulated) |
| `config_n_cells_test.yml` | Cell number variation (simulated) |
| `LueckenTheis_config.yml` | Luecken, Theis NatMethods 2021 |
| `NeurIPS_CITE-seq.yml` | NeurIPS 2021 CITE-seq (full) |
| `NeurIPS_CITE-seq_Tcells.yml` | NeurIPS 2021 CITE-seq (T cells) |
| `NeurIPS_CITE-seq_filtered_no_ribo.yml` | NeurIPS 2021 CITE-seq without ribosomal genes |
| `Theis_dropout_simulations.yml` | Simulations based on Luecken, Theis Splatter parameters|

### Key Config Parameters

```yaml
name: my_experiment
samplesheet: path/to/samplesheet.tsv   # columns: sample_name, file_path (.h5ad)
output_dir: path/to/output/
batch_key: Batch                        # obs column for batch labels
label_key: Group                        # obs column for cell type labels
integration_methods: [scvi, scanorama]
gene_selection: [hvg, angl]
anglemania_mode: cosine                 # cosine | spearman | diem
permutation_function: sample            # sample | permute_nonzero
n_genes: 2000
```

## Pipeline Overview

```
Input .h5ad files
       │
       ▼
┌─────────────────────┐
│  1. Preprocess      │  prepare_inputs.R
│  Gene selection     │  → {sample}_{gene_selection}.tsv
└─────────────────────┘
       │
       ▼
┌─────────────────────┐
│  2. Integrate       │  integration.py / seurat_integration.R
│  Batch correction   │  → embedding/{sample}/{method}/*.tsv
└─────────────────────┘
       │
       ▼
┌──────────────────────────────────────────┐
│  3. Metrics (run in parallel)            │
│   • scib_metrics.py  → scib metrics      │
│   • bNMI.py          → balanced NMI      │
│   • cms.R            → CMS, BRAS         │
│   • ldfDiff.R        → ldfDiff           │
└──────────────────────────────────────────┘
       │
       ▼
┌─────────────────────┐
│  4. Combine         │  combine_metrics.py
│  Final scores       │  → *_combined_metrics.tsv
└─────────────────────┘
```

## Output Structure

```
output_dir/
├── preprocessed/                          # {sample}_{gene_selection}.tsv  (gene lists)
├── embedding/
│   └── {sample}/{method}/                 # {sample}_{gene_selection}.tsv  (cell embeddings)
└── metrics/
    ├── {sample}/{method}/                 # per-sample metric TSVs
    └── {name}_combined_metrics.tsv        # final benchmark results
```

The combined metrics TSV has columns: `sample`, `integration_method`, `gene_selection`, all individual metrics, and summary scores (`Bio conservation`, `Batch correction`, `Total`).

## Analysis Notebooks

Notebooks are in `scripts/notebooks/` and generate HTML reports and plots under `output/plots/`.

| Notebook | Purpose |
|----------|---------|
| `bm_simulated_visualization.Rmd` | Main benchmark results on simulated data |
| `bm_simulated_visualization_dropout.Rmd` | Benchmark results with dropout effects |
| `bm_LueckenTheisNatMethods2021_visualization.Rmd` | Benchmark on Luecken & Theis real datasets |
| `compare_gene_selections.Rmd` | Overlap and comparison of selected genes |
| `btvr_sanity_check_visualization.Rmd` | BTVR metric validation |
| `simulate_single_cell.rmd` | Generate splatter-simulated datasets |
| `simulate_single_cell_with_dropout.rmd` | Generate dropout-simulated datasets |
| `explore_anglemania_genes.Rmd` | Exploratory analysis of anglemania gene selections |

## Key Dependencies

- [anglemania](https://github.com/omnideconv/anglemania) — the gene selection method being benchmarked
- [scib-metrics](https://github.com/YosefLab/scib-metrics) — batch integration benchmarking metrics
- [scvi-tools](https://github.com/scverse/scvi-tools) — scVI/scANVI integration
- [CellMixS](https://bioconductor.org/packages/CellMixS/) — CMS and ldfDiff metrics
- [balanced-clustering](https://github.com/scverse/balanced-clustering) — balanced NMI
- [Seurat v5](https://satijalab.org/seurat/) — HVG selection and integration
- [Snakemake](https://snakemake.readthedocs.io/) — workflow management
