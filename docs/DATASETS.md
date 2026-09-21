# Datasets Used in the Benchmark Pipeline

Inventory of the datasets the Snakemake pipeline reads. Datasets enter the
pipeline only through a config's `samplesheet:` field
(`scripts/config_files/*.yml`), so this document is organised by samplesheet.

Sizes measured with `du` on 2026-09-17; dimensions read from the `.h5ad` files.

---

## 1. Main Datasets

The set the benchmark actually runs on going forward: the two real-data
benchmarks, and three slices of the splatter simulations (the `de.facLoc 0.1`
signal level, the mixed-effects stability set, and the dropout scenarios).

**Disk footprint — `.h5ad` datasets only:**

| # | Group | Samples | Size |
|---|---|---:|---:|
| 1.1 | Luecken & Theis 2021 benchmark | 6 | **17 G** |
| 1.2 | NeurIPS 2021 CITE-seq BMMC | 9 | **7.4 G** |
| 1.3 | Splatter `de.facLoc 0.1` — `facScale0.1` | 54 | **28 G** |
| 1.3 | Splatter `de.facLoc 0.1` — `facScale0.2` | 54 | **28 G** |
| 1.4 | Splatter mixed-effects stability | 11 + baseline | **5.7 G** |
| 1.5 | Splatter dropout — own sim | 54 | **41 G** |
| 1.5 | Splatter dropout — Theis-parameterised | 20 | **23 G** |
| | **Total** | **209** | **≈ 150 G** |

Only `facScale0.1` is used by the active configs; drop the `facScale0.2` slice
and the total is **≈ 122 G**. Derived artefacts (gene-pair stat tables,
pipeline outputs, SEACells) are listed separately in §4 and are not counted above.

---

### 1.1 Luecken & Theis, Nature Methods 2021 — benchmark datasets

The field-standard integration benchmark datasets. **17 G** (6 files).

- **Directory:** `/local/Projects/AAkalin_Neuroblastoma/data/Luecken_Theis_NatMethods_2021_benchmarking/datasets/`
- **Samplesheet:** `datasets/samplesheet.tsv`
- **Configs:** `LueckenTheis_config.yml` (all 6), `config_genelevel_pancreas.yml` and `config_genelevel_pancreas_cosine_diag.yml` (pancreas only, via `samplesheet_genelevel_pancreas.tsv`)
- **Keys:** `batch_key: batch`, `label_key: celltype`
- **Raw counts:** in `layers["counts"]`; `X` is normalised

| Sample | Cells | Genes | Batches | Cell types | Size | Notes |
|---|---:|---:|---:|---:|---:|---|
| `human_pancreas_norm_complexBatch` | 16 382 | 19 093 | 9 | 14 | 3.5 G | Human pancreas; batches are protocols (`celseq`, `celseq2`, `fluidigmc1`, `inDrop1-4`, `smarter`, `smartseq2`). The main real-data workhorse for the per-gene score work. |
| `Immune_ALL_human` | 33 506 | 12 303 | 10 | 16 | 4.7 G | Human immune cells across studies (`10X`, `Freytag`, `Oetjen_A/P/U`, `Sun_sample1-4`, `Villani`). Labels also in `final_annotation`. |
| `Lung_atlas_public` | 32 472 | 15 148 | 16 | 17 | 1.5 G | Human lung atlas, donor/site batches. Labels also in `cell_type`. |
| `sim1_raw` | 12 097 | 9 979 | 6 | 7 | 1.6 G | Splatter simulation shipped with the paper. Also has `Batch`/`Group`. |
| `sim2_raw` | 19 318 | 10 000 | 4 | 4 | 3.0 G | As above; `obs` additionally has `Sub`/`SubBatch`. |
| `sim2_subbatch_raw` | 19 318 | 10 000 | 16 | 4 | 3.0 G | Nested sub-batch variant: `batch` = `Batch{1-4}Sub{1-4}`, `Batch` still 4 levels. Investigated in `scripts/notebooks/investigate_sim2_subbatch_anglemania.Rmd`. |

### 1.2 NeurIPS 2021 Open Problems — CITE-seq bone marrow mononuclear cells

GEO accession **GSE194122**. Multi-site, multi-donor CITE-seq; the GEX and ADT
modalities are also split into separate objects so each can be integrated on its
own. **7.4 G** across the 9 objects referenced by samplesheets.

- **Directory:** `/data/akalin/akollot/anglemania_benchmark/data/NeurIPS_2021/CITE-seq_CMBC/`
- **Keys:** `batch_key: batch` (12 site×donor combinations `s1d1`…`s4d9`), `label_key: cell_type`
- **Other useful `obs` columns:** `Site` (4), `DonorID` (9), `Samplename`, `Modality`, `is_train`
- **Raw counts:** in `layers["counts"]`; precomputed `obsm`: `GEX_X_pca`, `GEX_X_umap`, `ADT_X_pca`, `ADT_X_umap`, `ADT_isotype_controls`

**a) Full object** — `NeurIPS_CITE-seq.yml`, samplesheet `samplesheet_no_Tcells.tsv` — **2.7 G**

| Sample | Cells | Features | Cell types | Size |
|---|---:|---:|---:|---:|
| `GSE194122_openproblems_neurips2021_cite_BMMC_clean` | 90 261 | 14 087 (GEX + ADT) | 45 | 1.4 G |
| `cite_BMMC_gex` | 90 261 | 13 953 genes | 45 | 1.3 G |
| `cite_BMMC_adt` | 90 261 | 134 antibodies | 45 | 112 M |

**b) T-cell subset** — `NeurIPS_CITE-seq_Tcells.yml`, samplesheet `samplesheet_Tcells.tsv` — **601 M**

| Sample | Cells | Features | Cell types | Size |
|---|---:|---:|---:|---:|
| `cite_BMMC_Tcells` | 24 675 | 14 087 | 13 | 300 M |
| `cite_BMMC_Tcells_gex` | 24 675 | 13 953 | 13 | 271 M |
| `cite_BMMC_Tcells_adt` | 24 675 | 134 | 13 | 31 M |

A deliberately hard case: all labels are closely related T-cell subtypes, so
biological conservation has to resolve fine-grained structure.

**c) Custom-filtered, ribosomal genes removed** — `NeurIPS_CITE-seq_filtered_no_ribo.yml`, samplesheet `filtered_custom/samplesheet.tsv` — **4.1 G**

| Sample | Cells | Features | Cell types | Size |
|---|---:|---:|---:|---:|
| `GSE194122_openproblems_neurips2021_cite_BMMC_filtered` | 89 744 | 13 968 | 45 | 1.8 G |
| `cite_BMMC_gex_filtered` | 89 744 | 13 834 | 45 | 1.7 G |
| `cite_BMMC_adt_filtered` | 89 744 | 134 | 45 | 696 M |

QC metrics added in `obs` (`pct_counts_mt`, `pct_counts_ribo`, `n_genes_by_counts`, …).

### 1.3 Splatter — the `de.facLoc 0.1` signal level

The moderate-biological-signal slice of the batch × biological-variance grids:
`de.facLoc` fixed at 0.1, swept over `batch.facLoc` ∈ {0, 0.01, 0.05, 0.1, 0.3, 0.6}
× `nbatch` ∈ {2, 4, 5} × `ngroup` ∈ {2, 4, 5}, at `groupCells` 1000 and
`nGenes` 10 000 — 54 files per `de.facScale` condition.

| Condition | Directory | Files | Size | Config |
|---|---|---:|---:|---|
| `de.facScale 0.1` (mild) | `h5ad/facScale0.1/*de.facLoc0.1_*groupCells1000` | 54 | **28 G** | `config_varying_batch_biovar_facScale0.1.yml`, `..._facScale0.1_pyanglemania.yml` |
| `de.facScale 0.2` (stronger) | `h5ad/facScale0.2/*de.facLoc0.1_*groupCells1000` | 54 | **28 G** | `config_varying_batch_biovar_facScale0.2.yml` |

Base path: `/local/Projects/AAkalin_Neuroblastoma/Results/simulate_single_cell/`.
Keys: `batch_key: Batch`, `label_key: Group`. Individual files run 155 MB
(`nbatch2_ngroup2`) to 959 MB (`nbatch5_ngroup5`).

**Note:** each `facScale` directory also holds a `groupCells500` copy of the
whole grid (54 more files, 14 G, per `de.facLoc` level), but no samplesheet
references it. Both the `groupCells500` copies and the other `de.facLoc` levels
were pruned from disk; see §4.

Three files from this slice also drive the batch-strength sweep for the per-gene
score work (`config_genelevel_batchsweep.yml`, samplesheet
`samplesheet_genelevel_batchsweep.tsv`): `batch.facLoc` ∈ {0.1, 0.3, 0.6} at
`nbatch4_ngroup4_groupCells1000`, 615 MB each.

### 1.4 Splatter — mixed-effects stability set

The set behind `docs/ANGLEMANIA_INSTABILITY_ANALYSIS.md`. One baseline
simulation plus 11 perturbations, so gene-set stability can be measured along
four axes: leave-one-group-out (LOGO), leave-one-batch-out (LOBO), sparsity, and
composition. **5.7 G** (5.1 G for the 11 perturbations + 615 M baseline).

- **Baseline:** `batch.facLoc0.3_de.facLoc0.1_nbatch4_ngroup4_groupCells1000` — 16 000 cells × 10 000 genes, 4 batches, 4 groups, in `h5ad/facScale0.1/` (615 M)
- **Perturbations:** `h5ad/mixed_effects/`, all suffixed onto the baseline name

| Suffix | Axis | n | Size each |
|---|---|---:|---:|
| `_dropGroupGroup{1-4}` | LOGO — one cell type removed | 4 | ~460 M |
| `_dropBatchBatch{1-4}` | LOBO — one batch removed | 4 | ~460 M |
| `_dropout.mid3.5_dropout.shape-1` | sparsity | 1 | 1.2 G |
| `_composition_skewed10_20_30_40` | composition — skewed group proportions | 1 | 155 M |
| `_composition_rareGroup1_1pct` | composition — one rare (1 %) group | 1 | 155 M |

| Config | Samplesheet | Samples | Gene selections |
|---|---|---:|---|
| `config_mixed_effects.yml` | `h5ad/mixed_effects/samplesheet.tsv` | 12 (baseline + 11) | `hvg`, `angl` |
| `config_mixed_effects_meanonly.yml` | `samplesheet_meanonly_lobo.tsv` | 5 (baseline + 4 LOBO) | `anglmean` |
| `config_mixed_effects_altranking.yml` | `samplesheet_altranking_lobo.tsv` | 5 (baseline + 4 LOBO) | `anglmod1`, `anglmod05`, `anglsiggene` |
| `config_genelevel_sim.yml` | `samplesheet_genelevel_sim.tsv` | 1 (baseline) | the 11 `anglgene_*` per-gene arms |

The last three configs re-use the same `.h5ad` files; their gene lists are
written directly into `preprocessed/` by helper scripts in `scripts/R/` rather
than by `prepare_inputs.py`.

### 1.5 Splatter — dropout scenarios

Simulations with logistic dropout, for testing gene selection under realistic
sparsity. **64 G** total.

| Set | Directory | Files | Size | Config / samplesheet |
|---|---|---:|---:|---|
| Own grid | `h5ad/with_dropout/own_sim/` | 54 | **41 G** | `config_varying_batch_biovar_dropout_sample.yml` / `own_sim/samplesheet.tsv` |
| Theis-parameterised | `h5ad/with_dropout/Theis_sim/` | 20 | **23 G** | `Theis_dropout_simulations.yml` / `Theis_sim/samplesheet.json` |

**Own grid** (54 files, ~460 MB each): `batch.facLoc` ∈ {0.05, 0.1, 0.3} ×
`de.facLoc` ∈ {0.05, 0.1, 0.3} × `dropout.mid` ∈ {1, 1.5, 2, 2.5, 3, 3.5}, at
`dropout.shape` −1, `nbatch` 4, `ngroup` 3, 1000 cells/group.

**Theis-parameterised** (20 files, 530 MB – 1.5 GB): splatter runs parameterised
to match the Luecken & Theis `sim1`/`sim2` designs, then swept over dropout
severity. Uses a **JSON samplesheet**, which carries the full per-sample splatter
parameter vectors (`batchCells`, `batch.facLoc`, `batch.facScale`, `de.prob`,
`group.prob`, `de.facLoc`, `de.facScale`, `sparsity`, `lib.loc`, `seed`).

- `sim1`: 6 batches × 7 groups, 10 000 cells/batch, 12 097-cell scale
- `sim2`: 4 batches × 4 groups
- `dropout.mid` sweeps 0.5 → 3.5 in 10 even steps, at `dropout.shape` −5, `dropout.type` `experiment`

---

## 2. Splatter Parameter Naming

Simulated sample names encode their parameters, e.g.
`batch.facLoc0.3_de.facLoc0.1_nbatch4_ngroup4_groupCells1000`:

| Parameter | Meaning |
|---|---|
| `batch.facLoc` | batch-effect location factor — batch effect strength |
| `batch.facScale` | batch-effect scale factor |
| `de.facLoc` | DE location factor — biological (between-group) signal strength |
| `de.facScale` | DE scale factor |
| `nbatch` / `ngroup` | number of batches / cell-type groups |
| `groupCells` | cells per group |
| `dropout.mid` / `dropout.shape` | logistic dropout midpoint / steepness |

All own simulations use `batch_key: Batch`, `label_key: Group`, and are generated
by the notebooks in `scripts/notebooks/` (see `simulate_mixed_effects.Rmd`).

---

## 3. Notes

- **Samplesheet format:** TSV/CSV with at least `sample_name` and `file_path`;
  extra columns (the splatter parameters) are carried along and ignored by the
  pipeline. JSON is also accepted (`Theis_dropout_simulations.yml`).
- **Active comparison:** `hvg` vs `angl` gene selection with `scvi` and
  `scanorama` integration. `full` / `rand` are commented out in the configs.
- The Luecken & Theis and NeurIPS objects keep raw counts in `layers["counts"]`
  while `X` holds normalised values, which is what `prepare_inputs.py` and the
  metrics scripts expect.

---

## 4. Pruned and Retained Non-Main Data

### 4.1 Pruned simulations

The simulation tree was pruned to the main set: **1197 `.h5ad` files, 438 G
removed**, leaving the 193 files listed in §1. Run via
`scratchpad/prune_simulations.sh`.

| Removed | Files | Size |
|---|---:|---:|
| `h5ad/facScale0.1/` — other `de.facLoc` levels, `groupCells1000` | 270 | 137 G |
| `h5ad/facScale0.1/` — all `groupCells500` | 324 | 82 G |
| `h5ad/facScale0.2/` — other `de.facLoc` levels, `groupCells1000` | 270 | 137 G |
| `h5ad/facScale0.2/` — all `groupCells500` | 324 | 82 G |
| `simulate_single_cell/test_numcells/h5ad/` — cell-number sweep | 9 | 1.7 G |
| **Total** | **1197** | **438 G** |

Each `facScale` directory originally held 648 files (246 G): a full 6 × 6
`batch.facLoc` × `de.facLoc` grid × 9 `nbatch`/`ngroup` combos × 2 `groupCells`
levels. The main set keeps 54 per directory. The `groupCells500` copies were
referenced by no samplesheet at all.

**Configs removed with the data:**

| Config | Why |
|---|---|
| `config_btvr_sanity_check.yml` | needed all 6 `de.facLoc` levels at `nbatch4_ngroup4` (36 files) |
| `config_n_cells_test.yml` | its `test_numcells/` dataset was pruned |
| `config_varying_batch_biovar_facScale0.1_pyanglemania.yml` | a subset of the `facScale0.1` grid; redundant once that config was trimmed to the same slice |
| `snakemake_config.yml` | already dead — pointed at an 800-sample grid whose `.h5ad` files no longer existed |

`config_varying_batch_biovar_facScale0.1.yml` and `..._facScale0.2.yml` were
retargeted from the 324-row grid samplesheets to
`samplesheet_facScale{0.1,0.2}_deLoc0.1.tsv` (54 rows each), so they now run
exactly the §1.3 main slice.

Three notebooks still reference the removed configs and will not run as-is:
`notebooks/btvr_sanity_check_visualization.Rmd`,
`notebooks/bm_simulated_visualization_n_cells_test.Rmd`, and
`notebooks/bm_simulated_visualization.Rmd` (reads the `_pyanglemania` config).

### 4.2 Retained, not counted as datasets

Derived artefacts and archives left in place. These now dominate the remaining
footprint — the `Theis_sim` support directories alone are 354 G, far more than
the 23 G of dropout datasets they were built from.

| Path | Size | What |
|---|---:|---|
| `h5ad/with_dropout/Theis_sim/SEACells/` | 284 G | SEACells metacell outputs from an earlier line of work |
| `h5ad/with_dropout/Theis_sim/tmp_archive_sim2/` | 61 G | temporary archive |
| `h5ad/with_dropout/Theis_sim/preprocessed/` | 9 G | gene lists / intermediates |
| `h5ad/mixed_effects/anglemania_stats/` | 16 G | per-gene-pair stat tables from `scripts/python/stability_anglemania_stats.py`; **still read** by the `scripts/R/stability_*.R` scripts |
| `h5ad/mixed_effects/pipeline_output/` | 740 M | pipeline results for the mixed-effects + genelevel runs |
| `h5ad/mixed_effects/group_composition_effect/` | 68 M | tables from `stability_group_composition_effect.py` |
| `Luecken_Theis_.../datasets/all_datasets.zip` | 19 G | original download archive; the 6 `.h5ad` are already extracted |
| `NeurIPS_2021/.../GSE194122_..._cite_BMMC.h5ad` | 1.4 G | pre-cleaning original, superseded by the `_clean` object; in no samplesheet |

Everything above except `anglemania_stats/` is a candidate for removal if more
space is needed; `SEACells/` + `tmp_archive_sim2/` alone would free 345 G.
