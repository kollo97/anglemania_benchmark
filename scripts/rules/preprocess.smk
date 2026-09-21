rule preprocess_cpu:
    input:
        original_h5ad=lambda wildcards: file_paths[wildcards.sample]
    output:
        outfile = PREPROCESS_OUTPUT
    wildcard_constraints:
        gene_selection = "hvg|hvg_intersect|full|rand|topbtvr|bottombtvr"
    params:
        batch_key = config["batch_key"],
        anglemania_mode = config["anglemania_mode"],
        n_genes = config["n_genes"],
        permutation_function = config["permutation_function"],
        normalization_mode = config["normalization_mode"]
    resources:
        cpus_per_task = 2,
        mem_mb = 16000,
        runtime = 60
    shell:
        """
        conda run -n pyanglemania python3 pipeline_scripts/prepare_inputs.py \
            --infile {input.original_h5ad} \
            --outfile {output.outfile} \
            --batch_key {params.batch_key} \
            --gene_selection {wildcards.gene_selection} \
            --anglemania_mode {params.anglemania_mode} \
            --n_genes {params.n_genes} \
            --permutation_function {params.permutation_function} \
            --normalization_mode {params.normalization_mode}
        """

rule preprocess_gpu:
    input:
        original_h5ad=lambda wildcards: file_paths[wildcards.sample]
    output:
        outfile = PREPROCESS_OUTPUT
    wildcard_constraints:
        # longest alternative first, so `anglgene` cannot shadow `anglgene_anchor`
        gene_selection = "anglgene_biology_mod|anglgene_anchor_mod|anglgene_noisehi|anglgene_biology|anglgene_noiselo|anglgene_anchor|anglgene_nobin|anglgene_mod|hvg_rfilter|hvg_rweight|anglgene|angl"
    params:
        batch_key = config["batch_key"],
        anglemania_mode = config["anglemania_mode"],
        n_genes = config["n_genes"],
        permutation_function = config["permutation_function"],
        normalization_mode = config["normalization_mode"],
        # anglgene/anglgene_nobin only; ignored by the pair-ranking "angl" path
        n_bins = config.get("anglgene_n_bins", 20),
        weight_signal = config.get("anglgene_weight_signal", 0.5),
        weight_R = config.get("anglgene_weight_R", 0.5),
        # hvg_rfilter / hvg_rweight only
        r_keep_frac = config.get("r_keep_frac", 0.5),
        r_hvg_weight = config.get("r_hvg_weight", 1.0),
        # anglgene_* only: the two residual axes and the moderated-S prior trend
        loess_span = config.get("anglgene_loess_span", 0.5),
        moderate_n_bins = config.get("anglgene_moderate_n_bins", 10),
        # Every anglgene_* arm is a different ranking of the SAME cross-batch
        # M/S, so point them at one precomputed score table instead of paying
        # for the GPU pass once per arm. Unset = recompute per arm.
        score_table = lambda wildcards: (
            f"--score_table {config['anglgene_score_table_dir']}/"
            f"{wildcards.sample}_gene_level_scores.tsv"
            if config.get("anglgene_score_table_dir")
            and wildcards.gene_selection.startswith("anglgene")
            else ""
        )
    resources:
        cpus_per_task = 2,
        mem_mb = 16000,
        runtime = 120,
        gres = "gpu:1",
        gpu_slots = 1
    shell:
        """
        conda run -n pyanglemania python3 pipeline_scripts/prepare_inputs.py \
            --infile {input.original_h5ad} \
            --outfile {output.outfile} \
            --batch_key {params.batch_key} \
            --gene_selection {wildcards.gene_selection} \
            --anglemania_mode {params.anglemania_mode} \
            --n_genes {params.n_genes} \
            --permutation_function {params.permutation_function} \
            --normalization_mode {params.normalization_mode} \
            --n_bins {params.n_bins} \
            --weight_signal {params.weight_signal} \
            --weight_R {params.weight_R} \
            --r_keep_frac {params.r_keep_frac} \
            --r_hvg_weight {params.r_hvg_weight} \
            --loess_span {params.loess_span} \
            --moderate_n_bins {params.moderate_n_bins} {params.score_table}
        """
