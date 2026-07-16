rule preprocess_cpu:
    input:
        original_h5ad=lambda wildcards: file_paths[wildcards.sample]
    output:
        outfile = PREPROCESS_OUTPUT
    wildcard_constraints:
        gene_selection = "hvg|full|rand|topbtvr|bottombtvr"
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
        gene_selection = "angl"
    params:
        batch_key = config["batch_key"],
        anglemania_mode = config["anglemania_mode"],
        n_genes = config["n_genes"],
        permutation_function = config["permutation_function"],
        normalization_mode = config["normalization_mode"]
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
            --normalization_mode {params.normalization_mode}
        """
