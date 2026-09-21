rule integrate:
    input:
        original_h5ad=lambda wildcards: file_paths[wildcards.sample],
        feature_subset = PREPROCESS_OUTPUT
    output:
        embedding = join(out_embedding, "{sample}/{integration_method}/{sample}_{gene_selection}.tsv")
    params:
        out_integration = out_integration,
        batch_key = config["batch_key"],
        label_key = config["label_key"]
    resources:
        cpus_per_task=4,
        mem_mb = 32000,
        runtime = 360
    shell:
        """
        if [[ "{wildcards.integration_method}" = "harmony" || "{wildcards.integration_method}" = "scvi" || "{wildcards.integration_method}" = "scanvi" || "{wildcards.integration_method}" = "scanorama" ]]; then
            python3 pipeline_scripts/integration.py \
                --infile={input.original_h5ad} \
                --feature_subset={input.feature_subset} \
                --outfile={output.embedding} \
                --integration_method={wildcards.integration_method} \
                --batch_key={params.batch_key} \
                --label_key={params.label_key}
        elif [[ "{wildcards.integration_method}" = "seurat" ]]; then
            Rscript pipeline_scripts/seurat_integration.R --infile={input.original_h5ad} \
                                        --feature_subset={input.feature_subset} \
                                        --outfile {output.embedding} \
                                        --batch_key {params.batch_key} \
                                        --label_key {params.label_key}
        fi

        """