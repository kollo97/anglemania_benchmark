rule integrate:
    input:
        original_h5ad=lambda wildcards: file_paths[wildcards.sample],
        feature_subset = rules.preprocess.output
    output:
        integrated = join(out_integration, "{sample}/{integration_method}/{sample}_{gene_selection}.h5ad")
    params:
        out_integration = out_integration,
        batch_key = config["batch_key"],
        label_key = config["label_key"]
    resources:
        # Dynamically assign GPU resource based on integration method
        gpu=lambda wildcards: 1 if wildcards.integration_method in ["scvi", "scanvi"] else 0
    shell:
        """
        export JAX_PLATFORMS=cpu
        if [[ "{wildcards.integration_method}" = "harmony" || "{wildcards.integration_method}" = "scvi" || "{wildcards.integration_method}" = "scanvi" ]]; then
            python3 integration.py \
                --infile={input.original_h5ad} \
                --feature_subset={input.feature_subset} \
                --outfile={output.integrated} \
                --integration_method={wildcards.integration_method} \
                --batch_key={params.batch_key} \
                --label_key={params.label_key}
        elif [[ "{wildcards.integration_method}" = "seurat" ]]; then
            Rscript seurat_integration.R --infile={input.original_h5ad} \
                                        --feature_subset={input.feature_subset} \
                                        --outfile {output.integrated} \
                                        --batch_key {params.batch_key} \
                                        --label_key {params.label_key}
        fi

        """