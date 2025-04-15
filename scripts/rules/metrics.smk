rule metrics:
    input:
        original_h5ad=lambda wildcards: file_paths[wildcards.sample],
        integrated_h5ad = rules.integrate.output
    output:
        metrics = join(out_metrics, "{sample}/{integration_method}/{sample}_{gene_selection}.tsv")
    params:
        out_metrics = out_metrics,
        batch_key = config["batch_key"],
        label_key = config["label_key"]
    # resources:
    #     # Dynamically assign GPU resource based on integration method
    #     gpu=lambda wildcards: 1 if wildcards.integration_method in ["scvi", "scanvi"] else 0
    shell:
        """
        python3 pipeline_scripts/metrics.py \
        --original_h5ad={input.original_h5ad} \
        --integrated_h5ad={input.integrated_h5ad} \
        --outfile={output.metrics} \
        --sample={wildcards.sample} \
        --gene_selection={wildcards.gene_selection} \
        --integration_method={wildcards.integration_method} \
        --batch_key={params.batch_key} \
        --label_key={params.label_key}
        """

rule combine:
    input:
        all_metrics = expand(rules.metrics.output, sample=SAMPLES, integration_method=INTEGRATION_METHODS, gene_selection=GENE_SELECTIONS)
    output:
        combined = join(out_metrics, config["name"] + "_combined_metrics.tsv")
    run:
        import os
        import pandas as pd
        files = input["all_metrics"]
        df_list = []
        for file in files:
            metrics = pd.read_csv(file, sep = "\t")
            # print(metrics)
            # Ensure there are no duplicated column names
            metrics = metrics.loc[:, ~metrics.columns.duplicated()]
            metrics = metrics.loc[:, ~metrics.columns.isin(["Unnamed: 0","cell_cycle_conservation","hvg_overlap", "trajectory"])] # old scib, not scib_metrics
            df_list.append(metrics)
        print(df_list[0])
        combined_df = pd.concat(df_list, axis = 0)
        print(combined_df.head())
        if os.path.exists(output[0]):
            # load previous metrics df if it exists, in order to update with new metrics and/or overwrite the samples that are already in there
            # it's like a custom append function
            metrics_df = pd.read_csv(output[0], sep = "\t")
            combined_df = pd.concat([metrics_df[~metrics_df["sample"].isin(combined_df["sample"])], combined_df], axis = 0)
        combined_df.to_csv(output[0], sep = "\t", index = False)
            
            