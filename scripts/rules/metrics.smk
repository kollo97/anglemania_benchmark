rule scib_metrics:
    input:
        original_h5ad=lambda wildcards: file_paths[wildcards.sample],
        embedding_tsv = rules.integrate.output.embedding
    output:
        metrics = join(out_metrics, "{sample}/{integration_method}/{sample}_{gene_selection}_scibmetrics.tsv")
    params:
        out_metrics = out_metrics,
        batch_key = config["batch_key"],
        label_key = config["label_key"]
    # resources:
    #     # Dynamically assign GPU resource based on integration method
    #     gpu=lambda wildcards: 1 if wildcards.integration_method in ["scvi", "scanvi"] else 0
    shell:
        """
        python3 pipeline_scripts/metrics-scib_metrics.py \
            --original_h5ad={input.original_h5ad} \
            --embedding_tsv={input.embedding_tsv} \
            --outfile={output.metrics} \
            --sample={wildcards.sample} \
            --gene_selection={wildcards.gene_selection} \
            --integration_method={wildcards.integration_method} \
            --batch_key={params.batch_key} \
            --label_key={params.label_key}
        """

rule bnmi:
    input:
        original_h5ad=lambda wildcards: file_paths[wildcards.sample],
        embedding_tsv = rules.integrate.output.embedding
    output:
        bnmi = join(out_metrics, "{sample}/{integration_method}/{sample}_{gene_selection}_bnmi.tsv")
    params:
        out_metrics = out_metrics,
        label_key = config["label_key"]
    shell:
        """
        python3 pipeline_scripts/metrics-bNMI.py \
            --original_h5ad={input.original_h5ad} \
            --embedding_tsv={input.embedding_tsv} \
            --outfile={output.bnmi} \
            --sample={wildcards.sample} \
            --gene_selection={wildcards.gene_selection} \
            --integration_method={wildcards.integration_method} \
            --label_key={params.label_key}
        """

rule cms:
    input:
        original_h5ad=lambda wildcards: file_paths[wildcards.sample],
        embedding_tsv = rules.integrate.output.embedding
    output:
        cms = join(out_metrics, "{sample}/{integration_method}/{sample}_{gene_selection}_cms.tsv")
    params:
        out_metrics = out_metrics,
        batch_key = config["batch_key"]
    shell:
        """
        Rscript pipeline_scripts/metrics-cms.R \
            --original_h5ad {input.original_h5ad} \
            --embedding_tsv {input.embedding_tsv} \
            --batch_key {params.batch_key} \
            --outfile {output.cms} \
            --sample {wildcards.sample} \
            --gene_selection {wildcards.gene_selection} \
            --integration_method {wildcards.integration_method}
        """

rule ldfdiff:
    input:
        original_h5ad=lambda wildcards: file_paths[wildcards.sample],
        embedding_tsv = rules.integrate.output.embedding
    output:
        ldfdiff = join(out_metrics, "{sample}/{integration_method}/{sample}_{gene_selection}_ldfDiff.tsv")
    params:
        out_metrics = out_metrics,
        batch_key = config["batch_key"],
        label_key = config["label_key"]
    shell:
        """
        Rscript pipeline_scripts/metrics-ldfDiff.R \
            --original_h5ad {input.original_h5ad} \
            --embedding_tsv {input.embedding_tsv} \
            --batch_key {params.batch_key} \
            --label_key {params.label_key} \
            --outfile {output.ldfdiff} \
            --sample {wildcards.sample} \
            --gene_selection {wildcards.gene_selection} \
            --integration_method {wildcards.integration_method}
        """

rule combine:
    input:
        scib_metrics = expand(rules.scib_metrics.output, sample=SAMPLES, integration_method=INTEGRATION_METHODS, gene_selection=GENE_SELECTIONS),
        bnmi = expand(rules.bnmi.output, sample=SAMPLES, integration_method=INTEGRATION_METHODS, gene_selection=GENE_SELECTIONS),
        cms = expand(rules.cms.output, sample=SAMPLES, integration_method=INTEGRATION_METHODS, gene_selection=GENE_SELECTIONS),
        ldfdiff = expand(rules.ldfdiff.output, sample=SAMPLES, integration_method=INTEGRATION_METHODS, gene_selection=GENE_SELECTIONS)
    output:
        combined = join(out_metrics, config["name"] + "_combined_metrics.tsv")
    shell:
        """
        python3 pipeline_scripts/combine_metrics.py \
            --scib_metrics {input.scib_metrics} \
            --bnmi {input.bnmi} \
            --cms {input.cms} \
            --ldfdiff {input.ldfdiff} \
            --outfile {output.combined}
        """
