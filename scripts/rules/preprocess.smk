rule preprocess:
    input:
        # file_path="/local/Projects/AAkalin_Neuroblastoma/Results/simulate_single_cell/h5ad/{sample}.h5ad"
        original_h5ad=lambda wildcards: file_paths[wildcards.sample]

    output:
        outfile = join(out_preprocessed, "{sample}_{gene_selection}.tsv")
        
    params:
        batch_key = config["batch_key"],
        anglemania_mode = config["anglemania_mode"]
    # log:
    #     "/home/akollot/projects/CRC1588/dataset_integration/anglemania_benchmark/output/snakemake_logs/preprocess/{sample}_{gene_selection}_log_test.txt"
    shell:
        """
        Rscript prepare_inputs.R \
            --infile {input.original_h5ad} \
            --outfile {output.outfile} \
            --batch_key {params.batch_key} \
            --gene_selection {wildcards.gene_selection} \
            --anglemania_mode {params.anglemania_mode} 
        """