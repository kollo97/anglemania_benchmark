datasets=(Hao_Satija_Cell_2021_PBMC pbmcsca)
gene_selections=(hvg anglemania)
integration_methods=(harmony scvi scanvi seurat)
for dataset in "${datasets[@]}"; do
    for gene_selection in "${gene_selections[@]}"; do
        for integration_method in "${integration_methods[@]}"; do
            echo ${dataset} ${gene_selection} ${integration_method}
            sbatch --export=ALL,dataset=${dataset},integration_method=${integration_method},hvg_anglemania=${gene_selection} slurm_calc_metrics.sh
        done
    done
done

# python3 metrics.py \
#     --original_h5ad="/data/akalin/akollot/anglemania_benchmark/data/${dataset}/adata_${gene_selection}.h5ad" \
#     --integrated_h5ad="/data/akalin/akollot/anglemania_benchmark/results/${dataset}/${gene_selection}/${integration_method}.h5ad" \
#     --outdir="/data/akalin/akollot/anglemania_benchmark/results/${dataset}/${gene_selection}/" \
#     --integration_method="${integration_method}" \
#     --batch_key="orig.ident" \
#     --label_key="celltype.l3"

