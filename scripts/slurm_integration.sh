#!/bin/bash
#SBATCH --output=/home/akollot/projects/CRC1588/dataset_integration/anglemania_benchmark/slurm_out/%j.out
#SBATCH --error=/home/akollot/projects/CRC1588/dataset_integration/anglemania_benchmark/slurm_out/%j.out
#SBATCH --gpus=1
#SBATCH --mem=124G                     
#SBATCH --ntasks=1 --cpus-per-task=10
#SBATCH --time=96:00:00

cd /home/akollot/projects/CRC1588/dataset_integration/anglemania_benchmark/scripts/

#---------------------------------------------------------------------------------#
# INPUTS
# ---------------------------------------------------------------------------------#
## Example:
dataset=Hao_Satija_Cell_2021_PBMC
hvg_anglemania=anglemania
integration_method=seurat

if [ "$dataset" == "Hao_Satija_Cell_2021_PBMC" ]; then
    batch_key="orig.ident"
    label_key="celltype.l3"
fi

if [ "$dataset" == "pbmcsca" ]; then
    batch_key="batch"
    label_key="CellType"
fi

original_h5ad=/data/akalin/akollot/anglemania_benchmark/data/${dataset}/adata_${hvg_anglemania}.h5ad
infile="/data/akalin/akollot/anglemania_benchmark/data/${dataset}/adata_${hvg_anglemania}_${integration_method}.h5ad"
cp ${original_h5ad} ${infile}

outdir=/data/akalin/akollot/anglemania_benchmark/results/${dataset}/${hvg_anglemania}/

if [[ "$integration_method" == "harmony" || "$integration_method" == "scvi" || "$integration_method" == "scanvi" ]]; then
    python3 integration.py \
        --infile="/data/akalin/akollot/anglemanuel_benchmark/data/${dataset}/adata_${hvg_anglemania}.h5ad" \
        --outdir="/data/akalin/akollot/anglemanuel_benchmark/results/${dataset}/${hvg_anglemania}/" \
        --integration_method="${integration_method}" \
        --batch_key=${batch_key} \
        --label_key=${label_key}
elif [[ "$integration_method" == "seurat" ]]; then
    Rscript seurat_integration.R --infile ${infile} \
                                 --outdir ${outdir} \
                                 --batch_key ${batch_key} \
                                 --label_key ${label_key}
fi

rm ${infile}
