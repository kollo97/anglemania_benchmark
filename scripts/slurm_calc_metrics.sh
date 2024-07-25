#!/bin/bash
#SBATCH --output=/home/akollot/projects/CRC1588/dataset_integration/anglemania_benchmark/slurm_out/%j.out
#SBATCH --error=/home/akollot/projects/CRC1588/dataset_integration/anglemania_benchmark/slurm_out/%j.out
#SBATCH --mem=64G                     
#SBATCH --ntasks=1 --cpus-per-task=10
#SBATCH --time=96:00:00

cd /home/akollot/projects/CRC1588/dataset_integration/anglemania_benchmark/scripts/
#---------------------------------------------------------------------------------#
# INPUTS
# ---------------------------------------------------------------------------------#
## Example:
# dataset=pbmcsca
# gene_selection=anglemania
# integration_method=harmony

if [ "$dataset" == "Hao_Satija_Cell_2021_PBMC" ]; then
    batch_key="orig.ident"
    label_key="celltype.l3"
fi

if [ "$dataset" == "pbmcsca" ]; then
    batch_key="batch"
    label_key="CellType"
fi

# ---------------------------------------------------------------------------------#
# RUN
# ---------------------------------------------------------------------------------#
original_h5ad="/data/akalin/akollot/anglemania_benchmark/data/${dataset}/adata_${hvg_anglemania}.h5ad"
tmp_h5ad="/data/akalin/akollot/anglemania_benchmark/data/${dataset}/adata_${hvg_anglemania}_${integration_method}.h5ad"
cp ${original_h5ad} ${tmp_h5ad}
outfile="/data/akalin/akollot/anglemania_benchmark/results/${dataset}/${hvg_anglemania}/${integration_method}_metrics.tsv"
# if the outfile exists, stop the job
if [ -f "${outfile}" ]; then
    echo "File ${outfile} already exists. Skipping..."
    exit 0
fi
echo "Calculating metrics, final output: ${outfile}"
python3 metrics.py \
    --original_h5ad="${tmp_h5ad}" \
    --integrated_h5ad="/data/akalin/akollot/anglemania_benchmark/results/${dataset}/${hvg_anglemania}/${integration_method}.h5ad" \
    --outfile="${outfile}" \
    --integration_method="${integration_method}" \
    --batch_key="${batch_key}" \
    --label_key="${label_key}"

rm ${tmp_h5ad}