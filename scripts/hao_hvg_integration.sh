

# if tcells_hao_pbmc_hvg.tsv does not exist
# run
gene_selection="hvg"
samplename="cd4_tcells"
out_dir="/data/akalin/akollot/anglemania_benchmark/results/Hao_Satija_Cell_2021_PBMC/tcells/cd4/"

out_gene_selection="${out_dir}/${samplename}_${gene_selection}.tsv"
out_integrated_h5ad="${out_dir}/${samplename}_${gene_selection}_scvi.h5ad"
out_metrics="${out_dir}/${samplename}_${gene_selection}_metrics.tsv"

if [ ! -f "${out_gene_selection}" ]; then
    echo "running prepare_inputs.R with ${gene_selection}"
    Rscript prepare_inputs.R \
        --infile "${out_dir}/${samplename}.h5ad" \
        --outfile "${out_gene_selection}" \
        --batch_key "orig.ident" \
        --gene_selection "${gene_selection}" || exit 1
fi

# if tcells_hao_pbmc_angl_scvi.h5ad does not exist
# run
if [ ! -f "${out_integrated_h5ad}" ]; then
    echo "running scvi integration with ${gene_selection}"
    python3 integration.py \
        --infile="${out_dir}/${samplename}.h5ad" \
        --feature_subset="${out_gene_selection}" \
        --outfile="${out_integrated_h5ad}" \
        --integration_method="scvi" \
        --batch_key="orig.ident" \
        --label_key="celltype.l2" || exit 1

fi

if [ ! -f "${out_metrics}" ]; then
    echo "running metrics"
    python3 metrics.py \
        --original_h5ad="${out_dir}/${samplename}.h5ad" \
        --integrated_h5ad="${out_integrated_h5ad}" \
        --outfile="${out_metrics}" \
        --sample="${samplename}" \
        --gene_selection="${gene_selection}" \
        --integration_method="scvi" \
        --batch_key="orig.ident" \
        --label_key="celltype.l2" || exit 1
fi