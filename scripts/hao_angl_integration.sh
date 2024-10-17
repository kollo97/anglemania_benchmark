

# if tcells_hao_pbmc_hvg.tsv does not exist
# run
gene_selection="angl"
anglemania_mode="pearson"
samplename="cd4_tcells"
angl_threshold=2
out_dir="/data/akalin/akollot/anglemania_benchmark/results/Hao_Satija_Cell_2021_PBMC/tcells/cd4/"
out_gene_selection="${out_dir}/${samplename}_${gene_selection}_${anglemania_mode}_t${angl_threshold}.tsv"
out_integrated_h5ad="${out_dir}/${samplename}_${gene_selection}_${anglemania_mode}_t${angl_threshold}_scvi.h5ad"
out_metrics="${out_dir}/${samplename}_${gene_selection}_${anglemania_mode}_t${angl_threshold}_metrics.tsv"

if [ ! -f "${out_gene_selection}" ]; then
    echo "running prepare_inputs.R with ${gene_selection} and anglemania mode: ${anglemania_mode}"
    Rscript prepare_inputs.R \
        --infile "${out_dir}/${samplename}.h5ad" \
        --outfile "${out_gene_selection}" \
        --batch_key "orig.ident" \
        --gene_selection "${gene_selection}" \
        --anglemania_mode "${anglemania_mode}" \
        --anglemania_threshold "${angl_threshold}" || exit 1
fi

# if tcells_hao_pbmc_angl_scvi.h5ad does not exist
# run
if [ ! -f "${out_integrated_h5ad}" ]; then
    echo "running scvi integration with ${gene_selection} and anglemania mode: ${anglemania_mode}"
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
