# See tutorial at: http://tiny.cc/snakemake_tutorial
import os
import pandas as pd
from os.path import join
from snakemake.utils import validate
# configfile: "/home/akollot/projects/CRC1588/dataset_integration/anglemania_benchmark/scripts/snakemake_config.yml"
#### COMMON RUN SPECIFICS ####
out_preprocessed = join(config["output_dir"], config["out_preprocessed"])
out_integration = join(config["output_dir"], config["out_integration"])
out_metrics = join(config["output_dir"], config["out_metrics"])


############## INPUT #################
SAMPLESHEET = pd.read_csv(config["samplesheet"], sep = "\t")
# SAMPLESHEET = SAMPLESHEET[SAMPLESHEET.ngroup <= 4] # 2 or 4 cell types per batch
# SAMPLESHEET = SAMPLESHEET[SAMPLESHEET.nbatch <= 5] # 2,3,4, or 5 batches
# SAMPLESHEET = SAMPLESHEET.iloc[1:10,:]
# validate(SAMPLESHEET, "samplesheet_schema.yaml") #FIXME: no schema yet. check samplesheet format. needs sample_name and original_h5ad columns
file_paths = SAMPLESHEET.set_index("sample_name")["file_path"].to_dict()

SAMPLES = file_paths.keys()
INTEGRATION_METHODS = config["integration_methods"]
GENE_SELECTIONS = config['gene_selection']


include: "rules/preprocess.smk"
include: "rules/integrate.smk"
include: "rules/metrics.smk"

rule all:
    input:
        expand(rules.preprocess.output, sample=SAMPLES, gene_selection=GENE_SELECTIONS),
        expand(rules.integrate.output.integrated, sample=SAMPLES, integration_method=INTEGRATION_METHODS, gene_selection=GENE_SELECTIONS),
        expand(rules.metrics.output, sample=SAMPLES, integration_method=INTEGRATION_METHODS, gene_selection=GENE_SELECTIONS),
        rules.combine.output