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
# check if samplesheet is tsv or json
file_ext = os.path.splitext(config["samplesheet"])[1].lower()
if file_ext == ".tsv":
    SAMPLESHEET = pd.read_csv(config["samplesheet"], sep = "\t")
elif file_ext == ".json":
    SAMPLESHEET = pd.read_json(config["samplesheet"])
elif file_ext == "csv":
    SAMPLESHEET = pd.read_csv(config["samplesheet"], sep = ",")
else:
    raise ValueError("samplesheet must be csv, tsv or json") 
    # FIXME: add config schema. it should check samplesheet format. needs sample_name and file_path columns 
    # validate(SAMPLESHEET, "samplesheet_schema.yaml") 
# SAMPLESHEET = SAMPLESHEET[SAMPLESHEET.ngroup <= 4] # 2 or 4 cell types per batch
# SAMPLESHEET = SAMPLESHEET[SAMPLESHEET.nbatch <= 5] # 2,3,4, or 5 batches
# SAMPLESHEET = SAMPLESHEET.iloc[1:2,:]
# SAMPLESHEET = SAMPLESHEET[SAMPLESHEET.sample_name.str.contains("sim2")]

file_paths = SAMPLESHEET.set_index("sample_name")["file_path"].to_dict()

SAMPLES = file_paths.keys()
INTEGRATION_METHODS = config["integration_methods"]
GENE_SELECTIONS = config['gene_selection']

include: "rules/preprocess.smk"
include: "rules/integrate.smk"
include: "rules/metrics.smk"

rule all:
    input:
        # preprocess
        expand(rules.preprocess.output, 
               sample=SAMPLES, 
               gene_selection=GENE_SELECTIONS),

        # integrate
        expand(rules.integrate.output.integrated, 
               sample=SAMPLES, 
               integration_method=INTEGRATION_METHODS, 
               gene_selection=GENE_SELECTIONS),

        # metrics
        expand(rules.metrics.output, 
               sample=SAMPLES, 
               integration_method=INTEGRATION_METHODS, 
               gene_selection=GENE_SELECTIONS),

        # combine
        rules.combine.output