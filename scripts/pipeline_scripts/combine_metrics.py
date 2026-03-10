import os
import pandas as pd
from functools import reduce
import csv
import argparse


#------------------------------------------------------------------------------#
# COMMAND LINE ARGUMENTS
#------------------------------------------------------------------------------#
parser = argparse.ArgumentParser()
parser.add_argument("--scib_metrics", nargs="+", required=True)
parser.add_argument("--bnmi", nargs="+", required=True)
parser.add_argument("--cms", nargs="+", required=True)
parser.add_argument("--ldfdiff", nargs="+", required=True)
parser.add_argument("--outfile", required=True)

args = parser.parse_args()
#------------------------------------------------------------------------------#
# MAIN FUNCTION
#------------------------------------------------------------------------------#
def combine_metrics(scib_metrics, bnmi, cms, ldfdiff, outfile):
    # HELPER FUNCTION TO GET DELIMITER, IF NOT TSV (write.table from R has ' ' by default)
    def get_delimiter(file_path, bytes = 4096):
        sniffer = csv.Sniffer()
        data = open(file_path, "r").read(bytes)
        delimiter = sniffer.sniff(data).delimiter
        return delimiter
    result_list = [scib_metrics, bnmi, cms, ldfdiff]
    metric_type = {
        "Bio conservation" : [
            "Isolated labels", # isolated label F1
            "Silhouette label", # isolated label ASW
            "bNMI", # balanced NMI
            "cLISI",
            "ldfDiff", # local density factor difference from CellMixS
            "Graph connectivity"
        ],
        "Batch correction" : [
            "PCR comparison",
            "iLISI",
            "cms", # cell-specific mixing score from CellMixS
            "BRAS"
        ]
    }
    final_df_list = []
    for files in result_list:
        tmp_df_list = []
        for file in files:
            delim = get_delimiter(file)
            metrics = pd.read_csv(file, sep = delim)
            # print(metrics)
            # Ensure there are no duplicated column names
            metrics = metrics.loc[:, ~metrics.columns.duplicated()]
            metrics = metrics.loc[:, ~metrics.columns.isin(["Unnamed: 0","cell_cycle_conservation","hvg_overlap", "trajectory"])] # old scib, not scib_metrics
            tmp_df_list.append(metrics)

        combined_df = pd.concat(tmp_df_list, axis = 0)
        final_df_list.append(combined_df)
    final_df = reduce(lambda left, right: pd.merge(left, right, on = ["sample", "integration_method", "gene_selection"], how = "outer"), final_df_list)

    # re-calculate Bio conservation and Batch correction scores including the new metrics
    # that are not implemented in scib metrics
    # 1. subset the dictionary to the intersection with the columns of the final_df
    # 2. calculate the scores => weighted mean of the mean of the bio and batch metrics
    # 3. replace the scores in the final_df: Bio conservation, Batch correction, Total
    metric_type["Batch correction"] = list(set(metric_type["Batch correction"]).intersection(set(final_df.columns)))
    metric_type["Bio conservation"] = list(set(metric_type["Bio conservation"]).intersection(set(final_df.columns)))
    for idx, row in final_df.iterrows():
        bio_conservation = row[metric_type["Bio conservation"]].mean()
        batch_correction = row[metric_type["Batch correction"]].mean()
        final_df.loc[idx, "Bio conservation"] = bio_conservation
        final_df.loc[idx, "Batch correction"] = batch_correction
        final_df.loc[idx, "Total"] = 0.6 * bio_conservation + 0.4 * batch_correction

    final_df.to_csv(outfile, sep = "\t", index = False)

if __name__ == "__main__":
    combine_metrics(
        scib_metrics=args.scib_metrics,
        bnmi=args.bnmi,
        cms=args.cms,
        ldfdiff=args.ldfdiff,
        outfile=args.outfile
    )

