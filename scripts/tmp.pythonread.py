import pandas as pd
import glob
import os

magnusdb = pd.read_csv("db/2024-11-05.query_result_full.csv", low_memory=False)
sequenced = pd.read_csv("files/testfiles.csv", header=None, names=["basedir", "basename"])
sequenced["archive_id"] = sequenced["basename"].str.split("_").str[0]
sequenced["library_id"] = sequenced["basename"].str.split("_").str[1]



merged_data = pd.merge(sequenced, magnusdb, on="library_id", how="left")

outdir = "tmp/"

def find_bam_path(row):
    basedir = row["basedir"]
    bam_pattern = f"{basedir}/results/align/*/*.bam"
    bam_paths = glob.glob(bam_pattern)
    return bam_paths[0] if len(bam_paths)==1 else None

merged_data["bam_path"] = merged_data.apply(find_bam_path, axis=1)


print(merged_data["bam_path"])
