import pandas as pd
import glob
import os

magnusdb = pd.read_csv("db/2024-11-05.query_result_full.csv", low_memory=False)
sequenced = pd.read_csv("files/241105.files.csv", header=None, names=["basedir", "basename"])
sequenced["archive_id"] = sequenced["basename"].str.split("_").str[0]
sequenced["library_id"] = sequenced["basename"].str.split("_").str[1]



merged_data = pd.merge(sequenced, magnusdb, on="library_id", how="left")


def find_bam_file(basedir):
    bam_path_pattern = os.path.join(basedir, "results", "align", "*", "*.bam")
    bam_files = glob.glob(bam_path_pattern)
    return bam_files[0] if bam_files else None

merged_data["bam_file"] = merged_data["basedir"].apply(find_bam_file)
print(merged_data.head())
print("Dimensions of sequenced DataFrame:", sequenced.shape)
print("Dimensions of magnusdb DataFrame:", magnusdb.shape)
print("Dimensions of merged_data DataFrame:", merged_data.shape)