import pandas as pd
import glob

sequenced = pd.read_csv("files/testfiles.csv", header=None, names=["basedir", "basename"], skipinitialspace=True)
sequenced["archive_id"] = sequenced["basename"].str.split("_").str[0]
sequenced["library_id"] = sequenced["basename"].str.split("_").str[1]

bam_dict = {}
for _, row in sequenced.iterrows():
    basedir = row["basedir"]
    archive_id = row["archive_id"]
    bam_pattern = f"{basedir}/results/align/*/*.bam"
    bam_files = glob.glob(bam_pattern)
    bam_dict.setdefault(archive_id, []).extend(bam_files)

all_archive_ids = list(bam_dict.keys())

rule all:
    input:
        expand("merged_bams/{archive_id}.bam", archive_id=all_archive_ids)

rule merge_bams:
    input:
        lambda wildcards: bam_dict[wildcards.archive_id]
    output:
        "merged_bams/{archive_id}.bam"
    shell:
        """
        samtools merge {output} {input}
        """