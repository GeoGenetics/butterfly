import pandas as pd
import glob

"""
supersimple automation using snakemake
- takes SMDB and a list of finished libs (these are automatically generated/downloaded every day)
- runs the PMP (post mapping pipeline)
- starts the R script that makes fancy plots

run like this:
    snakemake -s scripts/pmp.smk -j10 -pn --config smdb=db/2024-11-05.query_result_full.csv finished_libs=files/testfiles.csv

problems: 
- we need to specify ds or ss in metadmg (currently hardcoded ds) but the samples could be a mixture, might need a new arg 
- bamfilter cant run on these big ass files, need to reimplement
- need to change base_outdir to somwhere sensible 
- need to remove FakeCountry and FakeCore from data[], this is just for testing 
- a lot of hardcoded vars in metadmg/bamfilter should be moved up 
- butteryfly R script needs to be updated to take args
"""

# ----------------- global vars ----------------- #

nodes = "/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp"
names = "/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp"
acc2tax = "/projects/caeg/data/db/mikkels/combined_accession2taxid_20221112.gz"
base_outdir = "tmp/" # to change to something like "/projects/caeg/data/pmp/"


# ----------------- read in metadata ----------------- #

smdb_file = config["smdb"]
finished_file = config["finished_libs"]

# get SMDB and the list of finished libraries as pandas df 
magnusdb = pd.read_csv(smdb_file, low_memory=False)
finished = pd.read_csv(finished_file, header=None, names=["basedir", "basename"], skipinitialspace=True)
finished["archive_id"] = finished["basename"].str.split("_").str[0]
finished["library_id"] = finished["basename"].str.split("_").str[1]

# merge 
data = pd.merge(finished, magnusdb, on="library_id", how="left") 

# find the path of the final output bam 
def find_bam_path(row):
    basedir = row["basedir"]
    bam_pattern = f"{basedir}/results/align/*/*.bam" # this is needed as the format is not fixed (:)
    bam_paths = glob.glob(bam_pattern)
    return bam_paths[0] if len(bam_paths)==1 else None # it should always be 1, otherwise terrible things 

data["bam_path"] = data.apply(find_bam_path, axis=1)

# just keep cols we need, and no dups 
data = data[["library_id", "bam_path", "archive_id", "country_ocean", "basename", "field_sample_parent_id"]]
data = data.drop_duplicates(keep="first")

# add some data here for testing 
data["country_ocean"] = "FakeCountry" # this is needed for testing
data["field_sample_parent_id"] = "FakeCore"

# outdir depends on the country and master core id 
data["outdir"] = data.apply(lambda row: f"{base_outdir}{row['country_ocean']}/{row['field_sample_parent_id']}", axis=1)


# ----------------- functions for input ----------------- #

def get_bam(wildcards):
    bam_file = data.loc[data['basename'] == wildcards.filebase, 'bam_path'].values[0]
    return bam_file

archive_group = data.groupby(['archive_id', 'outdir'])['basename'].apply(list).to_dict()
def get_filtered_bams(wildcards):
    filebases = archive_group[(wildcards.archive_id, wildcards.outdir)]
    return [f"{wildcards.outdir}/library/bamfilter/{filebase}.filtered.bam" for filebase in filebases]


# ----------------- snake ----------------- #

rule all:
    input:
        expand("{outdir}/stats/metadmg/aggregate/{archive_id}.stat.gz", zip, outdir=data["outdir"], archive_id=data["archive_id"]),
        expand("{outdir}/bamfilter/{archive_id}.stats.tsv", zip, outdir=data["outdir"], archive_id=data["archive_id"]),
        expand("{outdir}/bamfilter/{archive_id}.stats_filtered.tsv", zip, outdir=data["outdir"], archive_id=data["archive_id"])

rule reassign:
    input:
        get_bam
    output:
        "{outdir}/library/bamfilter/{filebase}.reassign.bam"
    params:
        tmp_dir = "tmp/",
    threads: 10
    shell:
        """
        filterBAM reassign --threads {threads} --bam {input} --tmp-dir {params.tmp_dir} --out-bam {output} --iters 0 --min-read-ani 94 --min-read-count 3
        """

rule filterbam:
    input:
        bam = "{outdir}/library/bamfilter/{filebase}.reassign.bam"
    output:
        bam_filtered = "{outdir}/library/bamfilter/{filebase}.filtered.bam",
        stats = "{outdir}/library/bamfilter/{filebase}.stats.tsv",
        stats_filtered = "{outdir}/library/bamfilter/{filebase}.stats_filtered.tsv"
    params:
        tmp_dir="tmp/"
    threads: 10
    shell:
        """
        filterBAM filter --threads {threads} --bam {input.bam} --tmp-dir {params.tmp_dir} --bam-filtered {output.bam_filtered} --stats {output.stats} --stats-filtered {output.stats_filtered} --min-read-ani 94 --min-read-count 3 --min-expected-breadth-ratio 0.5 --min-normalized-entropy auto --min-normalized-gini auto --min-breadth 0 --min-avg-read-ani 90 --min-coverage-evenness 0.4 --min-coverage-mean 0 --include-low-detection
        """

rule merge_bams:
    input:
        get_filtered_bams
    output:
        "{outdir}/merged_bams/{archive_id}.bam"
    shell:
        """
        samtools merge {output} {input}
        """

rule sort_bam:
    input:
        merged_bam="{outdir}/merged_bams/{archive_id}.bam"
    output:
        sorted_bam="{outdir}/merged_bams/{archive_id}.sorted.bam"
    threads: 8
    shell:
        "samtools sort {input.merged_bam} -n -m10G -@{threads} -o {output.sorted_bam}"


rule metadmg_lca:
    input:
        sorted_bam="{outdir}/merged_bams/{archive_id}.sorted.bam"
    output:
        bdamage="{outdir}/results/metadmg/lca/{archive_id}.bdamage.gz",
        lca_stat="{outdir}/results/metadmg/lca/{archive_id}.stat.gz",
        lca="{outdir}/results/metadmg/lca/{archive_id}.lca.gz"
    threads: 32
    params:
        "{outdir}/results/metadmg/lca/{archive_id}"
    shell:
        """
        /projects/caeg/apps/metaDMG-cpp/metaDMG-cpp lca --threads {threads} \
        --bam {input.sorted_bam} \
        --nodes {nodes} \
        --names {names} \
        --acc2tax {acc2tax} \
        --weight_type 1 --fix_ncbi 0 --sim_score_low 0.95 --how_many 30 --temp /projects/caeg/people/dlm551/tmp \
        --out_prefix {params}
        """

rule metadmg_getdmg:
    input:
        sorted_bam="{outdir}/merged_bams/{archive_id}.sorted.bam"
    output:
        "{outdir}/results/metadmg/damage/{archive_id}.bdamage.gz"
    params:
        "{outdir}/results/metadmg/damage/{archive_id}"
    threads: 32
    shell:
        """
         /projects/caeg/apps/metaDMG-cpp/metaDMG-cpp getdamage --threads {threads} --run_mode 0 --min_length 30 --print_length 30 \
        --out_prefix {params} {input.sorted_bam}
        """

rule metadmg_dfit:
    input:
        bdamage="{outdir}/results/metadmg/lca/{archive_id}.bdamage.gz"
    output:
        dfit="{outdir}/results/metadmg/dfit/{archive_id}.dfit.gz"
    params:
        "{outdir}/results/metadmg/dfit/{archive_id}"
    threads: 32
    shell:
        """
        /projects/caeg/apps/metaDMG-cpp/metaDMG-cpp dfit {input.bdamage} --threads {threads} \
        --nodes {nodes} \
        --names {names} \
        --lib ds --nopt 10 --doboot 1 --nbootstrap 20 --showfits 2 --seed 31924 \
        --out_prefix {params}
        """

rule metadmg_aggregate:
    input:
        bdamage="{outdir}/results/metadmg/lca/{archive_id}.bdamage.gz",
        lca_stat="{outdir}/results/metadmg/lca/{archive_id}.stat.gz",
        dfit="{outdir}/results/metadmg/dfit/{archive_id}.dfit.gz"
    output:
        aggregate_stat="{outdir}/stats/metadmg/aggregate/{archive_id}.stat.gz"
    params:
        "{outdir}/stats/metadmg/aggregate/{archive_id}"
    shell:
        """
        /projects/caeg/apps/metaDMG-cpp/metaDMG-cpp aggregate {input.bdamage} \
        --nodes {nodes} \
        --names {names} \
        --lcastat {input.lca_stat} --dfit {input.dfit} \
        --out_prefix {params}
        """

rule filterbam_getstats:
    input:
        sorted_bam="{outdir}/merged_bams/{archive_id}.sorted.bam"
    output:
        stats = "{outdir}/bamfilter/{archive_id}.stats.tsv",
        stats_filtered = "{outdir}/bamfilter/{archive_id}.stats_filtered.tsv"
    params:
        tmp_dir="tmp/"
    threads: 10
    shell:
        """
        filterBAM filter --threads {threads} --bam {input.sorted_bam} --tmp-dir {params.tmp_dir} --stats {output.stats} --stats-filtered {output.stats_filtered}
        """