import os
import pandas as pd


filelist = config["files"]
outdir = config["out"]

df = pd.read_csv(filelist)

def get_bam(wildcards):
    bam_file = df.loc[df['Filebase'] == wildcards.filebase, 'Basedir'].values[0]
    return f"{bam_file}/results/align/merge_alns/{wildcards.filebase}_collapsed.bam"


counts = df['ArchiveSampleID'].value_counts().to_dict()
archive_group = df.groupby('ArchiveSampleID')['Filebase'].apply(list).to_dict()

def get_filtered_bams(wildcards):
    filebases = archive_group[wildcards.archive_id]
    return [f"{outdir}/bamfilter/{filebase}_collapsed.filtered.bam" for filebase in filebases]



rule all:
    input:
        #expand(outdir + "/stats/metadmg/aggregate/{archive_id}.{n}.stat.gz",zip, archive_id=archive_group.keys(), n=[counts[aid] for aid in archive_group.keys()]),
        #expand(outdir + "/bamfilter/{archive_id}.{n}_collapsed.stats_filtered.tsv",zip, archive_id=archive_group.keys(), n=[counts[aid] for aid in archive_group.keys()]),
        outdir + "/bamfilter.tsv",
        outdir + "/metadamage.tsv"

wildcard_constraints:
    filebase = '|'.join(list(df['Filebase'])),


rule reassign:
    input:
        get_bam
    output:
        bam = outdir + "/bamfilter/{filebase}_collapsed.bam"
    params:
        tmp_dir = "tmp/",
        tmp_input = "tmp/{filebase}.bam"
    threads: 10
    shell:
        """
        if [ ! -f {params.tmp_input} ]; then
            ln -s {input} {params.tmp_input}
        fi
        filterBAM reassign --threads {threads} --bam {params.tmp_input} --tmp-dir {params.tmp_dir} --out-bam {output.bam} --iters 0 --min-read-ani 94 --min-read-count 3
        """


rule filterbam:
    input:
        bam = outdir + "/bamfilter/{filebase}_collapsed.bam"
    output:
        bam_filtered = outdir + "/bamfilter/{filebase}_collapsed.filtered.bam",
        stats = outdir + "/bamfilter/{filebase}_collapsed.stats.tsv",
        stats_filtered = outdir + "/bamfilter/{filebase}_collapsed.stats_filtered.tsv"
    params:
        tmp_dir="tmp/"
    threads: 10
    shell:
        """
        filterBAM filter --threads {threads} --bam {input.bam} --tmp-dir {params.tmp_dir} --bam-filtered {output.bam_filtered} --stats {output.stats} --stats-filtered {output.stats_filtered} --min-read-ani 94 --min-read-count 3 --min-expected-breadth-ratio 0.5 --min-normalized-entropy auto --min-normalized-gini auto --min-breadth 0 --min-avg-read-ani 90 --min-coverage-evenness 0.4 --min-coverage-mean 0 --include-low-detection
        """

rule merge_bam:
    input:
        get_filtered_bams
    output:
        merged_bam=outdir + "/merged/{archive_id}.{n}.bam"
    params:
        n=lambda wildcards: counts[wildcards.archive_id]
    threads: 8
    shell:
        "samtools merge {output.merged_bam} {input} -@{threads}"


rule sort_bam:
    input:
        merged_bam=outdir + "/merged/{archive_id}.{n}.bam"
    output:
        sorted_bam=outdir + "/merged/{archive_id}.{n}.sorted.bam"
    threads: 8
    shell:
        "samtools sort {input.merged_bam} -n -m10G -@{threads} -o {output.sorted_bam}"


rule metadmg_lca:
    input:
        sorted_bam=outdir + "/merged/{archive_id}.{n}.sorted.bam"
    output:
        bdamage=outdir + "/results/metadmg/lca/{archive_id}.{n}.bdamage.gz",
        lca_stat=outdir + "/results/metadmg/lca/{archive_id}.{n}.stat.gz",
        lca=outdir + "/results/metadmg/lca/{archive_id}.{n}.lca.gz"
    threads: 32
    params:
        outdir + "/results/metadmg/lca/{archive_id}.{n}"
    shell:
        """
        /projects/caeg/apps/metaDMG-cpp/metaDMG-cpp lca --threads {threads} \
        --bam {input.sorted_bam} \
        --nodes /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp \
        --names /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp \
        --acc2tax /projects/caeg/data/db/mikkels/combined_accession2taxid_20221112.gz \
        --weight_type 1 --fix_ncbi 0 --sim_score_low 0.95 --how_many 30 --temp /projects/caeg/people/dlm551/tmp \
        --out_prefix {params}
        """

rule metadmg_getdmg:
    input:
        sorted_bam=outdir + "/merged/{archive_id}.{n}.sorted.bam"
    output:
        outdir + "/results/metadmg/damage/{archive_id}.{n}.bdamage.gz"
    params:
        outdir + "/results/metadmg/damage/{archive_id}.{n}"
    threads: 32
    shell:
        """
         /projects/caeg/apps/metaDMG-cpp/metaDMG-cpp getdamage --threads {threads} --run_mode 0 --min_length 30 --print_length 30 \
        --out_prefix {params} {input.sorted_bam}
        """

rule metadmg_dfit:
    input:
        bdamage=outdir + "/results/metadmg/lca/{archive_id}.{n}.bdamage.gz"
    output:
        dfit=outdir + "/results/metadmg/dfit/{archive_id}.{n}.dfit.gz"
    params:
        outdir + "/results/metadmg/dfit/{archive_id}.{n}"
    threads: 32
    shell:
        """
        /projects/caeg/apps/metaDMG-cpp/metaDMG-cpp dfit {input.bdamage} --threads {threads} \
        --names /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp \
        --nodes /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp \
        --lib ds --nopt 10 --doboot 1 --nbootstrap 20 --showfits 2 --seed 31924 \
        --out_prefix {params}
        """

rule metadmg_aggregate:
    input:
        bdamage=outdir + "/results/metadmg/lca/{archive_id}.{n}.bdamage.gz",
        lca_stat=outdir + "/results/metadmg/lca/{archive_id}.{n}.stat.gz",
        dfit=outdir + "/results/metadmg/dfit/{archive_id}.{n}.dfit.gz"
    output:
        aggregate_stat=outdir + "/stats/metadmg/aggregate/{archive_id}.{n}.stat.gz"
    params:
        outdir + "/stats/metadmg/aggregate/{archive_id}.{n}"
    shell:
        """
        /projects/caeg/apps/metaDMG-cpp/metaDMG-cpp aggregate {input.bdamage} \
        --nodes /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp \
        --names /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp \
        --lcastat {input.lca_stat} --dfit {input.dfit} \
        --out_prefix {params}
        """

rule filterbam_getstats:
    input:
        sorted_bam=outdir + "/merged/{archive_id}.{n}.sorted.bam"
    output:
        stats = outdir + "/bamfilter/{archive_id}.{n}_collapsed.stats.tsv",
        stats_filtered = outdir + "/bamfilter/{archive_id}.{n}_collapsed.stats_filtered.tsv"
    params:
        tmp_dir="tmp/"
    threads: 10
    shell:
        """
        filterBAM filter --threads {threads} --bam {input.sorted_bam} --tmp-dir {params.tmp_dir} --stats {output.stats} --stats-filtered {output.stats_filtered}
        """

rule collate_bf:
    input:
        lambda wildcards: [ outdir + "/bamfilter/" + aid + "." + str(n) + "_collapsed.stats.tsv" for aid, n in zip(archive_group.keys(),[counts[aid] for aid in archive_group.keys()])],
    output:
         outdir + "/bamfilter.tsv"
    shell:
        "bash butterfly/scripts/collate_bf.sh  {output} {input}"


rule collate_md:
    input:
        lambda wildcards: [outdir + "/stats/metadmg/aggregate/" + aid + "." + str(n) + ".stat.gz" for aid, n in zip(archive_group.keys(),[counts[aid] for aid in archive_group.keys()])],
    output:
         outdir + "/metadamage.tsv"
    shell:
        "bash butterfly/scripts/collate_bf.sh  {output} {input}"
