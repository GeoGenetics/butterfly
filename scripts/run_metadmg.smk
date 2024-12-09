import os
import pandas as pd

# ----------------- Functions ----------------- #
def get_bam_list(bamlist_file):
    bamlist = pd.read_csv(bamlist_file, low_memory=False, header=None)
    return bamlist[0].tolist()  

bamlist_file = config['list']
bam_files = get_bam_list(bamlist_file)
base_outdir = config["base_outdir"]
IDENTITY_CUTOFF = 0.90
SORT_BAMS = False 
sample_to_bam = {os.path.basename(bam).replace('.bam', ''): bam for bam in bam_files}

nodes = config["nodes"]
names = config["names"]
acc2tax = config["acc2tax"]

# ----------------- Rules ----------------- #
rule all:
    input:
        expand(f"{base_outdir}/{{sample}}.agg.stat.gz", sample=[os.path.basename(bam).replace('.bam', '') for bam in bam_files])

if SORT_BAMS:
    rule sort_bam:
        input:
            bam=lambda wildcards: sample_to_bam[wildcards.sample]
        output:
            sorted_bam=f"{base_outdir}/{{sample}}.sorted.bam"
        params:
            sample=lambda wildcards: os.path.basename(wildcards.bam_file).replace('.bam', ''),
        threads: 4
        shell:
            "samtools sort {input.bam} -n -m10G -@{threads} -o {output.sorted_bam}"

rule metadmg_lca:
    input:
        bam=lambda wildcards: (
            f"{base_outdir}/{wildcards.sample}.sorted.bam"
            if SORT_BAMS
            else sample_to_bam[wildcards.sample]
        )
    output:
        bdamage=f"{base_outdir}/{{sample}}.bdamage.gz",
        lca_stat=f"{base_outdir}/{{sample}}.stat.gz",
        lca=f"{base_outdir}/{{sample}}.lca.gz"
    threads: 4
    params:
        f"{base_outdir}/{{sample}}"
    shell:
        """
        /projects/caeg/people/dlm551/metaDMG-cpp/metaDMG-cpp lca --threads {threads} \
        --bam {input.bam} \
        --nodes {nodes} \
        --names {names} \
        --acc2tax {acc2tax} \
        --weight_type 1 --fix_ncbi 0 --sim_score_low {IDENTITY_CUTOFF} --how_many 30 --temp /projects/caeg/people/dlm551/tmp \
        --out_prefix {params}
        """

rule metadmg_dfit:
    input:
        bdamage=f"{base_outdir}/{{sample}}.bdamage.gz"
    output:
        dfit=f"{base_outdir}/{{sample}}.dfit.gz"
    params:
        f"{base_outdir}/{{sample}}"
    threads: 4
    shell:
        """
        /projects/caeg/people/dlm551/metaDMG-cpp/metaDMG-cpp dfit {input.bdamage} --threads {threads} \
        --nodes {nodes} \
        --names {names} \
        --lib mix --nopt 10 --doboot 1 --nbootstrap 20 --showfits 2 --seed 31924 \
        --out_prefix {params}
        """

rule metadmg_aggregate:
    input:
        bdamage=f"{base_outdir}/{{sample}}.bdamage.gz",
        lca_stat=f"{base_outdir}/{{sample}}.stat.gz",
        dfit=f"{base_outdir}/{{sample}}.dfit.gz"
    output:
        aggregate_stat=f"{base_outdir}/{{sample}}.agg.stat.gz"
    params:
        f"{base_outdir}/{{sample}}.agg"
    shell:
        """
        /projects/caeg/people/dlm551/metaDMG-cpp/metaDMG-cpp aggregate {input.bdamage} \
        --nodes {nodes} \
        --names {names} \
        --lcastat {input.lca_stat} --dfit {input.dfit} \
        --out_prefix {params}
        """
