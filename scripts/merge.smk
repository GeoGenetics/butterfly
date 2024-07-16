import os

filelist = "pre-filt/list" #"/maps/projects/wintherpedersen/people/lnc113/COREX/DONE.txt"
outdir = "pre-filt/"


def get_cgg_bam_counts():
    bam_counts = {}
    with open(filelist) as f:
        lines = f.readlines()
        for line in lines:
            if "20231201" not in line:
                cgg = line.split("/")[-1].split("_")[0] #may have to constantly change this 
                if cgg not in bam_counts:
                    bam_counts[cgg] = 0
                bam_counts[cgg] += 1
    return bam_counts

cgg_bam_counts = get_cgg_bam_counts()
cggs = list(cgg_bam_counts.keys())

rule all:
    input:
        #expand("merged/{cgg}.{n}.DS.sorted.bam",zip, cgg=cgg_bam_counts.keys(), n=cgg_bam_counts.values())
        expand(outdir + "stats/metadmg/aggregate/{cgg}.{n}.DS.stat.gz", cgg=cgg_bam_counts.keys(), n=cgg_bam_counts.values()),
        #expand("new/results/bamfilter/{cgg}.{n}.comp.reassign2.filtered.bam", cgg=cgg_bam_counts.keys(), n=cgg_bam_counts.values())


rule create_bam_list:
    output:
        bamlist="tmp/{cgg}.txt"
    shell:
        """
        grep {wildcards.cgg} {filelist} | grep -v 20231201 > {output.bamlist}
        """

rule merge_bam:
    input:
        bamlist="tmp/{cgg}.txt"
    output:
        merged_bam=outdir + "merged/{cgg}.{n}.DS.bam"
    params:
        n=lambda wildcards: cgg_bam_counts[wildcards.cgg]
    threads: 8
    shell:
        "samtools merge {output.merged_bam} -b {input.bamlist} -@{threads}"

rule sort_bam:
    input:
        merged_bam=outdir + "merged/{cgg}.{n}.DS.bam"
    output:
        sorted_bam=outdir + "merged/{cgg}.{n}.DS.sorted.bam"
    threads: 8
    shell:
        "samtools sort {input.merged_bam} -n -m10G -@{threads} -o {output.sorted_bam}"

rule metadmg_lca:
    input:
        sorted_bam=outdir + "merged/{cgg}.{n}.DS.sorted.bam"
    output:
        bdamage=outdir + "results/metadmg/lca/{cgg}.{n}.DS.bdamage.gz",
        lca_stat=outdir + "results/metadmg/lca/{cgg}.{n}.DS.stat.gz",
        lca=outdir + "results/metadmg/lca/{cgg}.{n}.DS.lca.gz"
    threads: 32
    params:
        outdir + "results/metadmg/lca/{cgg}.{n}.DS"
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
        sorted_bam=outdir + "merged/{cgg}.{n}.DS.sorted.bam"
    output:
        outdir + "results/metadmg/damage/{cgg}.{n}.DS.bdamage.gz"
    params:
        outdir + "results/metadmg/damage/{cgg}.{n}.DS"
    threads: 32
    shell: 
        """
         /projects/caeg/apps/metaDMG-cpp/metaDMG-cpp getdamage --threads {threads} --run_mode 0 --min_length 30 --print_length 30 \
        --out_prefix {params} {input.sorted_bam}
        """

rule metadmg_dfit:
    input:
        bdamage=outdir + "results/metadmg/lca/{cgg}.{n}.DS.bdamage.gz"
    output:
        dfit=outdir + "results/metadmg/dfit/{cgg}.{n}.DS.dfit.gz"
    params:
        outdir + "results/metadmg/dfit/{cgg}.{n}.DS"
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
        bdamage=outdir + "results/metadmg/lca/{cgg}.{n}.DS.bdamage.gz",
        lca_stat=outdir + "results/metadmg/lca/{cgg}.{n}.DS.stat.gz",
        dfit=outdir + "results/metadmg/dfit/{cgg}.{n}.DS.dfit.gz"
    output:
        aggregate_stat=outdir + "stats/metadmg/aggregate/{cgg}.{n}.DS.stat.gz"
    params:
        outdir + "stats/metadmg/aggregate/{cgg}.{n}.DS"
    shell:
        """
        metaDMG-cpp/metaDMG-cpp aggregate {input.bdamage} \
        --nodes /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp \
        --names /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp \
        --lcastat {input.lca_stat} --dfit {input.dfit} \
        --out_prefix {params}
        """

rule filter_bam:
    input:
        sorted_bam=outdir + "merged/{cgg}.{n}.DS.sorted.bam"
    output:
        filtered_bam=outdir + "new/results/bamfilter/{cgg}.{n}.comp.reassign2.filtered.bam"
    threads: 4
    params:
        stats = outdir + "results/bamfilter/{cgg}.comp.reassign2.stats.tsv.gz",
        filtered = outdir + "results/bamfilter/{cgg}.comp.reassign2.stats-filtered.tsv.gz"
    shell:
        """
        filterBAM filter -m 16G -t {threads} -N --bam {input.sorted_bam} \
        --stats {params.stats} \
        --stats-filtered {params.filtered}  \
        --bam-filtered {output.filtered_bam}
        """
