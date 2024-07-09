

while read bam 
do
	outname=merged/$(basename $bam .bam).sorted.bam
	if [ ! -f $outname ]; then 
		echo "samtools sort $bam -n -m10G -@4 -o $outname"
	fi
done < negs.txt


for bam in $(ls merged/Lib*.sorted.bam); do
	if [ ! -f new/results/metadmg/damage/$(basename $bam|cut -f1-3 -d.).stat ]; then 
		echo "/projects/caeg/apps/metaDMG-cpp/metaDMG-cpp lca --threads 8 \
		--bam $bam \
		--nodes /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp \
		--names /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp \
		--acc2tax /projects/caeg/data/db/mikkels/combined_accession2taxid_20221112.gz \
		--weight_type 1 --fix_ncbi 0 --sim_score_low 0.95 --how_many 30 --temp /projects/caeg/people/dlm551/tmp \
		--out_prefix new/results/metadmg/lca/$(basename $bam|cut -f1-3 -d.)"

		echo "/projects/caeg/apps/metaDMG-cpp/metaDMG-cpp getdamage --threads 8 --run_mode 0 --min_length 30 --print_length 30 \
		--out_prefix new/results/metadmg/damage/$(basename $bam|cut -f1-3 -d.) $bam"
	fi
done < negs.txt
 
for bdamage in $(ls new/results/metadmg/lca/Lib*bdamage.gz)
do
	if [ ! -f new/results/metadmg/dfit/$(basename $bdamage .bdamage.gz).dfit.gz ] ; then 
	echo "/projects/caeg/apps/metaDMG-cpp/metaDMG-cpp dfit $bdamage --threads 20 \
	--names /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp \
	--nodes /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp \
	--lib ds --nopt 10 --doboot 1 --nbootstrap 20 --showfits 2 --seed 31924 \
	--out_prefix new/results/metadmg/dfit/$(basename $bdamage .bdamage.gz) "
fi


	if [ ! -f new/stats/metadmg/aggregate/$(basename $bdamage .bdamage.gz).stat.gz ]; then #re pull this and make sure it gets damage est 
		echo "metaDMG-cpp/metaDMG-cpp aggregate $bdamage \
		--nodes /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp \
		--names /projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp \
		--lcastat new/results/metadmg/lca/$(basename $bdamage .bdamage.gz).stat.gz --dfit new/results/metadmg/dfit/$(basename $bdamage .bdamage.gz).dfit.gz \
		--out_prefix new/stats/metadmg/aggregate/$(basename $bdamage .bdamage.gz) "
	fi
done 

for f in $(ls merged/Lib*.sorted.bam); 
do 
 	echo filterBAM filter -t 10 -N --bam $f --stats results/bamfilter/$(basename $f).comp.reassign2.stats.tsv.gz --stats-filtered results/bamfilter/$(basename $f).comp.reassign2.stats-filtered.tsv.gz #--bam-filtered results/bamfilter/$(basename $f).comp.reassign2.filtered.bam
done 
