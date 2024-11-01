filelist=/maps/projects/wintherpedersen/people/lnc113/COREX/DONE.txt


# merge 
for cgg in $(cat $filelist|cut -d/ -f9|sort|uniq)
do
	domerge=false
	bamlist=tmp/$cgg
	grep $cgg $filelist|grep -v 20231201 > $bamlist 
	n=$(wc -l $bamlist|cut -f1 -d" ")

	# Merge if merge does not exist 
	outname=merged/$cgg.$n.DS.bam 
	if [ ! -f "$outname" ]; then domerge=true; fi 

	if [ "$domerge" = true ];
	then 
		echo "samtools merge $outname -b $bamlist -@4"
	fi
done 


# sort 
for f in $(ls merged/*DS.bam)
do
	outname=merged/$(basename $f .bam).sorted.bam
	if [ ! -f $outname ]; then 
		echo "samtools sort $f -n -m10G -@4 -o $outname"
	fi
done 



# metadmg 
for bam in $(ls merged/*DS.sorted.bam)
do
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
done |parallel 
 

for bdamage in $(ls new/results/metadmg/lca/*bdamage.gz)
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
done |parallel 

#bamfilter filter withuot filtering 
for f in $(ls merged/*DS.sorted.bam); 
do 
 	echo filterBAM filter -m 16G -t 5 -N --bam $f --stats results/bamfilter/$(basename $f).comp.reassign2.stats.tsv.gz --stats-filtered results/bamfilter/$(basename $f).comp.reassign2.stats-filtered.tsv.gz --bam-filtered results/bamfilter/$(basename $f).comp.reassign2.filtered.bam
done > bamfiltercmds

while read x
do
	 sbatch -W -t 24:00:00 -c 5 --mem=200GB -p compproduction  -o tmp/slurm-%j.out --wrap "$x" &
done < bamfiltercmds

