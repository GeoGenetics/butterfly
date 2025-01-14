
module load seqtk

for base in $(cut -f1-2 -d, ../250113.cambodia.final.csv|sed 1d)
do
	basedir=$(echo $base|cut -f1 -d,)
	basename=$(echo $base|cut -f2 -d,)
	lca=$basedir/results/metadmg/lca/${basename}_collapsed.lca.gz
	echo "zcat $lca|grep '9605:\"Homo\":\"genus\"'|cut -f1 > readids/$basename.txt"
done |parallel -j20 



for base in $(cut -f1-2 -d, ../250113.cambodia.final.csv|sed 1d)
do
	basedir=$(echo $base|cut -f1 -d,)
	basename=$(echo $base|cut -f2 -d,)
	ids=readids/$basename.txt
	fastq=$basedir/results/reads/low_complexity/${basename}_collapsed.fastq.gz
	echo "seqtk subseq $fastq $ids > fastqs/$basename.fastq"
done |parallel -j20


