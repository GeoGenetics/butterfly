
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


for fastq in $(ls fastqs/*)
do
	echo "bwa aln -l 512 -t 12 /projects/lundbeck/data/hs37d5/hs37d5.fa $fastq > mapping/$(basename $fastq|cut -f1 -d.).sai"
done |parallel -j3

for sai in $(ls mapping/*sai)
do
	fastq=fastqs/$(basename $sai|cut -f1 -d.).fastq
	echo "bwa samse -r \"@RG\\tID:$(basename $sai|cut -f1 -d.)\\tSM:$(basename $sai|cut -f1 -d.)\\tCN:CGG\\tPL:ILLUMINA\\tLB:$(basename $sai|cut -f1 -d.)\\tDS:seqcenter@sund.ku.dk\" \
	 /projects/lundbeck/data/hs37d5/hs37d5.fa $sai $fastq > mapping/$(basename $sai|cut -f1 -d.).sam"
done |parallel -j20

for sam in $(ls mapping/*sam)
do
	echo "samtools sort $sam -m4G | samtools view -F4 -b - > bams/$(basename $sam .sam).bam"
done |parallel -j10

samtools merge bams/merged.bam $(ls bams/*bam|grep -v merged) 

for bam in $(ls bams/*bam)
do
	echo "bash /projects/lundbeck/scratch/lundbeck-pipeline/analyses/haplogrep.sh -b $bam -f /projects/lundbeck/data/hs37d5/hs37d5.fa -o tmp/$(basename $bam)"
done 
