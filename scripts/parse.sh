> Tjornen.production_data2.tsv
for file in  stats/metadmg/aggregate/*.stat.gz;  do 
	echo "$file"; zcat "$file" | awk -v filename="$file" '{print filename "\t" $0}'|sed 1d >> Tjornen.production_data2.tsv ; 
done
printf "%s\t%s\n" "$(echo "file")" "$(zcat stats/metadmg/aggregate/*.stat.gz| head -1)" > Tjornen.metadmg_data.tsv
cat Tjornen.production_data2.tsv >> Tjornen.metadmg_data.tsv

>Tjornen.bam-filter_pre_stats2.tsv 
for file in  results/bamfilter/*stats.tsv.gz ;  do  	
	echo "$file"; zcat "$file" | awk -v filename="$file" '{print filename "\t" $0}'|sed 1d >>Tjornen.bam-filter_pre_stats2.tsv ;  
done
printf "%s\t%s\n" "$(echo "file")" "$(zcat results/bamfilter/*stats.tsv.gz | head -1 )" > Tjornen.bam-filter_stats.tsv
cat Tjornen.bam-filter_pre_stats2.tsv >> Tjornen.bam-filter_stats.tsv


# negatives 





# repliciates 
first=1

bffile=reps.Tjornen.bam-filter_stats.tsv
mdfile=reps.Tjornen.metadmg_data.tsv

while read lib 
do
	base=$(dirname $lib|cut -f1-10 -d/)
	agg=$base/stats/metadmg/aggregate/*.stat.gz
	bf=$base/results/bamfilter/*stats.tsv.gz
	if [ $first -eq 1 ]; then 
		printf "%s\t%s\n" "$(echo "$bf")" "$(zcat $bf | head -1 )" > $bffile
		printf "%s\t%s\n" "$(echo "$agg")" "$(zcat $agg | head -1 )" > $mdfile
		first=0 
	fi 
	echo "$bf"; zcat "$bf" | awk -v filename="filename" '{print filename "\t" $0}'|sed 1d >> $bffile; 
	echo "$agg"; zcat "$agg" | awk -v filename="filename" '{print filename "\t" $0}'|sed 1d >> $mdfile; 




done < library_list.txt 