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