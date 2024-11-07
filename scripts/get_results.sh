#!/usr/bin/bash

today=$(/usr/bin/date +%y%m%d)

fcdirs=/maps/projects/caeg/scratch/dlm551/robotbutterfly/files/$today.fcdirs.txt 
files=/maps/projects/caeg/scratch/dlm551/robotbutterfly/files/$today.files.csv 

>$fcdirs 
>$files


for f in $(/usr/bin/find /projects/caeg/data/production -maxdepth 1 -mindepth 1 -type d|grep -v others|grep -v _)
do
  dir=$(echo $f|rev|cut -f1 -d/|rev)
  c=$(echo -n $dir|wc -c)
  if [ $c -eq 4 ] 
  then 
    fcdir=$(find $f -maxdepth 3 -mindepth 3 -type d)
  elif [ $c -eq 3 ]
  then 
    fcdir=$(find $f -maxdepth 4 -mindepth 4 -type d)
  fi 
  echo "${fcdir}" >> $fcdirs
done 

while read fc 
do
  agg=$(find $fc -type f -path "*/stats/metadmg/aggregate/*" -name "*_collapsed.stat.gz" -mmin +30 2>/dev/null)
  if [ ! -z "${agg}" ]
  then 
    basefile=$(basename $agg _collapsed.stat.gz)
    echo $fc, $basefile >> $files 
  fi 
done < $fcdirs


#/home/dlm551/.conda/envs/why/bin/snakemake -s /projects/caeg/apps/automated_plots/helpers/launch.smk --config filelist=/maps/projects/caeg/apps/automated_plots/files/$today.csv -p -j20 > /home/dlm551/allout.txt 2>&1 



