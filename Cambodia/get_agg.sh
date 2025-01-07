nodes="/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp"
names="/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp"

tmp=$(mktemp)
cat $1 | sed 1d > $tmp 
while read line 
do
	master=$(echo $line|cut -f1 -d" ")
	basedir=$(echo $line|cut -f8 -d" ")
    
    bdamage=$(ls $basedir/results/metadmg/lca/*.bdamage.gz)
    lca_stat=$(ls $basedir/results/metadmg/lca/*.stat.gz 2>/dev/null)
    dfit=$(ls $basedir/results/metadmg/dfit/*.dfit.gz)

    aggregate_stat=output/$master/metadmg/$(basename $dfit .dfit.gz)
    mkdir -p $(dirname $aggregate_stat)

    if [ -z $lca_stat ]; then 
        agg=$(ls $basedir/stats/metadmg/aggregate/*gz)
        echo ln -s $agg $aggregate_stat.stat.gz
    else
        echo /projects/caeg/people/dlm551/metaDMG-cpp/metaDMG-cpp aggregate ${bdamage} \
        --nodes ${nodes} \
        --names ${names} \
        --lcastat ${lca_stat} --dfit ${dfit} \
        --out_prefix $aggregate_stat
    fi
done < $tmp 
rm $tmp 