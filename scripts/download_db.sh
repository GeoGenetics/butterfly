
check_file () {
	if [ ! -f $1 ]
	then 
		echo "cant find"
		exit 
	fi
}

today=$(/usr/bin/date +%y%m%d)
~/.conda/envs/caeg/bin/python /maps/projects/caeg/scratch/dlm551/robotbutterfly/scripts/download.py > /home/dlm551/thing.txt 2>&1 
echo "Database downloaded"
db=$(find /maps/projects/caeg/scratch/dlm551/robotbutterfly/downloads/ -name "*tsv.gz" -mmin -20)
check_file $db



zcat $db > /maps/projects/caeg/scratch/dlm551/robotbutterfly/db/$today.megatable.tsv 
rm /maps/projects/caeg/scratch/dlm551/robotbutterfly/downloads/* 
echo "Database saved in db/$today.megatable.tsv"