
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
db=$(find /maps/projects/caeg/scratch/dlm551/robotbutterfly/downloads/ -name "*zip" -mmin -20)
check_file $db

unzip -j $db "*csv" -d /maps/projects/caeg/scratch/dlm551/robotbutterfly/tmp/
csv=$(find /maps/projects/caeg/scratch/dlm551/robotbutterfly/tmp/ -name "*csv" -mmin -20)
check_file $csv 
echo "Database unzipped"

mv $csv /maps/projects/caeg/scratch/dlm551/robotbutterfly/db/$today.megatable.csv 
rm /maps/projects/caeg/scratch/dlm551/robotbutterfly/downloads/* 
echo "Database saved in db/$today.megatable.csv"