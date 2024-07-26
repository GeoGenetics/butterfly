output_file=$1
shift
files=("$@")

first=true
prog=cat 
for file in "${files[@]}"; do
    if $first_file; then
        echo "$file"; 
        if [ $(basename $file|rev|cut -f1 -d.|rev) == "gz" ]; then 
            prog=zcat 
        fi
        $prog "$file" | awk -v filename="$file" '{print filename "\t" $0}' > $output_file
        first_file=false
    else
        echo "$file"; 
        $prog "$file" | awk -v filename="$file" '{print filename "\t" $0}'|sed 1d >> $output_file
    fi
done
