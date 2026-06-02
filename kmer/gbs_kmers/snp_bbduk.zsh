#!/bin/bash

bbduk='/pub/soft/bbmap/bbduk.sh'

ref=$1
shift
dir=$1

output_file="$dir/${dir:t}_wgs_count.tsv.gz"
error_file="$dir/${dir:t}_wgs_count.err"

# 检查输出文件是否存在且大小大于0
if [ -s "$output_file" ]; then
    echo "Output file $output_file already exists and is not empty. Skipping..."
    exit 0
fi

{
    find $dir | egrep 'fq.gz$|fastq.gz$' | xargs gzip -cd \
    | $bbduk -Xmx10G -k=31 threads=4 ignorebadquality=t maskmiddle=f overwrite=t int=f ref=$ref in=stdin.fq out=/dev/null rpkm=/dev/stdout nzo=f \
    | grep -v '^#' | cut -f 1,5 | awk 'NR % 2 {printf $0"\t"; next} 1' | tr ':' '\t' \
    | cut -f 1,2,4,8 | gzip > "$output_file"
} 2> "$error_file"
