awk '{print "samtools faidx ../Avena_atlantica.fasta.gz " $1 " > "$1".fasta &"}' chr.list > exseq.sh
