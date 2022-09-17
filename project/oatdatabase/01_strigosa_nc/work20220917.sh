samtools faidx A.strigosa.fasta.gz Chr01 > Chr01.fasta &
samtools faidx A.strigosa.fasta.gz Chr02 > Chr02.fasta &
samtools faidx A.strigosa.fasta.gz Chr03 > Chr03.fasta &
samtools faidx A.strigosa.fasta.gz Chr04 > Chr04.fasta &
samtools faidx A.strigosa.fasta.gz Chr05 > Chr05.fasta &
samtools faidx A.strigosa.fasta.gz Chr06 > Chr06.fasta &
seqkit seq Chr03.fasta -u -o chr3A.fasta &
seqkit seq Chr04.fasta -u -o chr4A.fasta &
seqkit seq Chr05.fasta -u -o chr5A.fasta &
seqkit seq Chr06.fasta -u -o chr6A.fasta &
seqkit seq Chr07.fasta -u -o chr7A.fasta &
seqkit seq ./contig/chrUn.fasta -u -o chrUn.fasta
cat chr* > genome20220917/A_strigosa_nc_genome.fasta &
# gff rename chr
awk '{gsub(/^Chr01/,"chr1A");print $0}' ../A.strigosa.chr.HCC.gff | awk '{gsub(/^Chr02/,"chr2A");print $0}' | awk '{gsub(/^Chr03/,"chr3A");print $0}' | awk '{gsub(/^Chr04/,"chr4A");print $0}' | awk '{gsub(/^Chr05/,"chr5A");print $0}' | awk '{gsub(/^Chr06/,"chr6A");print $0}' | awk '{gsub(/^Chr07/,"chr7A");print $0}' | grep -v '^Chr' > A_strigosa_nc_sorted_chr.gff3 &
# 
grep '^chrUn' ../genome20220304/A.strigosa.gene.gff3 > A_strigosa_nc_sorted_chrUn.gff3
# merge gff file
cat A_strigosa_nc_sorted_chr.gff3 A_strigosa_nc_sorted_chrUn.gff3 > A_strigosa_nc_sorted.gff3
# try to rename transcript and cds and exon name
python transcript_rename.py
# change cds and protain name replace '*' 
awk '{if (/^>/) {gsub(/$/, .1); print} else if (/^.*$/){ print }}' A.strigosa.chr.HCC.cds > A_strigosa_nc_cds.fasta
awk '{if (/^>/) {gsub(/$/, .1); print} else if (/^.*$/){ print }}' ../A.strigosa.chr.HCC.pep | sed 's#\*##g' > A_strigosa_nc_pep.fasta
samtools faidx A_strigosa_nc_genome.fasta.gz
md5sum A_strigosa_nc_genome.fasta.gz > A_strigosa_nc_genome.fasta.gz.md5
bgzip A_strigosa_nc_cds.fasta
bgzip A_strigosa_nc_pep.fasta
bgzip A_strigosa_nc_sorted.gff3
