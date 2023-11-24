## Salmon Pipline for OatBioDB

### Introduction

salmon is a tool for quantifying the expression of transcripts using RNA-seq data.

### user manual

This pipline will generate salmon quantification results for all samples in the input file(By three genomes).

#### Input

genome_list.txt: a list of three genomes, sep is tab
example:

```
/public/home/acfaa2ssz7/genome/sfs/salmon/ACD_sativa_sfs_cds.fasta.index	/public/home/acfaa2ssz7/genome/sfs/ACD_sativa_sfs.gtf	sfs
/public/home/acfaa2ssz7/genome/sang/salmon/ACD_sativa_sang_cds.fasta.index	/public/home/acfaa2ssz7/genome/sang/ACD_sativa_sang.gtf	sang
/public/home/acfaa2ssz7/genome/OT3098v2/salmon/ACD_sativa_OT3098v2_cds.fasta.index	/public/home/acfaa2ssz7/genome/OT3098v2/ACD_sativa_OT3098v2.gtf	otv2
```

#### salmon_samples.txt

a list of samples, sep is tab
example:

```
SRR22937062_1.fastq.gz	SRR22937062_2.fastq.gz	by01_1
SRR22937061_1.fastq.gz	SRR22937061_2.fastq.gz	by01_2
SRR22937058_1.fastq.gz	SRR22937058_2.fastq.gz	by01_3
SRR22937057_1.fastq.gz	SRR22937057_2.fastq.gz	by02_1
SRR22937056_1.fastq.gz	SRR22937056_2.fastq.gz	by02_2
SRR22937055_1.fastq.gz	SRR22937055_2.fastq.gz	by02_3
SRR22937054_1.fastq.gz	SRR22937054_2.fastq.gz	yzy01_1
SRR22937053_1.fastq.gz	SRR22937053_2.fastq.gz	yzy01_2
SRR22937052_1.fastq.gz	SRR22937052_2.fastq.gz	yzy01_3
SRR22937051_1.fastq.gz	SRR22937051_2.fastq.gz	yzy02_1
SRR22937060_1.fastq.gz	SRR22937060_2.fastq.gz	yzy02_2
SRR22937059_1.fastq.gz	SRR22937059_2.fastq.gz	yzy02_3

```

#### Command

```bash
python /public/home/acfaa2ssz7/soft/RNAseqScript/RNASeqSalmon.py -g /public/home/acfaa2ssz7/soft/RNAseqScript/genome_list.txt -s salmon_samples.txt -t 32
```
