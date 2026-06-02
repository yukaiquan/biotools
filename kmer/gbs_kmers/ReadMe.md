# k-mers for evaluating sequencing datasets

**⚠️ Note: This repository is cloned from the original source: https://bitbucket.org/ipk_dg_public/gbs_kmers/src/main/**

---

## Overview

This tool is designed to:
- Identify the level of heterozygosity
- Detect sample mix-ups
- Compare samples across sequencing platforms

**Core functionality**: Creation of unique SNP-informed k-mers (31-mer)

**Prerequisites**: Generate a high-quality SNP matrix from a raw VCF file (filtration criteria: biallelic, missingness < 10%, heterozygosity < 10%)

---

## Workflow

### 1. Create a bed interval file from filtered SNPs

The 31-mer interval is created with 16 and 15 nucleotides before and after SNP positions respectively.

```bash
# Define input file
vcf='240131_GBS_biallelic_sel_samples.vcf.gz'  # filtered VCF file path

# Generate bed file
zcat $vcf | grep -v '^#' | cut -f1-5 | \
  awk '{print $1"\t"($2-16)"\t"($2+15)"\t"$1":"$2":"$4":"$5}' > GBS_SNPs.bed
```

---

### 2. Generate k-mer (31-mer) sequence for reference and alternate alleles

#### 2.1 Set required program paths

```bash
# Required programs
bedtools='/opt/Bio/bedtools/2.30.0/bin/bedtools'
samtools='/opt/Bio/samtools/1.16.1/bin/samtools'
bwa='/opt/Bio/bwa/0.7.17/bin/bwa'

# Input files
bed='GBS_SNPs.bed'
ref='230103_Vfaba_parts_hedin_v2+Unanchored_contigs.fasta'
```

#### 2.2 Generate k-mer sequences

```bash
$bedtools getfasta -name -bed $bed -fi $ref | \
  cut -d : -f -4 | tr : '\t' | \
  awk 'NR % 2 == 1 {printf $0"\t"; next} 1' | \
  awk '{print $1":"$2":"$3"\n"$5"\n"$1":"$2":"$4"\n"substr($5, 1, 15) $4 substr($5, 17)}' > ${bed:r}.fasta

# Create index
$samtools faidx ${bed:r}.fasta
```

---

### 3. Identify unique SNP-kmers after mapping to reference genome

```bash
# Define input files
ref='230103_Vfaba_parts_hedin_v2+Unanchored_contigs.fasta'
qry='GBS_SNPs.fasta'

# BWA alignment
$bwa aln $ref $qry > ${qry:r}.sai
$bwa samse $ref ${qry:r}.sai $qry | \
  samtools view -SF 2324 -q 20 | \
  cut -f 1,3,4 | \
  tr ':' '\t' | \
  awk '$2 - $5 == 15' | \
  cut -f -2 | \
  tr '\t' ':' | \
  uniq -c | \
  awk '$1 == 2 {print $2}' | \
  grep -Ff - $qry.fai | \
  cut -f 1 | \
  xargs $samtools faidx $qry > ${qry:r}_uniq.fasta
```

---

### 4. Count SNP-informed k-mers in raw sequencing read files

Supports GBS, WGS, HiFi, or other sequencing platforms.

```bash
# Define reference k-mer file
ref='GBS_SNPs_uniq.fasta'

# Generate sample list (one sample directory per line)
find /filer/agruppen/jayakodi/core50_GBS -type d > sample_list.txt

# Run with multiple CPUs
cat sample_list.txt | parallel --will-cite -j 10 ./snp_bbduk.zsh $ref '{}'
```

**Example sample directory format**:
```
/filer-dg/agruppen/dg4/jayakodi/SHAPE_II_WGS_mapping_2020/core50_GBS/HOR_10886
```

---

### 5. Identify homozygosity & heterozygosity rate

**Environment**: R (version 3.5.1)

#### 5.1 Define R function

```r
# Function for calculating homo- and heterozygous rate
read_kmer_count <- function(info, mincov=1, cores=1) {
  f <- info
  rbindlist(mclapply(mc.cores=cores, 1:nrow(f), function(i) {
    fread(f[i, path], col.names=c("chromosome", "position", "dr", "dv"), head=F) -> z
    z[, c("dp", "dr") := list(dr+dv, NULL)]
    z[, variant.id := paste0(chromosome, ":", position)]
    z[dp >= mincov, r := dv/dp]
    z[r <= 0.1, gt := 0] 
    z[r >= 0.9, gt := 2] 
    z[r >= 0.4 & r <= 0.6, gt := 1] 
    z[!is.na(gt)] -> z
    z[, accession := f[i, accession]]
  })) 
}
```

#### 5.2 Read k-mer count tables

```r
# Read k-mer count tables (need to be adjusted based on user's folder path and name)
fread("find /filer-dg/agruppen/dg7/jayakodi/Vf_panAssemblies/ | grep '_gbs_count.tsv.gz$'", 
      head=F, col.names="path") -> f

# Extract sample names (need to be adjusted based on user's directory and file name)
f[, accession := sub("/ccs_reads/ccs_reads_gbs_count.tsv.gz$", "", 
                     sub("^.*_assembly/", "", path))]

f -> rds

# Import k-mer counts and call genotypes
read_kmer_count(info=rds, mincov=1, cores=1) -> gt
```

#### 5.3 Calculate homozygosity & heterozygosity rate

```r
# Calculate statistics
gt -> t3
table(t3$gt, t3$accession) -> c
as.data.table(c) -> c
dcast(c, V2 ~ V1, value.var='N') -> b

# Rename columns
setnames(b, c("accession", "ref", "het", "alt"))
b[, tot := (ref+het+alt)]

# Calculate rates
b[, pHet := (het/tot*100)]
b[, pHom := ((ref+alt)/tot*100)]

# Save results
fwrite(b, file="het_estimation_result.txt", sep="\t", quote=F, row.names=F)
```

---

### 6. Compare sample identity or mixups when sequenced with different platforms

**Prerequisite**: The `read_kmer_count` function should be defined

Create separate genotype tables for GBS and WGS datasets (any two sequencing datasets).

#### 6.1 GBS dataset processing

```r
# Read GBS data
fread("find /filer/agruppen/jayakodi/core50_GBS | grep '\\_gbs_count.tsv.gz$'", 
      head=F, col.names="path") -> gbs_f

# Extract sample names (needs to be adjusted according to user's directory and file name)
gbs_f[, accession := sub("_gbs_count.tsv.gz$", "", 
                         sub("^.*BRIDGE_HOR", "HOR", path))]

# Process data
read_kmer_count(info=gbs_f, mincov=1, cores=1) -> gbs_gt
gbs_gt[, .(variant.id, accession, gbsGT=gt)] -> gbs
```

#### 6.2 WGS dataset processing

```r
# Read WGS data
fread("find /filer/agruppen/jayakodi/WGS_panels/core50 | grep '\\_gbs_count.tsv.gz$'", 
      head=F, col.names="path") -> wgs_f

# Extract sample names (needs to be adjusted according to user's directory and file name)
wgs_f[, accession := sub("_gbs_count.tsv.gz$", "", 
                         sub("^.*core50:", "", path))]

# Process data
read_kmer_count(info=wgs_f, mincov=1, cores=1) -> wgs_gt
wgs_gt[, .(variant.id, accession, wgsGT=gt)] -> wgs
```

#### 6.3 Combine genotype tables and check correlation

```r
# Combine two genotype tables
wgs[gbs, on=c("variant.id", "accession"), nomatch=0] -> dt

# Check genotype correlation between GBS and WGS datasets
dt[, .(Cor = cor(wgsGT, gbsGT)), by=accession][order(Cor)] -> dtCor
```

#### 6.4 Calculate identity-by-state (IBS) statistics

```r
# Calculate IBS statistics
dt[, .(n=.N), key=.(accession, d=abs(wgsGT - gbsGT))] -> ab
dcast(ab, accession ~ d, value.var='n', fill=0) -> ab

# Rename columns
setnames(ab, paste(0:2), c("ibs2", "ibs1", "ibs0"))
ab[, n := ibs2+ibs0+ibs1]
ab[, p := ibs2 / (ibs2+ibs0)][] -> z  # rate of homozygous difference

# Sort and save
setorder(z, p)
write.table(z, file="WGS_GBS_IBS.tsv", sep="\t", quote=F, row.names=F)
```

---

### 7. Find nearest neighbour for mixed-up samples

#### 7.1 Define R function

```r
seq_NN <- function(gt, sample.id, gbs_mat, cores=1) {
  w <- gbs_mat
  
  setorder(rbindlist(mclapply(mc.cores=1, sample.id, function(idx) {
    gt[idx, on="sample.id", .(variant.id, wgs_sample=sample.id, wgsGT)] -> b
    w[b, allow.cartesian=T, nomatch=0, on="variant.id"][,
      .(n=.N), key=.(sample.id, wgs_sample, d=abs(gbsGT - wgsGT))] -> ab
    dcast(ab, sample.id + wgs_sample ~ d, value.var='n', fill=0) -> ab
    setnames(ab, paste(0:2), c("ibs2", "ibs1", "ibs0"))
    ab[, n := ibs2+ibs0+ibs1]
    ab[, p := ibs2 / (ibs2+ibs0)]
  })), -p) -> o
}
```

#### 7.2 Check input data format

The column names for the two datasets should match the examples below:

**GBS data format**:
```
     variant.id accession     sample.id gbsGT
1: chr1H:147087 HOR_10886 HOR_10886_GBS     2
2: chr1H:158716 HOR_10886 HOR_10886_GBS     0
3: chr1H:253058 HOR_10886 HOR_10886_GBS     0
4: chr1H:253403 HOR_10886 HOR_10886_GBS     0
5: chr1H:370764 HOR_10886 HOR_10886_GBS     0
6: chr1H:940241 HOR_10886 HOR_10886_GBS     2
```

**WGS data format**:
```
     variant.id accession     sample.id wgsGT
1: chr1H:147049 HOR_10886 HOR_10886_WGS     0
2: chr1H:147068 HOR_10886 HOR_10886_WGS     0
3: chr1H:147087 HOR_10886 HOR_10886_WGS     2
4: chr1H:147434 HOR_10886 HOR_10886_WGS     2
5: chr1H:158715 HOR_10886 HOR_10886_WGS     0
6: chr1H:158716 HOR_10886 HOR_10886_WGS     0
```

#### 7.3 Identify the genetically closest genotype for single mix-up sample

```r
# Find the closest genotype to sample "HOR_10886_WGS"
seq_NN(gt=wgs, sample.id='HOR_10886_WGS', gbs_mat=gbs) -> nn
nn[n > 1000][!duplicated(wgs_sample)] -> y
```

#### 7.4 Identify the genetically closest genotype for many mix-up samples

```r
# Read sample list
fread("WGS_NN_check_list.txt", header=F) -> lst
```

**Sample list example**:
```
              V1
1: HOR_10886_WGS
2: HOR_11013_WGS
3: HOR_11431_WGS
4: HOR_11922_WGS
5: HOR_12184_WGS
6: HOR_12791_WGS
```

```r
# Batch processing
rbindlist(mclapply(mc.cores=10, lst$V1, function(i) {
  seq_NN(gt=wgs, sample.id=i, gbs_mat=gbs) -> nn
  nn[n > 1000][!duplicated(wgs_sample)]
})) -> yy

# Save results
write.table(yy, file="wgs_mixup_NN.tsv", sep="\t", quote=F, row.names=F)
```

---

## Output Files

| File Name | Description |
|-----------|-------------|
| `GBS_SNPs.bed` | SNP interval file |
| `GBS_SNPs.fasta` | k-mer sequence file |
| `GBS_SNPs_uniq.fasta` | Unique k-mer sequence file |
| `het_estimation_result.txt` | Heterozygosity estimation results |
| `WGS_GBS_IBS.tsv` | IBS statistics results |
| `wgs_mixup_NN.tsv` | Sample mix-up nearest neighbour analysis results |

---

## Notes

1. File paths need to be adjusted according to your actual directory structure
2. Sample name extraction rules need to be adjusted according to your actual file naming format
3. R scripts depend on `data.table` and `parallel` packages
4. It is recommended to run on a multi-core server to improve processing speed

---

## Acknowledgments

This repository is cloned from the original project:
**https://bitbucket.org/ipk_dg_public/gbs_kmers/src/main/**
