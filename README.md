# biotools

Bioinformatics Toolkit - A collection of scripts for genome annotation, variant detection, and functional annotation

## 📚 Table of Contents

- [Synteny Analysis](#synteny-analysis)
- [BSA/BSR Analysis](#bsabsr-analysis)
- [Blast Tools](#blast-tools)
- [FASTA Sequence Processing](#fasta-sequence-processing)
- [GFF/GTF Processing](#gffgtf-processing)
- [Functional Annotation](#functional-annotation)
- [ncRNA Analysis](#ncrna-analysis)
- [Genome Annotation](#genome-annotation)
- [Variant Detection (GATK)](#variant-detection-gatk)
- [Other Tools](#other-tools)

---

## Synteny Analysis

### `syntenic/filterBlock.py`

**Function**: Genome synteny block filtering and density analysis  
**Purpose**: Calculate syntenic gene density in chromosome intervals  
**Input**: Synteny file, chromosome length file  
**Output**: Interval-based synteny density distribution  
**Features**: Supports 1Mb sliding window analysis, suitable for whole genome duplication studies

---

## BSA/BSR Analysis

### `BSA/snpIndex.py`

**Function**: BSA SNP-index calculation  
**Purpose**: Calculate sliding window SNP-index for BSA mapping analysis  
**Input**: SNP position file, chromosome length file  
**Output**: SNP-index distribution per window  
**Features**: Multi-threaded processing, customizable window size

### `BSR/bsrSeq.py`

**Function**: BSR-seq (Bulk Segregant RNA-seq) analysis pipeline  
**Purpose**: Complete BSR-seq workflow including QC, alignment, variant detection  
**Pipeline**: Fastp QC → Alignment → FeatureCounts → DEG → GATK variant calling  
**Features**: Integrated tools, highly automated

---

## Blast Tools

### `blast/blast_best.py`

**Function**: Extract best hit for each query from Blast results  
**Purpose**: Filter Blast results to retain best match per query  
**Input**: Blast outfmt6 file  
**Output**: Deduplicated best match results

### `blast/blastout2idstatus.py`

**Function**: Statistics of Blast result distribution across chromosomes  
**Purpose**: Analyze query sequence matching distribution on different chromosomes  
**Input**: Blast output file  
**Output**: Match count statistics per query per chromosome

### `blast/blastOutMerge.py`

**Function**: Merge bidirectional Blast results  
**Purpose**: Combine forward and reverse Blast results for reciprocal best hit analysis  
**Input**: Two Blast output files  
**Output**: Merged results file

### `blast/diamond_blast.py`

**Function**: Diamond Blast parallel workflow  
**Purpose**: Fast protein sequence alignment using Diamond  
**Features**: Multi-threading, self-blast mode, database building  
**Advantage**: 20,000x faster than traditional Blast

---

## FASTA Sequence Processing

### `fasta/02_splitepep.py`

**Function**: Large-scale FASTA file splitting  
**Purpose**: Split large protein sequence files into smaller files (45,000 sequences per file)  
**Application**: Suitable for parallel processing of large-scale sequence alignment

### `fasta/ex_pep_from_wgdi_gff.py`

**Function**: Extract longest transcript from WGDI results  
**Purpose**: Extract representative longest transcript sequence per gene based on GFF  
**Input**: GFF file, protein sequence file  
**Output**: Longest transcript sequences per gene

### `fasta/exq_Sub_gene.py`

**Function**: Extract gene sequences by subgenome  
**Purpose**: Extract genes from specific subgenomes based on chromosome naming  
**Input**: BED file, CDS/PEP sequences  
**Output**: Subgenome-specific sequence files

### `fasta/fastaN60.py`

**Function**: Calculate FASTA N60 statistics  
**Purpose**: Evaluate sequence assembly quality  
**Input**: FASTA file  
**Output**: N60 and other statistics

### `fasta/fix_fasta.py`

**Function**: FASTA format repair and standardization  
**Purpose**: Clean special characters from sequences, unify format  
**Processing**: Remove spaces, hyphens, terminators, etc.  
**Output**: Standardized FASTA file (60 chars per line)

### `fasta/longest_gene_get.py`

**Function**: Extract longest transcript per gene  
**Purpose**: Select longest representative sequence from multi-transcript genes  
**Features**: Supports various transcript naming formats  
**Application**: Gene family analysis, sequence preparation for functional annotation

### `fasta/rename_avena_fasta.py`

**Function**: Oat genome sequence renaming  
**Purpose**: Standardize oat genome sequence names  
**Application**: Species-specific sequence naming normalization

### `fasta/rename_Pangene_short.py`

**Function**: Pan-gene sequence name simplification  
**Purpose**: Shorten pan-genome sequence name identifiers  
**Application**: Sequence management in pan-genome analysis

---

## GFF/GTF Processing

### `gff/convertgff2augustgtf.py`

**Function**: Convert GFF to Augustus GTF format  
**Purpose**: Convert standard GFF to Augustus-compatible GTF format  
**Input**: GFF file  
**Output**: Augustus GTF file

### `gff/EVM_TEsorter_remove.py`

**Function**: Remove TE annotations from EVM results  
**Purpose**: Clean transposable element annotations from EvidenceModeler results  
**Input**: EVM annotation file, TE list  
**Output**: Cleaned gene annotation file

### `gff/gushr_gtf_2_gff.py`

**Function**: Convert GUSHR GTF to GFF3 format  
**Purpose**: Convert GUSHR prediction results to standard GFF3 format  
**Features**: Supports UTR addition, element sorting  
**Input**: GUSHR GTF, original GFF  
**Output**: Standard GFF3 with UTR

### `gff/lowcase_fasta_from_gff.py`

**Function**: Convert FASTA sequences to lowercase based on GFF  
**Purpose**: Mark specific regions (e.g., repeats)  
**Input**: FASTA file, GFF file  
**Output**: Marked sequence file

### `gff/sorted_gff.py`

**Function**: GFF file sorting  
**Purpose**: Sort GFF by chromosome position  
**Application**: GFF file standardization

### `gff/Transcript2GenefromGff.py`

**Function**: Generate gene GFF from transcript GFF  
**Purpose**: Convert transcript-level annotations to gene-level  
**Input**: Transcript GFF  
**Output**: Gene GFF

---

## Functional Annotation

Collection of tools for parsing and processing functional annotation results from multiple databases:

### `go/get_go_term_tsv.py`

Convert GO OBO file to TSV format

### `go/goEnrich.R`

GO enrichment analysis using R

### `kegg/kegg_kofam.py`

Process KEGG KofamScan results

### `kegg/kegg_k_map.py`

Map KO numbers to KEGG pathways

### `Interproscan/Interproscan2IPR.py`

Extract protein domain annotations from InterProScan

### `Interproscan/InterProscan2GO_ann.py`

Extract GO terms from InterProScan results

### `kog/kog_result_ex.py`

Add functional descriptions to KOG annotation results

### `ncbi/blastp_nr_result.py`

Parse NR database Blast results with functional annotations

### `ncbi/nr_fmt6_2annotation.py`

Convert NR Blast outfmt6 to annotation table

### `ncbi/merge_annotation.py`

Merge annotation results from multiple databases

**Supported Databases**: GO, KEGG, InterPro, KOG, NR

---

## ncRNA Analysis

### `ncRNA/CIRCexplorer2_merge_matrix.py`

**Function**: Merge CIRCexplorer2 results  
**Purpose**: Combine circular RNA detection results from multiple samples  
**Input**: List of CIRCexplorer2 output files  
**Output**: Circular RNA expression matrix  
**Features**: Calculate SRPBM normalized values

### `ncRNA/get_exp_gff.py`

**Function**: Filter GFF by expression level  
**Purpose**: Extract transcripts expressed in multiple samples  
**Input**: Expression matrix, GTF file  
**Output**: Filtered GTF file

### `ncRNA/find_id_from_map.py`

**Function**: ID mapping lookup  
**Purpose**: Find ID relationships based on mapping file  
**Input**: Mapping file, ID list  
**Output**: ID relationship table

### `ncRNA/tagetGene_get.py`

**Function**: Process miRNA target gene prediction results  
**Purpose**: Process outputs from miRNA target prediction tools  
**Supports**: miRanda, TargetScan, RNA22  
**Output**: miRNA-target gene relationship table

### `ncRNA/stringtie_find_new.py`

**Function**: Extract novel transcripts from StringTie results  
**Purpose**: Filter novel transcripts assembled by StringTie  
**Input**: StringTie GTF file  
**Output**: Novel transcript GTF file  
**Features**: Filter by class_code (u, x, etc.)

---

## Genome Annotation

### `genome/TE_annotation.py`

**Function**: Transposable Element (TE) annotation pipeline  
**Purpose**: Three-step whole-genome TE annotation workflow  
**Pipeline**:

1. MITE-Hunter for MITE detection
2. LTR_Finder + LTR_Harvest + LTR_Retriever for LTR library construction
3. RepeatModeler + RepeatMasker for genome-wide annotation  
   **Features**: Multi-threaded, supports chromosome-level splitting

### `genome/TR_annotation.py`

**Function**: Tandem Repeat (TR) annotation  
**Purpose**: Microsatellite marker development using GMATA and TRF  
**Pipeline**: GMATA annotation → TRF fine analysis → dta2bed conversion  
**Input**: Genome FASTA  
**Output**: TR site BED file  
**Application**: Molecular marker development, genetic map construction

### `genome/dta2bed.py`

**Function**: Convert TRF dat file to BED format  
**Purpose**: Convert Tandem Repeat Finder results to BED format  
**Input**: TRF .dat file  
**Output**: BED format file (with repeat unit information)

---

## Variant Detection (GATK)

### `gatk/extract_snps_cli.py`

**Function**: Extract specific SNPs from VCF  
**Purpose**: Extract specific variants from large VCF file based on SNP ID list  
**Input**: VCF file, SNP ID list  
**Output**: VCF file with target SNPs  
**Features**: Supports compressed formats, efficient indexing

### `gatk/gatk_gvcf_merge.py`

**Function**: GVCF merging pipeline  
**Purpose**: Merge multiple sample GVCF files  
**Pipeline**: CombineGVCFs → GenotypeGVCFs  
**Input**: Multiple GVCF files  
**Output**: Multi-sample VCF

### `gatk/gatk_vcf_filter.py`

**Function**: VCF variant filtering  
**Purpose**: Quality control using GATK VariantFiltration  
**Filter criteria**: QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum  
**Output**: Filtered high-quality variant set

### `gatk/gatk2vcfshell.py`

**Function**: GATK workflow wrapper script  
**Purpose**: Automated GATK variant detection workflow  
**Features**: Integrated command-line parameters, simplified invocation

### `gatk/gvcf2vcf_tools.py`

**Function**: GVCF to VCF toolkit  
**Purpose**: GVCF to VCF conversion and processing  
**Application**: Single-sample variant detection result processing

### `gatk/ngs_rmbam.py`

**Function**: BAM file cleaning tool  
**Purpose**: Remove PCR duplicates from BAM files  
**Pipeline**: MarkDuplicates → Index building  
**Input**: BAM file  
**Output**: Deduplicated BAM file

### `gatk/snpEff24Dtv.py`

**Function**: SNP effect prediction and annotation  
**Purpose**: Variant functional annotation using SnpEff  
**Output**: Prediction of variant effects on gene function

---

## Other Tools

### `Bed/splitSubGenomeFrombed.py`

**Function**: Split BED file by subgenome  
**Purpose**: Extract subgenome-specific regions from whole-genome BED file  
**Input**: BED file, subgenome identifier  
**Output**: Subgenome-specific BED file

### `image/cut_image.py`

**Function**: Batch image cropping and stitching  
**Purpose**: Process bioinformatics visualization result images  
**Features**: Image cropping, rotation, batch stitching  
**Application**: Synteny plots, Circos plot post-processing

### `excel/excelSheet2file.py`

**Function**: Excel sheet splitting  
**Purpose**: Export multiple Excel sheets as separate files  
**Input**: Excel file  
**Output**: Multiple independent table files

### `jcvi/jcvi_files_product.py`

**Function**: JCVI toolkit helper script  
**Purpose**: Prepare input files for MCscan synteny analysis  
**Pipeline**: GFF to BED → Extract IDs → Get sequences  
**Input**: GFF, CDS, PEP files  
**Output**: Input files for JCVI analysis

### `kmer/kmc_fastq.py`

**Function**: K-mer analysis pipeline  
**Purpose**: K-mer spectrum analysis using KMC  
**Application**: Genome size estimation, heterozygosity assessment

### `pfm/pfm_scan.py`

**Function**: Transcription factor binding site prediction  
**Purpose**: Scan promoter regions for TF binding sites  
**Input**: Promoter sequences, PFM matrix  
**Output**: Binding site prediction results

### `mysql/pymysql_comm.py`

**Function**: MySQL database connection tool  
**Purpose**: Wrapper for pymysql operations, simplified database access  
**Features**: Context manager support, automatic connection management

### `mysql/nr_mysql_getann.py`

**Function**: Retrieve annotations from MySQL NR database  
**Purpose**: Query local NR database for protein functional annotations  
**Input**: Protein ID  
**Output**: Functional annotation information

---

## 🚀 Quick Start

### Requirements

- Python 3.7+
- Main dependencies:
  - numpy, pandas (data processing)
  - pysam (VCF processing)
  - biopython (sequence processing)
  - tqdm (progress bar)
  - Pillow (image processing)

### Install Dependencies

```bash
pip install numpy pandas pysam biopython tqdm Pillow
```

---

## 📖 Usage Examples

### Example 1: BSA SNP-index Analysis

```bash
python BSA/snpIndex.py snp_position.txt chr_length.txt 1000000 output.txt 8 snp
```

### Example 2: Extract Longest Transcripts

```bash
python fasta/longest_gene_get.py -i input.fasta -o longest.fasta
```

### Example 3: TE Annotation Pipeline

```bash
python genome/TE_annotation.py -i genome.fa -o te_annotation -t 16
```

### Example 4: Extract Specific SNPs from VCF

```bash
python gatk/extract_snps_cli.py -i input.vcf -s snp_list.txt -o output.vcf
```

---

## 📧 Contact

Author: yukaiquan  
Email: 1962568272@qq.com  
Project under continuous development...

---

## 📄 License

This project is licensed under Apache License 2.0

---

## 🙏 Acknowledgments

Thanks to the developers of:

- JCVI/MCscan (synteny analysis)
- GATK (variant detection)
- InterProScan (functional annotation)
- Diamond (sequence alignment)
- RepeatModeler/RepeatMasker (repeat annotation)
