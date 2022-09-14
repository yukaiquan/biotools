# RNAseq Script



The script comes from https://www.genek.cn/, thanks to the big guy for providing such a useful tool.

# Install

install R and mamba

```shell
conda activate base
conda install -c conda-forge mamba
mamba create -n R
mamba activate R
```

install Rsubread limma edgeR

```shell
mamba install -c bioconda bioconductor-rsubread
mamba install -c conda-forge r-argparser
mamba install -c bioconda bioconductor-limma
mamba install -c bioconda bioconductor-edger
```

if edgeR is fail ,you can input R

```R
# install.packges('BiocManager')
library(BiocManager)
BiocManager::install('edgeR')
```

# run-featurecounts

```shell
$ Rscript run-featurecounts.R --help
usage: run-featurecounts.R [--] [--help] [--opts OPTS] [--bam BAM]
       [--gtf GTF] [--threads THREADS] [--output OUTPUT]

run featureCounts and calculate FPKM/TPM

flags:
  -h, --help     show this help message and exit

optional arguments:
  -x, --opts     RDS file containing argument values
  -b, --bam      input: bam file
  -g, --gtf      input: gtf file
  -t, --threads  threads
  -o, --output   output prefix
  
$ Rscript run-featurecounts.R -b test1.bam -g ACD.test.gtf -t 1 -o test1
```



# run-merge matrix

```shell
$ perl script/abundance_estimates_to_matrix.pl -h

############################################################
#
# Usage:  script/abundance_estimates_to_matrix.pl --est_method <method>  sample1.results sample2.results ...
#
#      or  script/abundance_estimates_to_matrix.pl --est_method <method> --quant_files file.listing_target_files.txt
#
#      Note, if only a single input file is given, it's expected to contain the paths to all the target abundance estimation files.
#
# Required:
#
#  --est_method <string>           featureCounts|RSEM|eXpress|kallisto|salmon  (needs to know what format to expect)
#
# Options:
#
#  --cross_sample_norm <string>         TMM|UpperQuartile|none   (default: TMM)
#
#  --name_sample_by_basedir             name sample column by dirname instead of filename
#      --basedir_index <int>            default(-2)
#
#  --out_prefix <string>                default: 'matrix'
#
#  --quant_files <string>              file containing a list of all the target files.
#
############################################################
perl script/abundance_estimates_to_matrix.pl --est_method featureCounts --out_prefix test_matrix --quant_files input_bam.txt
```

--quant_files input_bam.txt 

content:

```
test1.bam
test2.bam
test3.bam
```

