#################################################R#################################################

library(data.table)
library(parallel)

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

# Define all directories for batch reading
target_dirs <- c(
  "../../01_fastq/GBS666",
)

# Batch collect all _count.tsv.gz file paths from all directories
gbs666_f <- rbindlist(lapply(target_dirs, function(dir) {
  # Build find command to search for all _count.tsv.gz files in specified directory
  find_cmd <- paste0("find ", shQuote(dir), " | grep '_count.tsv.gz$'")
  # Read file paths
  fread(find_cmd, head=F, col.names="path") -> tmp
  return(tmp)
}))

# Extract sample name accession (following original logic: remove path prefix and suffix)
# Optimized to be compatible with all directories: extract filename after last "/", then remove suffix
gbs666_f[, accession := sub("_count.tsv.gz$", "", basename(path))]

# 2. Remove everything after "_" -> get names like "sanfensan", "SRR32111801"
gbs666_f[, filename := sub("_.*", "", accession)]
# 3. Add _gbs666 suffix -> final result like "sanfensan_gbs666", "SRR32111801_gbs666"
gbs666_f[, accession := paste0(filename, "_gbs666")]
# Remove intermediate temporary column (optional)
gbs666_f[, filename := NULL]

length(gbs666_f$accession)

# Optional: remove duplicates (avoid duplicate files with same name in different directories)
gbs666_f <- unique(gbs666_f, by="accession")

write.table(gbs666_f, "all_samples.tsv", sep="\t", row.names=F, quote=F)

# Batch read and calculate genotypes (cores can be adjusted based on CPU cores, e.g., cores=8)
gbs666_gt <- read_kmer_count(info=gbs666_f, mincov=1, cores=64)

# Organize output results (keep needed columns)
gbs666 <- gbs666_gt[, .(variant.id, accession, gbs666GT=gt)]

# Optional: check basic result information
cat("Successfully read", length(unique(gbs666$accession)), "samples\n")
head(gbs666)


# Read GBS9000 data
fread("find ../01_GBS | grep '_gbs_count.tsv.gz$'", head=F, col.names="path") -> gbs9000_f
gbs9000_f[, accession := sub("_gbs_count.tsv.gz$", "", sub("^.*/", "", path))]  # needs to be adjusted according to user's directory and file name
gbs9000_f[, accession := paste0(accession, "_gbs9000")]

read_kmer_count(info=gbs9000_f, mincov=1, cores=64) -> gbs9000_gt
gbs9000_gt[, .(variant.id, accession, gbs9000GT=gt)] -> gbs9000

# Optional: check basic result information
cat("Successfully read", length(unique(gbs9000$accession)), "samples\n")
head(gbs9000)


# combine two genotype tables and check correlation
# gbs666[gbs9000, on=c("variant.id", "accession"), nomatch=0] -> dt
# dt[, .(Cor = cor(gbs666GT, gbs9000GT)), by=accession][order(Cor)] -> dtCor  # check genotype correlation between GBS9000 and GBS666 datasets


# Find nearest neighbour for mixed-up samples
seq_NN <- function(gt, sample.id, gbs_mat, cores=1) {
    w <- gbs_mat

    setorder(rbindlist(mclapply(mc.cores=1, sample.id, function(idx) {
        gt[idx, on="sample.id", .(variant.id, gbs666_sample=sample.id, gbs666GT)] -> b
        w[b, allow.cartesian=T, nomatch=0, on="variant.id"][,
        .(n=.N), key=.(sample.id, gbs666_sample, d=abs(gbs9000GT - gbs666GT))] -> ab
        dcast(ab, sample.id + gbs666_sample ~ d, value.var='n', fill=0) -> ab
        setnames(ab, paste(0:2), c("ibs2", "ibs1", "ibs0"))
        ab[, n := ibs2+ibs0+ibs1]
        ab[, p := ibs2 / (ibs2+ibs0)]
    })), -p) -> o
}

gbs9000[, sample.id := accession]
gbs666[, sample.id := accession]


head(gbs9000)
    variant.id               accession gbs9000GT               sample.id
          <char>                  <char>     <num>                  <char>
1: chr1A:1686299 TINKER2022_CN_20075_gbs9000         0 TINKER2022_CN_20075_gbs9000
2: chr1A:1686300 TINKER2022_CN_20075_gbs9000         0 TINKER2022_CN_20075_gbs9000
3: chr1A:1686307 TINKER2022_CN_20075_gbs9000         0 TINKER2022_CN_20075_gbs9000
4: chr1A:1686317 TINKER2022_CN_20075_gbs9000         0 TINKER2022_CN_20075_gbs9000
5: chr1A:1686360 TINKER2022_CN_20075_gbs9000         0 TINKER2022_CN_20075_gbs9000
6: chr1A:1805389 TINKER2022_CN_20075_gbs9000         0 TINKER2022_CN_20075_gbs9000

head(gbs666)
      variant.id accession gbs666GT sample.id
          <char>    <char>     <num>    <char>
1: chr1A:1686299  sanfensan_gbs666     0  sanfensan_gbs666
2: chr1A:1686300  sanfensan_gbs666     0  sanfensan_gbs666
3: chr1A:1686307  sanfensan_gbs666     0  sanfensan_gbs666
4: chr1A:1686317  sanfensan_gbs666     0  sanfensan_gbs666
5: chr1A:1686360  sanfensan_gbs666     0  sanfensan_gbs666
6: chr1A:1805389  sanfensan_gbs666     0  sanfensan_gbs666


# Identify the genetically closest genotype for single mix-up sample
seq_NN(gt=gbs666, sample.id='sanfensan_gbs666', gbs_mat=gbs9000) -> nn  # this code will search the closest genotype to the sample "sanfensan_gbs666"
nn[n > 1000][!duplicated(gbs666_sample)] -> y

# The following steps are time-consuming and can be split into HPC for acceleration
save(gbs666, gbs9000, seq_NN, file = "Kmer.RData")

#################################################shell#################################################

# run in linux 
# awk '{print $2}' all_samples.tsv | grep -v "accession" > GBS666_samples.txt
# Identify the genetically closest genotype for many mix-up samples

Rscript kmerNN.r -i samples.txt -d Kmer.RData -o gbs666_mixup_NN.tsv -t 64
# ls xa* | while read id; do echo -e "$id\tRscript kmerNN.r -i ${id} -d Kmer.RData -o ${id}_mixup_NN.tsv -t 120"; done > sbatch.tsv
# python /public/share/h13713/soft/sbatch_script.py -i sbatch.tsv -p com300 -m 500G -N 1 -n 120 -e python2.7
# https://github.com/yukaiquan/biotools/blob/main/slurm/sbatch_script.py