#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    library(data.table)
    library(parallel)
    library(argparse)
})

# Parse command line arguments
parser <- ArgumentParser(description = "WGS sample mix-up detection CLI tool (based on Nearest Neighbor algorithm)")

# Required arguments
parser$add_argument("-i", "--input",
    required = TRUE,
    help = "Input sample list file (format: one sample ID per line, no header)"
)
parser$add_argument("-d", "--data",
    required = TRUE,
    help = "K-mer data file path (.RData format, containing wgs, gbs, and other objects)"
)

# Optional arguments
parser$add_argument("-o", "--output",
    default = "wgs_mixup_NN.tsv",
    help = "Output result file path (default: wgs_mixup_NN.tsv)"
)
parser$add_argument("-t", "--threads",
    type = "integer", default = 64,
    help = "Number of parallel threads (default: 64, recommended not to exceed CPU cores)"
)
parser$add_argument("-m", "--min_n",
    type = "integer", default = 1000,
    help = "n value filter threshold (default: 1000)"
)

args <- parser$parse_args()

# Validate input file existence
if (!file.exists(args$input)) stop(paste("Input file does not exist:", args$input))
if (!file.exists(args$data)) stop(paste("Data file does not exist:", args$data))

# Load data
cat("Loading data file...\n")
load(args$data)

# Validate that required data objects exist
required_objects <- c("wgs", "gbs", "seq_NN")
missing_objects <- setdiff(required_objects, ls())
if (length(missing_objects) > 0) {
    stop(paste("Missing required objects in data file:", paste(missing_objects, collapse = ", ")))
}

# Read sample list
cat("Reading sample list...\n")
lst <- fread(args$input, header = FALSE)$V1 # Extract sample ID column
if (length(lst) == 0) stop("Sample list file is empty")

# Configure parallel parameters
threads <- min(args$threads, detectCores()) # Limit threads to not exceed CPU cores
cat(paste("Using", threads, "threads for parallel processing...\n"))

# Core processing logic
cat("Starting to process samples (total:", length(lst$V1), ")...\n")
yy <- rbindlist(mclapply(
    X = lst,
    mc.cores = threads,
    FUN = function(i) {
        # Call NN function
        nn <- seq_NN(gt = wgs, sample.id = i, gbs_mat = gbs)
        # Filter and remove duplicates
        nn[n > args$min_n]
        nn[n > args$min_n][!duplicated(wgs_sample)]
    }
))

# Output results
cat("Saving results to", args$output, "...\n")
write.table(
    yy,
    file = args$output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("Processing completed!\n")
