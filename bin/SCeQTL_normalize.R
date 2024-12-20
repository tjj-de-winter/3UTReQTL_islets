#!/usr/bin/env Rscript

### Description ###
# import a raw merged count table and normalize this using the default SCeQTL (according to DEseq method)

# Load necessary libraries
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("SCeQTL", quietly = TRUE)) devtools::install_github("XuegongLab/SCeQTL")

library(SCeQTL)

### Input variables ### 
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript normalize_counts.R <input_csv> <output_csv>", call. = FALSE)
}

input_file <- args[1]
output_file <- args[2]

### Code ###

# Import raw count data
raw_counts <- read.csv(input_file, row.names = 1)

# Apply normalization function to each gene
normalized_counts <- normalize(raw_counts)

# Save normalized data to a new CSV file
write.csv(normalized_counts, output_file, row.names = TRUE)