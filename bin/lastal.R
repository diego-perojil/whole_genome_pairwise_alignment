#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3)
  stop("Usage: pairwise_alignment.R 
        <assembly_fa> <query_fa> <output_folder> <dist_val>", call. = FALSE)

# Assign arguments to variables
assembly_fa <- args[1]
query_fa <- args[2]
dist_val <- args[3]

# Load necessary library
# Make sure the library 'CNEr' (or any other required libraries) is installed.
library(CNEr)

# Paths
reference_root <- tools::file_path_sans_ext(assembly_fa)
query_root <- tools::file_path_sans_ext(basename(query_fa))

# Run last aligner
output_maf <- paste0(reference_root, "_", query_root, ".maf")
lastal(db = reference_root, queryFn = query_fa,
       outputFn = output_maf,
       distance = dist_val, binary = "lastal", mc.cores = 100L) # mc.cores is 100 to make sure lastal uses the maximum possible number of cores, but will only ever use a number of cores = to number of contigs aligned. In the future, use process labels to efficiently manage core usage by lastal

# maf to psl conversion
output_psl <- paste0(reference_root, "_vs_", query_root, ".psl")
system2(command = "maf-convert", args = c("psl", output_maf), output_psl)
