#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3)
  stop("Usage: pairwise_alignment.R 
        <assembly_fa> <query_fa> <output_folder>", call. = FALSE)

# Assign arguments to variables
assembly_fa <- args[1]
query_fa <- args[2]
output_folder <- args[3]

# Load necessary library
# Make sure the library 'CNEr' (or any other required libraries) is installed.
library(CNEr)


# Paths

reference_root <- tools::file_path_sans_ext(basename(assembly_fa))
query_root <- tools::file_path_sans_ext(basename(query_fa))
lastdb_path <- reference_root



# Build the lastdb index
system2(command = "lastdb", args= c("-c",  lastdb_path, assembly_fa))

# Run last aligner

output_maf <- paste0(reference_root, "_", query_root, ".maf")
lastal(db = reference_root, queryFn = query_fa,
       outputFn = output_maf,
       distance = "near", binary = "lastal", mc.cores = 1L)

# maf to psl conversion
output_psl <- paste0(reference_root, "_", query_root, ".psl")
system2(command = "maf-convert", args = c("psl", output_maf), output_psl)
