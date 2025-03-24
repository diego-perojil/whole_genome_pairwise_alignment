#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1)
  stop("Usage: pairwise_alignment.R 
        <assembly_fa>", call. = FALSE)

# Assign arguments to variables
assembly_fa <- args[1]

# Load necessary library
# Make sure the library 'CNEr' (or any other required libraries) is installed.
library(CNEr)

# Paths
reference_root <- tools::file_path_sans_ext(basename(assembly_fa))
lastdb_path <- reference_root

# Build the lastdb index
system2(command = "lastdb", args= c("-c","-P 4",  lastdb_path, assembly_fa))