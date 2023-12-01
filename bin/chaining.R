#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments
if (length(args) != 3) {
    stop("Usage: process_psl_to_chain.R psl_file 
          assemblyTarget assemblyQuery", call. = FALSE)
}

# Assign arguments to variables
psl_file <- args[1]
assemblyTarget <- args[2]
assemblyQuery <- args[3]

# Load the CNEr library
library(CNEr)

# Join close alignments
chains <- axtChain(psl_file, assemblyTarget = assemblyTarget,
                   assemblyQuery = assemblyQuery, distance = "near",
                   removePsl = FALSE, binary = "axtChain")

# Sort and combine
chainMergeSort(chains, assemblyTarget, assemblyQuery,
                           allChain = paste0(sub("\\.2bit$", "",
                                                 basename(assemblyTarget),
                                                 ignore.case=TRUE),
                                             "_",
                                             sub("\\.2bit$", "",
                                                 basename(assemblyQuery),
                                                 ignore.case = TRUE),
                                             ".all.chain"),
                           removeChains = FALSE,
                           binary = "chainMergeSort")
