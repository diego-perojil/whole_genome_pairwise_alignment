#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments
if (length(args) != 4) {
    stop("Usage: process_psl_to_chain.R psl_file 
          assemblyTarget assemblyQuery", call. = FALSE)
}

# Assign arguments to variables

netSyntenicFile <- args[1]
allPreChain <- args[2]
assemblyTarget <- args[3]
assemblyQuery <- args[4]

# Load the CNEr library
library(CNEr)

netToAxt(netSyntenicFile, allPreChain, assemblyTarget, assemblyQuery,
         axtFile = paste0(sub("\\.2bit$", "",
                              basename(assemblyTarget),
                              ignore.case=TRUE),
                        "_",
                        sub("\\.2bit$", "",
                            basename(assemblyQuery),
                            ignore.case = TRUE),
                        ".net.axt"),
         removeFiles = FALSE,
         binaryNetToAxt = "netToAxt", binaryAxtSort = "axtSort")