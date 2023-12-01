#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments
if (length(args) != 3) {
    stop("Usage: process_psl_to_chain.R psl_file 
          assemblyTarget assemblyQuery", call. = FALSE)
}

# Assign arguments to variables
allChain <- args[1]
assemblyTarget <- args[2]
assemblyQuery <- args[3]

# Load the CNEr library
library(CNEr)

## Filtering out chains
allPreChain <- chainPreNet(allChain, assemblyTarget, assemblyQuery,
                           allPreChain = paste0(sub("\\.2bit$", "",
                                                  basename(assemblyTarget),
                                                  ignore.case=TRUE),
                                             "_",
                                             sub("\\.2bit$", "",
                                                 basename(assemblyQuery),
                                                 ignore.case = TRUE),
                                                 ".all.pre.chain"),
                           removeAllChain = FALSE, binary = "chainPreNet")

## Keep the best chain and add synteny information
chainNetSyntenic(allPreChain, assemblyTarget, assemblyQuery,
                     netSyntenicFile=paste0(sub("\\.2bit$", "",
                                                basename(assemblyTarget),
                                                ignore.case=TRUE),
                                             "_",
                                             sub("\\.2bit$", "",
                                                 basename(assemblyQuery),
                                                 ignore.case = TRUE),
                                               ".noClass.net"),
                     binaryChainNet="chainNet", binaryNetSyntenic="netSyntenic")
