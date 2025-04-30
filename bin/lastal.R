#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3)
    stop("Usage: pairwise_alignment.R <assembly_fa> <query_fa> <dist_val>", call. = FALSE)

# Assign arguments to variables
assembly_fa <- args[1]
query_fa    <- args[2]
dist_val    <- args[3]

# Load necessary library
library(CNEr)

# Prepare paths
reference_root <- tools::file_path_sans_ext(assembly_fa)
query_root     <- tools::file_path_sans_ext(basename(query_fa))
output_maf     <- paste0(reference_root, "_", query_root, ".maf")

# Multi-threaded LAST wrapper
lastalMultiTh <- function(db,
                          queryFn,
                          outputFn = sub("\\.(fa|fasta)$", ".maf",
                                         paste(basename(db), basename(queryFn), sep = "_"),
                                         ignore.case = TRUE),
                          distance = c("far","medium","near"),
                          threads  = 100L,
                          binary   = "lastal",
                          echoCommand = FALSE) {
  distance <- match.arg(distance)
  if (tools::file_ext(outputFn) != "maf")
    stop("The outputFn must end in .maf", call. = FALSE)

  # Build substitution matrix
  matrixFile <- tempfile(fileext = ".lastzMatrix")
  on.exit(unlink(matrixFile))
  write.table(scoringMatrix(distance),
              file      = matrixFile,
              quote     = FALSE,
              sep       = " ",
              row.names = TRUE,
              col.names = TRUE)

  # Distance‐specific options
  lastOptions <- list(
    near   = paste("-a 600 -b 150 -e 3000 -p", matrixFile, "-s 2"),
    medium = paste("-a 400 -b  30 -e 4500 -p", matrixFile, "-s 2"),
    far    = paste("-a 400 -b  30 -e 6000 -p", matrixFile, "-s 2")
  )

  # Add -P<threads>
  opts <- paste(lastOptions[[distance]], paste0("-P", threads))

  # Build and run the command
  cmd <- paste(binary, opts, "-f 1", db, queryFn, ">", outputFn)
  message("Running lastalMultiTh with -P", threads, " threads...")
  if (echoCommand) {
    message(cmd)
  } else {
    status <- system(cmd, intern = FALSE)
    if (status != 0) stop("lastal failed with exit code ", status)
  }

  invisible(outputFn)
}

# Invoke the wrapper
lastalMultiTh(db       = reference_root,
              queryFn  = query_fa,
              outputFn = output_maf,
              distance = dist_val)

# Convert MAF to PSL
output_psl <- paste0(reference_root, "_vs_", query_root, ".psl")
system2("maf-convert", args = c("psl", output_maf), stdout = output_psl)
