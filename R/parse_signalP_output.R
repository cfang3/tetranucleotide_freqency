# Read signalP output and call the function "calc_tetra_nuc.R" on a whole list of files 

# Load tidyverse packages, including dplyr, readr, and ggplot2
library(seqinr)
library(tidyverse)

# Load .out.neg files 
signalP.output <- list.files(path = "data", recursive = TRUE, pattern = "out.neg$")

# Load the calc_tetranucleotide_freqs function
source("R/calc_tetra_nuc.R")

# Run one at a time to make function is working properly
fn <- "data/GCF_000003135.1_ASM313v1_genomic.aa.fsa.out.neg"
system.time({
  test_obj <- calc_tetranucleotide_freqs(fn, discard.data = FALSE)
})

# Get a character vector of all paths to the signalP output files (.neg files)
signalP.output.path <- paste0("data/", signalP.output)

# Run this for every signalP file you have and applies the function over the list of files.
lapply(signalP.output.path, calc_tetranucleotide_freqs, save.file = TRUE, discard.data = FALSE)
