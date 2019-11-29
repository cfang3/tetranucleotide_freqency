# Read Cameron's signalP output (from???)

# Load tidyverse packages, including dplyr, readr, and ggplot2
library(seqinr)
library(tidyverse)

# Load the data
fn <- "data/GCF_000003645.1_ASM364v1_genomic.aa.fsa.out.neg"


# Load teh calc_tetranucleotide_freqs function
source("R/calc_tetra_nuc.R")

# Run one at a time
system.time({
  test_obj <- calc_tetranucleotide_freqs(fn, discard.data = FALSE)
})


# Get a character vector of all the signalP output files
all.filenames <- paste0("data/", dir("data/"))

# Run this for every signalP file you have
vapply(all.filenames, calc_tetranucleotide_freqs, save.data = TRUE, discard.data = TRUE)

