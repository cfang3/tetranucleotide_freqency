### Read signalP output and call the function "calc_tetra_nuc.R" on a whole list of files. 
### Use the tetranucleotide freqncies from a whole list of files and perfrom a pca on them.

## Calculate tetranucleotide frequencies on all files
# Load necessary packages
library(seqinr)
library(tidyverse)
library(vegan)

# Load necessary functions
source("R/calc_tetra_nuc.R")
source('~/Documents/steen_lab/R/get_gg_title.R')
source('~/Documents/steen_lab/R/pca.R')

# Load .out.neg files 
signalP.output <- list.files(path = "data", recursive = TRUE, pattern = "out.neg$")

# Run one at a time to make function is working properly
fn <- "data/GCF_000003135.1_ASM313v1_genomic.aa.fsa.out.neg"
system.time({
  test_obj <- calc_tetra_nuc(fn, discard.data = FALSE)
})

# Get a character vector of all paths to the signalP output files (.neg files)
signalP.output.path <- paste0("data/", signalP.output)

# Run this for every signalP file you have and applies the function over the list of files.
lapply(signalP.output.path, calc_tetra_nuc)

## Perform PCA on the tetranuclotide frequency Rdata files and save a graph of it.
# Get path for Rdata files
Rdata.fns <- paste0("data/", dir("data/", "\\.Rdata$"))

# Run PCA on all saved Rdata files that contain tetranucleotide frequencies
system.time({
  graph_list <- mapply(Rdata.fns, FUN = pca, signalP.out.fn = signalP.output.path)
})



