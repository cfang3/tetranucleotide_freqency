Tetranucleotide Frequency and Principal Component Analysis
==========================================================

This is R Markdown document contains an explanation of a total of four
.R files of which three are functions (`pca()`, `calc_tetra_nuc()`, and
`get_gg_title()`). The last is a .R file that contains the necessary
packages and sources the beforementioned functions. The goal of the code
is to measure tetranucleotide frequencies of extracellular and
intracellular enzymes of an organism and perform a PCA to determine if
extracellular and intracellular enzymes are markedly different. This
should be run on a whole list of organisms from RefSeq.

Directory structure and the files that it contains
--------------------------------------------------

The directory structure starts with a directory called "/Documents."
Under that is a directory called "/steen\_lab"" which contains four
folders labeled "data," "plots," "misc" and "R". Under "data" contains 3
types of files: fna (has whole shotgun genome sequence), fsa (which
contains the proteins that the genomes encode), and .out.neg (signalP
output file that contains the sequences of which enzymes are
intracellular and extracellular)

Loading the necessary functions and packages
--------------------------------------------

This .R file loads the functions and packages that is used within the
three created functions. It also utlizes mapply and lapply to apply the
code on many file paths

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

Caclulating tetrancleotide frequencies
--------------------------------------

This part relies on applying the function `calc_tetra_nuc()` on a whole
list of file paths, using `lapply()`. This path name is the only
arguement for `calc_tetra_nuc()` and leads to files that contain the
outputs from signalP, a protein analysis program.

    signalP.output <- list.files(path = "data", recursive = TRUE, pattern = "out.neg$")

    signalP.output.path <- paste0("data/", signalP.output)

    signalP.output.path

    ## [1] "data/GCF_000003135.1_ASM313v1_genomic.aa.fsa.out.neg"
    ## [2] "data/GCF_000003645.1_ASM364v1_genomic.aa.fsa.out.neg"
    ## [3] "data/GCF_000003925.1_ASM392v1_genomic.aa.fsa.out.neg"

    lapply(signalP.output.path, calc_tetra_nuc)

### Inside `calc_tetra_nuc()`

The calc\_tetra\_nuc function takes the .out.neg files, fna files and
fsa file and reads them.

    calc_tetra_nuc <- function(signalP.out.fn, save.file = TRUE, discard.data = FALSE) {
      

      # Read the SignalP output file = signalP.out.fn
      colnames <- c("name", "Cmax", "Cmax.pos", "Ymax", "Ymax.pos", "Smax", "Smax.pos", "Smean", "D", "prediction", "Dmaxcut", "networks.used")
      signalP_results <- read_fwf(file = signalP.out.fn, col_positions = fwf_empty(signalP.out.fn, skip = 1, col_names = colnames), skip = 2)
      
      # Strip the excess information off the sequence name to pull out the root of the fasta file name to get the paths of the fna and fsa file
      fsa.path <- str_extract(signalP.out.fn, ".+(?=(\\.out\\.neg)$)") # regular expression: match one or more characters of everything before .out.neg
      fna.path <- str_replace(fsa.path, "(aa\\.fsa)$", "fna") # regular expression: take characters aa and every character following fsa to the end of the string and replace it with fna
      
      # Read two different fasta files: the fna file, which has DNA sequences with multiple genes (whole genome shotgun sequence) 
      # AND the fsa file, with amino acid sequences of single genes (which is what we fed to signalP)
      curr_fna <- read.fasta(fna.path) # curr_fasta means "current fasta" because we'll call this function on a whole list of fasta files
      curr_fsa_file <- seqinr::read.fasta(fsa.path)  # Open the amino acid sequence file that has the start and end address of the DNA sequence 

Next, the function makes empty collumns to fill in later with
tetranucleotide frequencies

      # For each sequence, calculate tetranucleotide frequencies
      freq_list <- vector("list", nrow(signalP_results)) # "pre-allocate" a list a number of elements equal to the number of rows of the SignalP output
      freq_results_names <- vector("character", nrow(signalP_results))

The function then takes every row in the signalP output file and takes
the name of the protein sequence. It then extracts the attribute "Annot"
from the fsa file protein name. This attribute contains the beginning
and end numbers of the nucleotide sequence that corresponds to the
protein name.

      for(i in 1:nrow(signalP_results)) {
        
        # Extract the current sequence name from the signalP output; get rid the underscore and everything after it
        curr.aa.seq.name <- as.character(signalP_results[i, "name"]) # annoyingly, if we don't wrap signalP_results[i, "name"] in as.character(), we'll get back a tibble 
        
        # Find the header of that sequence in the curr_fsa_file
        curr.fsa.header <- attr(curr_fsa_file[[curr.aa.seq.name]], "Annot")
        
        # Extract the two numbers after the two # signs; these are the first and last addresses of the gene sequence that we're looking at
        curr.dna.start.address <- as.numeric(str_extract(curr.fsa.header, "(?<=#\\s)\\d+(?=\\s#)")) # regular expression: match everything behind the first # and space and dont include it; followed by any number of digits; match everything after the next space and # and dont include it
        curr.dna.end.address <- as.numeric(str_extract(curr.fsa.header, "(?<=#\\s\\d{1,10}\\s#\\s)\\d+(?=\\s#)")) # regular expression: look behind the first #, space, number, space, #, space and dont include it; followed by any amount of digits; then match everything after space and # and dont include it

Using those start and end numbers, the gene sequence, containing many
nuclotides, is extracted and put into "curr.gene.seq." Then using a
package called seqinr, the tetranucleotide frequencies are calculated.

        # Get the DNA sequence we're dealing with & extract the gene sequence that corresponds to the AA sequence that signalP looked at
        curr.dna.seq.name <- str_extract(curr.aa.seq.name, ".+(?=(_(\\d+)$))")
        
        # Extract the actual DNA gene sequence from the DNA fna file, based on the name of the sequence and its start and end address
        curr.dna.seq <- curr_fna[[curr.dna.seq.name]]
        curr.gene.seq <- curr.dna.seq[curr.dna.start.address : curr.dna.end.address]  # gets the actual nuclotide sequence based on the start and end address in which the tetranucleotide frequency will be measured
        curr.tet.freqs <- seqinr::count(curr.gene.seq, alphabet = , wordsize = 4) # calculate the frequencies 

The tetranucleotide frequencies are put in a list and added as a column
to the signalP output data table. These R objects are then saved into a
Rdata file in the /data folder.

        # Put the tetranucleotide freqencies in the list, and put the name of the gene in the vector of names
        freq_list[[i]] <- list(name = curr.aa.seq.name, 
                               freq.table = curr.tet.freqs)
        
      }
      
      # Convert the list of tetranucleotide frequencies into a data frame with a list-column and join with the signalP output
      freq_df <- freq_list %>%
        tibble(name = map_chr(., "name"),
               freqs = map(., "freq.table")) %>%
        select(-.)
      signalP_with_freqs <- signalP_results %>%
        left_join(freq_df, by = "name")
      
      # Save the file in the current directory and output a message for when it is saved
      if(save.file) {
        save.file.fn <- paste0(signalP.out.fn, ".Rdata")
        save(signalP_with_freqs, file = save.file.fn) # note that it will save the file in whatever directory the signalP output file was in
        message(paste("Saving output data frame as", save.file.fn))
      }
      
      if(discard.data) {
        return(NULL)
      } else {
        return(freq_df)
      }
      
    }

PCA output
----------

This part relies on the function `pca()` being applied on a list of two
types of files (Rdata files and signalP output files) using `mapply()`.
The function `pca()` itself relies on two arguements: **1**. The file
path to the .Rdata file **2**. The file path to the signalP output file.

    Rdata.fns <- paste0("data/", dir("data/", "\\.Rdata$"))

    Rdata.fns

    ## [1] "data/GCF_000003135.1_ASM313v1_genomic.aa.fsa.out.neg.Rdata"
    ## [2] "data/GCF_000003645.1_ASM364v1_genomic.aa.fsa.out.neg.Rdata"
    ## [3] "data/GCF_000003925.1_ASM392v1_genomic.aa.fsa.out.neg.Rdata"

    signalP.output.path

    ## [1] "data/GCF_000003135.1_ASM313v1_genomic.aa.fsa.out.neg"
    ## [2] "data/GCF_000003645.1_ASM364v1_genomic.aa.fsa.out.neg"
    ## [3] "data/GCF_000003925.1_ASM392v1_genomic.aa.fsa.out.neg"

    graph_list <- mapply(Rdata.fns, FUN = pca, signalP.out.fn = signalP.output.path)

### Inside `pca()`

The function `pca()` itself relies on two arguements: **1**. The file
path to the .Rdata file **2**. The file path to the signalP output file.

    pca <- function(signalP.freqs, signalP.out.fn, save.file = TRUE, discard.data = FALSE)

This function takes the first arguement, signalP.freqs = Rdata.fns, and
loads it into the global environment. An empty matrix is formed equal to
the number of rows from signalP.freqs with 256 columns (256 because that
is the number of tetranuceotide rearrangements that can occur, 4^4)

    # Load signalP.freqs into global environment
      load(signalP.freqs)
      
      # Need to combine all of these into one matrix
      # First make an empty matrix
      tet.freqs <- matrix(data = NA_integer_, nrow = nrow(signalP_with_freqs), ncol = 256)

Using a for loop, we fill the empty matrix with the tetranucleotide
frequencies

     # Use a loop to convert each tetranucleotide data table into a vector and put it into the matrix we'll use for PCA
      system.time({
        for(i in 1:nrow(signalP_with_freqs)) {
          tet.freqs[i, ] <- as.vector(signalP_with_freqs[i, "freqs"][[1]][[1]])
        }
      })

Take the names of the tetranucelotides (ie aaaa, aaac, etc.) and put
them as the column headers in the matrix that we created before

    # Now pull out the 'dimension names' for the tetranucleotide frequency table to get the tetranucleotides in question
      freq.names <- attr(signalP_with_freqs[1, "freqs"][[1]][[1]], "dimnames")[[1]]
      colnames(tet.freqs) <- freq.names

Standardize the matrix and use `rda()` from the vegan package to perform
the PCA

      # Standardize the matrix
      tet.freqs.std <- vegan::decostand(tet.freqs, method = "total", MARGIN = 1)
      
      # Use vegan::rda() to perform PCA and export data
      system.time({
        tet.freqs.std.rda <- vegan::rda(X=tet.freqs.std)
      })

Take the prediction values (tells us if a protein is extracelular vs
intracellular) and add it to a new dataframe that contains PCA values

    ## Graph the intracellular against extracellular enzymes with PCA values

    # Pull out the 'prediction' column of the data frame - this is the 'labels' we'll use to tell rda() which vectors are exported
      prediction <- signalP_with_freqs$prediction
      
      # Create data frame of PCA values and add prediction from SignalP output to compare intracellular and extracellular enzymes
      curr_PCA <- as.data.frame(tet.freqs.std.rda$CA$u) %>%
        mutate(exported = signalP_with_freqs$prediction)

Get the eigenvalues to eventually put into the ggplot graph to visualize
how much variance is captured

    # Calcualte variance explained
      eigenvalues <- tet.freqs.std.rda$CA$eig
      var.explained <- eigenvalues / sum(eigenvalues)

The beginnings of the ggplot are made by mapping aesthetics and
`get_gg_title()` is used to extract the organism name. This is
eventually used as the main title of the ggplot.

      # Create ggplot 
      p <- ggplot(curr_PCA, aes(x=PC1, y=PC2, colour = exported)) + 
        geom_point(alpha = 0.3)   # make the points a bit transparent 
      
       # Get organism name to use as ggplot title
      get_gg_title(signalP.out.fn)

The final ggplot is made with new x and y axis labels, new main title,
and color mapping between extracellular and intracellular enzymes.
Finally the plot is saved as a png in plots folder

      graph <- p +
        xlab(paste0("PC1, ", round(var.explained["PC1"], digits = 2)*100, "%")) +
        ylab(paste0("PC2, ", round(var.explained["PC2"], digits = 2)*100, "%")) + 
        scale_color_manual(values = c("red", "black")) + 
        ggtitle(ggplot.name) +
        theme_minimal() + 
        theme(text = element_text(size = 24)) 
      print(graph)
      
      # Save plot under plots folder
      ggsave(paste0("plots/", ggplot.name), device = "png", height = 8, width = 10, units = "in", dpi = 300)

### Inside `get_gg_title()`

This function takes the second arguement of `pca()`, signalP.out.fn =
signalP.output.path, and reads the fna file like in `calc_tetra_nuc()`.
It then extracts the organism name by taking the attribute "Annot" from
the first item in the list of the fna file. A regular expression is then
used to isolate the name.

    # Extract the organism name from the current fna file
    get_gg_title <- function(signalP.out.fn) {
      
      # Extract the root of the signalP.out.fn file name to get the fna path
      fsa.path <- str_extract(signalP.out.fn, ".+(?=(\\.out\\.neg)$)") # directory path to fsa file using stringr
      fna.path <- str_replace(fsa.path, "(aa\\.fsa)$", "fna") # directory path to fna file using stringr
      
      # Read the fna file that has the organism name under attr("Annot")
      curr_fna <- seqinr::read.fasta(fna.path) # curr_fasta means "current fasta" because we'll call this function on a whole list of fasta files
      
      # Extract the organism's name from the attr of the fna file to be used in the PCA plot
      ggplot.name <<- str_extract(attr(curr_fna[[1]],"Annot"),"(?<=\\s).{1,1000}(?=,\\s)")
    }

    source('~/Documents/steen_lab/R/get_gg_title.R')
    get_gg_title(signalP.output.path[[1]])
    ggplot.name

    ## [1] "Bifidobacterium longum subsp. longum ATCC 55813 SCAFFOLD1"

PCA Plots
---------

    Rdata.fns <- paste0("data/", dir("data/", "\\.Rdata$"))

    # Run PCA on all saved Rdata files that contain tetranucleotide frequencies
    system.time({
      graph_list <- mapply(Rdata.fns, FUN = pca, signalP.out.fn = signalP.output.path)
    })

![](Final_proj_files/figure-markdown_strict/unnamed-chunk-23-1.png)![](Final_proj_files/figure-markdown_strict/unnamed-chunk-23-2.png)![](Final_proj_files/figure-markdown_strict/unnamed-chunk-23-3.png)

    ##    user  system elapsed 
    ##  90.147   5.883 111.668

Analyzing the Results
---------------------

Based on the preliminary results in does not appear that extracellular
enzymes are markedly more diverse than non-secreted enzymes. There does
not seem to be a seperation of black points and red points indicating
such a difference. In addition, there is unfortuately not much variance
captured within the plot, with the highest amount of variance captured
being 27%.

Future Work
-----------

While the intital results seem to disprove the hypothesis that
extracellular enzymes are more diverse, this code needs to be run on
thousands of more organisms. One area that needs to be improved is that
there needs to be code added to do permanova. This will give us an
objective rating as to whether or not extracellular enzymes are that
much different than intracellular ones
