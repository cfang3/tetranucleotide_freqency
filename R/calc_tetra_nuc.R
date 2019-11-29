##' This takes the filename of a signalP output and returns a list of tetranucleotide frequencies of the associated DNA sequences

calc_tetranucleotide_freqs <- function(signalP.out.fn, save.file = TRUE, discard.data = FALSE) {
  
  
  
  
  # Strip the excess information off the sequence name to pull out the root of the fasta file name
  # Use regular expressions and the r package called stringr to do this
  # see https://www.rstudio.com/resources/cheatsheets/
  #fasta.core.name <- str_extract(seq.name, "^.{1,}_genomic\\.aa\\.fsa\\.out\\.neg$") # What ive tried to do is: "start of string, followed by 1 or more of any character, followed by a period
  # fasta.core.name <- str_extract(seq.name, "(?<=^(data/))(.+)(?=_genomic\\.aa\\.fsa\\.out\\.neg$)")

  # Read the SignalP output file = signalP.out.fn
  colnames <- c("name", "Cmax", "Cmax.pos", "Ymax", "Ymax.pos", "Smax", "Smax.pos", "Smean", "D", "prediction", "Dmaxcut", "networks.used")
  signalP_results <- read_fwf(file = fn, col_positions = fwf_empty(fn, skip = 1, col_names = colnames), skip = 2)
 # browser()
  # Read two different fasta files: the fna file, which has DNA sequences with multiple genes (either contigs or closed genomes, I think) 
  #    AND the faa file, with amino acid sequences of signle genes (which is what we fed to signalP)
  
  fsa.path <- paste0("data1/proteins/", proteins)
  fna.path <- paste0("data1/nucleotides/", nucleotides)
  

  curr_fna <- read.fasta("data1/nucleotides/GCF_000003135.1_ASM313v1_genomic.fna") # curr_fasta means "current fasta" because we'll call this function on a whole list of fasta files
  
  # Open the amino acid sequence file to get the start and end address of hte DNA sequence
  curr_fsa_file <- seqinr::read.fasta("data1/proteins/GCF_000003135.1_ASM313v1_genomic.aa.fsa")

  # # For each sequence, calculate tetranucleotide frequencies
  # # TESTING HERE: FOR JUST THE FIRST LINE
  # seq.name.deleteme <- fasta_results[1, "name"]
  # # Get the sequence from curr_fasta
  # test.freqs <- seqinr::count(curr_fasta[["NZ_GG666849.1_1"]], wordsize = 4)
  
  freq_list <- vector("list", nrow(signalP_results)) # "pre-allocate" a list a number of elements equal to the number of rows fo the 
  freq_results_names <- vector("character", nrow(signalP_results))
  for(i in 1:nrow(signalP_results)) {
    
    # Extract the current sequence name from the signalP output; get rid the underscore and everything after it
    curr.aa.seq.name <- as.character(signalP_results[i, "name"]) # annoyingly, if we don't wrap signalP_results[i, "name"] in as.character(), we'll get back a tibble (a kind of data frame)
    
    # Find the header of that sequence in the curr_fsa_file
    curr.fsa.header <- attr(curr_fsa_file[[curr.aa.seq.name]], "Annot")
    # Extract the two numbers after the two # signs; these are the first and last addresses of the gene sequence that we're looking at
    curr.dna.start.address <- as.numeric(str_extract(curr.fsa.header, "(?<=#\\s)\\d+(?=\\s#)"))
    curr.dna.end.address <- as.numeric(str_extract(curr.fsa.header, "(?<=#\\s\\d{1,10}\\s#\\s)\\d+(?=\\s#)"))
    
    # Get the DNA sequence we're dealing with & extract the gene sequence that corresponds to teh AA sequence that signalP looked at
    curr.dna.seq.name <- str_extract(curr.aa.seq.name, ".+(?=(_(\\d+)$))")
    
    # Extract the actual DNA gene sequence from the DNA fasta file, based on the name of the sequence and its' start and end address
    curr.dna.seq <- curr_fna[[curr.dna.seq.name]]
    curr.gene.seq <- curr.dna.seq[curr.dna.start.address : curr.dna.end.address]
    curr.tet.freqs <- seqinr::count(curr.gene.seq, alphabet = , wordsize = 4) # calculate the frequencies
    
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
