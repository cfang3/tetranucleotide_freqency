# This takes the filename of a signalP output and returns a list of tetranucleotide frequencies of the associated DNA sequences

calc_tetra_nuc <- function(signalP.out.fn, save.file = TRUE, discard.data = FALSE) {
  

  # Read the SignalP output file = signalP.out.fn
  colnames <- c("name", "Cmax", "Cmax.pos", "Ymax", "Ymax.pos", "Smax", "Smax.pos", "Smean", "D", "prediction", "Dmaxcut", "networks.used")
  signalP_results <- read_fwf(file = signalP.out.fn, col_positions = fwf_empty(signalP.out.fn, skip = 1, col_names = colnames), skip = 2)
  
  # Strip the excess information off the sequence name to pull out the root of the fasta file name to get the paths of the fna and fsa file
  fsa.path <- str_extract(signalP.out.fn, ".+(?=(\\.out\\.neg)$)") # regular expression: match one or more characters of everything before .out.neg
  fna.path <- str_replace(fsa.path, "(aa\\.fsa)$", "fna") # regular expression: take characters aa and every character following fsa to the end of the string and replace it with fna
  
  # Read two different fasta files: the fna file, which has DNA sequences with multiple genes (whole genome shotgun sequence) 
  # AND the fsa file, with amino acid sequences of single genes (which is what we fed to signalP)
  curr_fna <- seqinr::read.fasta(fna.path) # curr_fasta means "current fasta" because we'll call this function on a whole list of fasta files
  curr_fsa_file <- seqinr::read.fasta(fsa.path)  # Open the amino acid sequence file that has the start and end address of the DNA sequence 

  
  # For each sequence, calculate tetranucleotide frequencies
  freq_list <- vector("list", nrow(signalP_results)) # "pre-allocate" a list a number of elements equal to the number of rows of the SignalP output
  freq_results_names <- vector("character", nrow(signalP_results)) 
  for(i in 1:nrow(signalP_results)) {
    
    # Extract the current sequence name from the signalP output; get rid the underscore and everything after it
    curr.aa.seq.name <- as.character(signalP_results[i, "name"]) # annoyingly, if we don't wrap signalP_results[i, "name"] in as.character(), we'll get back a tibble 
    
    # Find the header of that sequence in the curr_fsa_file
    curr.fsa.header <- attr(curr_fsa_file[[curr.aa.seq.name]], "Annot")
    
    # Extract the two numbers after the two # signs; these are the first and last addresses of the gene sequence that we're looking at
    curr.dna.start.address <- as.numeric(str_extract(curr.fsa.header, "(?<=#\\s)\\d+(?=\\s#)")) # regular expression: match everything behind the first # and space and dont include it; followed by any number of digits; match everything after the next space and # and dont include it
    curr.dna.end.address <- as.numeric(str_extract(curr.fsa.header, "(?<=#\\s\\d{1,10}\\s#\\s)\\d+(?=\\s#)")) # regular expression: look behind the first #, space, number, space, #, space and dont include it; followed by any amount of digits; then match everything after space and # and dont include it
    
    # Get the DNA sequence we're dealing with & extract the gene sequence that corresponds to the AA sequence that signalP looked at
    curr.dna.seq.name <- str_extract(curr.aa.seq.name, ".+(?=(_(\\d+)$))") #regular expression: Take any number of character/numbers until there is a "_". Take the "_" and any number after it until the end of the string and leave it out.
    
    # Extract the actual DNA gene sequence from the DNA fna file, based on the name of the sequence and its start and end address
    curr.dna.seq <- curr_fna[[curr.dna.seq.name]]
    curr.gene.seq <- curr.dna.seq[curr.dna.start.address : curr.dna.end.address]  # gets the actual nuclotide sequence based on the start and end address in which the tetranucleotide frequency will be measured
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