# Extract the organism name from the current fna file
get_gg_title <- function(signalP.out.fn) {
  
  # Parse the signalP.out.neg file
  #colnames <- c("name", "Cmax", "Cmax.pos", "Ymax", "Ymax.pos", "Smax", "Smax.pos", "Smean", "D", "prediction", "Dmaxcut", "networks.used")
  #signalP_results <- read_fwf(file = signalP.out.fn, col_positions = fwf_empty(signalP.out.fn, skip = 1, col_names = colnames), skip = 2)
  
  # Extract the root of the signalP.out.fn file name to get the fna path
  fsa.path <- stringr::str_extract(signalP.out.fn, ".+(?=(\\.out\\.neg)$)") # directory path to fsa file using stringr
  fna.path <- stringr::str_replace(fsa.path, "(aa\\.fsa)$", "fna") # directory path to fna file using stringr
  
  # Read the fna file that has the organism name under attr("Annot")
  curr_fna <- seqinr::read.fasta(fna.path) # curr_fasta means "current fasta" because we'll call this function on a whole list of fasta files
  
  # Extract the organism's name from the attr of the fna file to be used in the PCA plot
  ggplot.name <<- stringr::str_extract(attr(curr_fna[[1]],"Annot"),"(?<=\\s).{1,1000}(?=,\\s)")
 

  
  
  
  
  
  
  
  
  
  
  
  
  
}