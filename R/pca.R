# Gives PCA of each Rdata file produced from calc_tetra_nuc

pca <- function(signalP.freqs, signalP.out.fn, save.file = TRUE, discard.data = FALSE) {
  browser()
  
  # Load signalP.freqs into global environment
  load(signalP.freqs)
  
  # Need to combine all of these into one matrix
  # First make an empty matrix
  tet.freqs <- matrix(data = NA_integer_, nrow = nrow(signalP_with_freqs), ncol = 256)
  
  # Use a loop to convert each tetranucleotide data table into a vector and put it into the matrix we'll use for PCA
  system.time({
    for(i in 1:nrow(signalP_with_freqs)) {
      tet.freqs[i, ] <- as.vector(signalP_with_freqs[i, "freqs"][[1]][[1]])
    }
  })
  
  # Now pull out the 'dimension names' for the tetranucleotide frequency table to get the tetranucleotides in question
  freq.names <- attr(signalP_with_freqs[1, "freqs"][[1]][[1]], "dimnames")[[1]]
  colnames(tet.freqs) <- freq.names
  
  # Pull out the 'prediction' column of the data frame - this is the 'labels' we'll use to tell rda() which vectors are exported
  prediction <- signalP_with_freqs$prediction
  
  # Standardize the matrix
  tet.freqs.std <- vegan::decostand(tet.freqs, method = "total", MARGIN = 1)
  
  # Use vegan::rda() to perform PCA and export data
  system.time({
    tet.freqs.std.rda <- vegan::rda(X=tet.freqs.std)
  })
  
  ## Graph the intracellular against extracellular enzymes with PCA values
  
  # Create data frame of PCA values and add prediction from SignalP output to compare intracellular and extracellular enzymes
  curr_PCA <- as.data.frame(tet.freqs.std.rda$CA$u) %>%
    mutate(exported = signalP_with_freqs$prediction)
  
  # Calcualte variance explained
  eigenvalues <- tet.freqs.std.rda$CA$eig
  var.explained <- eigenvalues / sum(eigenvalues)
  
  # Get organism name to use as ggplot title
  get_gg_title(signalP.out.fn)
  
  # Create ggplot 
  p <- ggplot(curr_PCA, aes(x=PC1, y=PC2, colour = exported)) + 
    geom_point(alpha = 0.3)   # make the points a bit transparent 
  
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
  
  
  
  # Do the permanova
  #raw.dist <- vegdist(tet.freqs.std, method = "euclidian")
  #adonis_obj <- adonis(raw.dist ~ prediction)
  
  
}
