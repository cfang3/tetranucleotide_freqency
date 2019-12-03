# Gives PCA of each Rdata file produced from calc_tetra_nuc

pca <- function(SignalP.freqs, save.file = TRUE, discard.data = FALSE) {
  
  
  
  # Need to combine all of these into one matrix
  # First make an empty matrix
  tet.freqs <- matrix(data = NA_integer_, nrow = nrow(SignalP.freqs), ncol = 256)
  
  # Use a loop to convert each tetranucleotide data table into a vector and put it into the matrix we'll use for PCA
  system.time({
    for(i in 1:nrow(SignalP.freqs)) {
      tet.freqs[i, ] <- as.vector(SignalP.freqs[i, "freqs"][[1]][[1]])
    }
  })
 
  # Now pull out the 'dimension names' for the tetranucleotide frequency table to get the tetranucleotides in question
  freq.names <- attr(SignalP.freqs[1, "freqs"][[1]][[1]], "dimnames")[[1]]
  colnames(tet.freqs) <- freq.names
  
  
  # Pull out the 'prediction' column of the data frame - this is the 'labels' we'll use to tell rda() which vectors are exported
  prediction <- SignalP.freqs$prediction

  # Use vegan::rda() to perform PCA and export data
  # First standardize the matrix
  library(vegan)
  tet.freqs.std <- decostand(tet.freqs, method = "total", MARGIN = 1)
  
  system.time({
    tet.freqs.std.rda <- vegan::rda(X=tet.freqs.std)
  })
  
  
  
  # See what happens when we do this the ggplot2 way
  curr_PCA <- as.data.frame(tet.freqs.std.rda$CA$u) %>%
    mutate(exported = SignalP.freqs$prediction)
  
  # Calcualte variance explained
  eigenvalues <- tet.freqs.std.rda$CA$eig
  var.explained <- eigenvalues / sum(eigenvalues)

  p <- ggplot(curr_PCA, aes(x=PC1, y=PC2, colour = exported)) + 
    geom_point(alpha = 0.3)   # make the points a bit transparent 
  
  graph <- p +
    xlab(paste0("PC1, ", round(var.explained["PC1"], digits = 2)*100, "%")) +
    ylab(paste0("PC2, ", round(var.explained["PC2"], digits = 2)*100, "%")) + 
    #scale_color_brewer(type = "qual") + 
    scale_color_manual(values = c("red", "black")) + 
    ggtitle("Ureaplasma parvum serovar 3 str. ATCC 700970") +
    theme_minimal() + 
    theme(text = element_text(size = 24)) 
  print(graph)
  browser()
  #scale_colour_manual(values = c())
  ggsave("plots/PCA_plot_row_1.png", height = 8, width = 10, units = "in", dpi = 300)
  
  
  
  # Do the permanova
  #raw.dist <- vegdist(tet.freqs.std, method = "euclidian")
  #adonis_obj <- adonis(raw.dist ~ prediction)
  
  
}