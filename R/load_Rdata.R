# Load .Rdata files

Rdata.fns <- paste0("data/trial/", dir("data/trial", "\\.Rdata$"))

#pca(signalP_with_freqs)
 
#load(Rdata.fns[1])
# x <- signalP_with_freqs

files <- Rdata.fns
  names <- substr(files,1,100)

for(i in c(1:length(names))){
  name <- names[i]    
  assign (name, load(files[i]))
        
}



# sapply(Rdata.fns, load)
