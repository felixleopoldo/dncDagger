# Export score tables to GOBNILP format


number2binary = function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  if (missing(noBits)) {
    return(binary_vector)
  } else {
    binary_vector[-(1:(length(binary_vector) - noBits))]
  }
}

#Bostonpar<-scoreparameters("bge",Boston)
#itfit<-iterativeMCMC(Bostonpar, chainout=TRUE, scoreout=TRUE)

scorefile <- "scorefile.scores"

tables <- itfit$scoretable$tables
tables[[1]]

p <- length(tables)
write(p, file = scorefile, append = FALSE)

labels <- colnames(itfit$scoretable$adjacency)
adjmat <- itfit$scoretable$adjacency

# For each node
for (i in seq(p)) {
  n_plus1 <- length(tables[[i]])
  nscores <- n_plus1 * length(tables[[i]][[1]]) # the first plus1 table

  print(labels[i])  


  write(paste(labels[i], nscores), file = scorefile, append = TRUE)
  
  potparents <- adjmat[i, adjmat[i, ] == 1]
    print(potparents)
    
  # for each plus1 table. TODO: Remove the current node
  for (j in seq(n_plus1)) {
    parent_scores <- tables[[i]][[j]]
    n_parsets <- length(parent_scores)
    plus1var <- j # TODO: find the correct one.
    
    for (k in seq(n_parsets)) {
      #print(n_parents)
      bitvec <- number2binary(k - 1, length(potparents))
      str <- paste(tables[[i]][[j]][k], sum(bitvec) + 1, plus1var)

      parents <- names(potparents)[as.logical(bitvec)]
      parentstr <- paste(parents, collapse = " ")
      str <- paste(str, parentstr)
      write(str, file = scorefile, append = TRUE)
    }
  }
}
