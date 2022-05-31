# Export score tables to GOBNILP format


number2binary = function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  if (missing(noBits)) {
    return(binary_vector)
  } else {
    binary_vector[-(1:(length(binary_vector) - noBits))]
  }
}

Bostonpar<-scoreparameters("bge",Boston)
itfit<-iterativeMCMC(Bostonpar, chainout=TRUE, scoreout=TRUE)

scorefile <- "scorefile.scores"

tables <- itfit$scoretable$tables
p <- length(tables)


write(p, file = scorefile, append = FALSE)

labels <- colnames(itfit$scoretable$adjacency)
adjmat <- itfit$scoretable$adjacency

# For each node
for (i in seq(p)) {
  n_plus1 <- length(tables[[i]])
  nscores <- n_plus1 * length(tables[[i]][[1]]) # the first plus1 table
  write(paste(labels[i], nscores), file = scorefile, append = TRUE)
  potparents <- adjmat[i, adjmat[i, ] == 1]
  
  # going through the plus1 table
  # First should always be included since it has no plus1 parents
  parent_scores <- tables[[i]][[1]]
  n_parsets <- length(parent_scores)
  for (k in seq(n_parsets)) {
    bitvec <- rev(number2binary(k - 1, length(potparents)))
    str <- paste(tables[[i]][[1]][k], sum(bitvec))
    
    parents <- names(potparents)[as.logical(bitvec)]
    parentstr <- paste(parents, collapse = " ")
    str <- paste(str, parentstr)
    write(str, file = scorefile, append = TRUE)
  }
  
  j = 1
  # for each plus1 table.
  for (node in seq(p)) {
    # Fix index so it matches the tables.
    # Skip the existing parents and it self
    if(adjmat[i, node] == 1 || i ==node){
      next
    } 

    parent_scores <- tables[[i]][[j]]
    n_parsets <- length(parent_scores)
    plus1var <- labels[node] 
    
    for (k in seq(n_parsets)) {
      bitvec <- rev(number2binary(k - 1, length(potparents)))
      str <- paste(tables[[i]][[j]][k], sum(bitvec) + 1, plus1var)

      parents <- names(potparents)[as.logical(bitvec)]
      parentstr <- paste(parents, collapse = " ")
      str <- paste(str, parentstr)
      write(str, file = scorefile, append = TRUE)
    }
    j <- j+1
  }
  
}
