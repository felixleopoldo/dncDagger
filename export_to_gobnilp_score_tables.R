# Export score tables to GOBNILP format
library("BiDAG")

number2binary = function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  if (missing(noBits)) {
    return(binary_vector)
  } else {
    binary_vector[-(1:(length(binary_vector) - noBits))]
  }
}


data2 <- read.delim("asia_100.dat", sep = " ", as.is = FALSE)
data <- dplyr::mutate_if(data2, is.factor, ~ as.numeric(.x)) - 1
colnames(data) <- seq(0, ncol(data)-1)


scores <- scoreparameters("bde",data, bdepar = list(chi = 1, edgepf = 1))
itfit <- iterativeMCMC(scores, chainout=TRUE, scoreout=TRUE)


scorefile <- "scorefile.scores"

tables <- itfit$scoretable$tables
p <- length(tables)
write(p, file = scorefile, append = FALSE)

labels <- colnames(itfit$scoretable$adjacency)
adjmat <- itfit$scoretable$adjacency

conv <- function(){
  
}

# For each node
for (i in seq(p)) {
  n_plus1 <- length(tables[[i]])
  nscores <- n_plus1 * length(tables[[i]][[1]]) # the first plus1 table
  write(paste(labels[i], nscores), file = scorefile, append = TRUE)
  potparents <- colnames(adjmat)[adjmat[i, ] == 1]
  
  # going through the plus1 table
  # First should always be included since it has no plus1 parents
  parent_scores <- tables[[i]][[1]]
  n_parsets <- length(parent_scores)
  for (k in seq(n_parsets)) {
    bitvec <- rev(number2binary(k - 1, length(potparents)))
    str <- paste(tables[[i]][[1]][k], sum(bitvec))
    
    parents <- potparents[as.logical(bitvec)]
    parentstr <- paste(parents, collapse = " ")
    str <- paste(str, parentstr)
    write(str, file = scorefile, append = TRUE)
  }

  # The below is for the plus1 tables
  j = 1
  for (node in seq(p)) {
    # Fix index so it matches the tables.
    # Skip the existing parents and it self
    if(adjmat[i, node] == 1 || i ==node){
      next
    } 

    parent_scores <- tables[[i]][[j+1]] # get the plus1 score table
    
    n_parsets <- length(parent_scores)
    plus1var <- labels[node] 
    
    for (k in seq(n_parsets)) {
      # get bitvector of 
      bitvec <- rev(number2binary(k - 1, length(potparents))) 
      str <- paste(parent_scores[k], sum(bitvec) + 1, plus1var)
      parents <- potparents[as.logical(bitvec)]
      parentstr <- paste(parents, collapse = " ")
      str <- paste(str, parentstr)
      write(str, file = scorefile, append = TRUE)
    }
    j <- j+1
  }
  
}
