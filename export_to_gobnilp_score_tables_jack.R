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

parent_table <- function(n, p) {
  rbind(rep(NA, n), BiDAG:::parentlistnonempty(1:p, n))[, 1:p, drop = FALSE]
}

set.seed(1)
filename <- "data/datap25d6.csv"
data <- read.csv(filename, check.names = FALSE)[-1,]
#data <- data[1:200,]
colnames(data) <- seq(0, ncol(data)-1)


scores <- scoreparameters("bde",data, bdepar = list(chi = 1, edgepf = 1))
itfit <- iterativeMCMC(scores, chainout=TRUE, scoreout=TRUE)

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
  potparents <- colnames(adjmat)[adjmat[, i] == 1]
  
  parenttable <- parent_table(p, length(potparents))
  
  # going through the plus1 table
  # First should always be included since it has no plus1 parents
  parent_scores <- tables[[i]][[1]]
  n_parsets <- length(parent_scores)
  for (k in seq(n_parsets)) {
    parvec <- parenttable[k, ]
    parvec <- parvec[which(!is.na(parvec))]
    str <- paste(tables[[i]][[1]][k], length(parvec))
    
    parents <- potparents[parvec]
    parentstr <- paste(parents, collapse = " ")
    str <- paste(str, parentstr)
    write(str, file = scorefile, append = TRUE)
  }

  # The below is for the plus1 tables
  j = 1
  for (node in seq(p)) {
    # Fix index so it matches the tables.
    # Skip the existing parents and it self
    if(adjmat[node, i] == 1 || i ==node){
      next
    } 

    parent_scores <- tables[[i]][[j+1]] # get the plus1 score table
    
    n_parsets <- length(parent_scores)
    plus1var <- labels[node] 
    
    for (k in seq(n_parsets)) {
      parvec <- parenttable[k, ]
      parvec <- parvec[which(!is.na(parvec))]
      str <- paste(parent_scores[k], length(parvec) + 1, plus1var)
      parents <- potparents[parvec]
      parentstr <- paste(parents, collapse = " ")
      str <- paste(str, parentstr)
      write(str, file = scorefile, append = TRUE)
    }
    j <- j+1
  }
  
}


