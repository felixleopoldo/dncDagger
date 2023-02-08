rm(list = ls())
library("Rcpp")
library("Jmisc")
sourceAll(path = "R")
sourceCpp("src/cppfns.cpp")

get_scores <- function(filename) {
  # Discrete data

  set.seed(1)
  #data <- read.delim(filename, sep = ",", as.is = FALSE)
  #data <- dplyr::mutate_if(data, is.factor, ~ as.numeric(.x)) - 1
  
  #data <- read.csv(filename, check.names = FALSE)[-1,]
  #data <- data[1:200,]
  
  #colnames(data) <- seq(0, ncol(data) - 1)
  #colnames(data) <- seq(ncol(data))-1

  #data <- data[,c(1:5)]
  #myscore <- scoreparameters(scoretype = "bde", data, bdepar = list(chi = 1, edgepf = 1))

   data <- read.csv(filename, check.names = FALSE)
   myscore <- scoreparameters(scoretype = "bge", data, bgepar = list(am = 0.1))

  MAP <- TRUE

  #print("iterativeMCMC scoretable for node 1 (0)")
  res <- iterativeMCMC(myscore, chainout = TRUE, scoreout = TRUE, MAP = MAP) # , startspace = startspace)

  # this changes the score tables for each plus1 iteration.
  print("iterativeMCMC max score")
  print(res$result$maxorder)
  print(res$result$score)
  #print(res$result$plus1lists)
  #print(res$result$scoretable$tables[[1]])
  #print(res$result$scoretable$adjacency)


  # print("Testing scores")
  # print(myscore$nsmall)

  # order <- as.integer(res$result$maxorder)+1
  # tmp <-list()
  # tmp$maxmatrix <- res$bannedscore
  # scores <- orderscorePlus1max(
  #                           myscore$n, # number of nodes
  #                           order, # score nodes
  #                           seq(myscore$nsmall), # score positions
  #                           res$ptab$parenttable, # parenttable
  #                           res$ptab$aliases, # alisases
  #                           res$ptab$numparents, #numparents
  #                           res$plus1lists, # OBS! Tho order is changedhere compared to orderscorePlus1
  #                           res$rowmaps, # OBS! Tho order is changedhere compared to orderscorePlus1
  #                           res$result$scoretable$table,
  #                           tmp,
  #                           order)
  # print("check node scores")
  # print(order-1)
  # print(scores$totscores)
  # print(sum(scores$totscores))
  # #assert(sum(scores$totscores) == res$result$score)

  ret <- list()
  ret$parenttable <- lapply(res$ptab$parenttable, function(a) {
    df <- data.frame(a)
    df[is.na(df)] <- 0
    df <- df - 1
    m <- as.matrix(df)
    m <- as.matrix(df, row.names = 0, col.names = 0)
    rownames(m) <- colnames(m) <- NULL
    return(m)
  })

  ret$aliases <- lapply(res$ptab$aliases, function(a) a - 1)
  ret$numparents <- res$ptab$numparents
  ret$rowmaps_backwards <- lapply(res$rowmaps, function(a) a$backwards - 1)
  ret$plus1listsparents <- lapply(res$plus1lists$parents, function(a) a - 1)
  ret$scoretable <- res$result$scoretable$table
  ret$bannedscore <- res$bannedscore
  ret$MAP <- MAP

  # scoretable is a list of lists of scores corresponing to combations of the 
  # possible parents and which plus1 parent that is included.
  # The first list has the scores of combinations where no 
  # plus1 parent is included.
  # the others include one plus1 parent each in additon to the 
  # combination of possible parents. The nunbering of the list is done
  # in some special way, that Jack knows.
  # eg.
  # tables[[1]] has the scores of parent cobinations with out any pllus1 parent
  # tables[[j]], j>1 has the scores of parent combinations with
  # a plus1 parent that Jack knows.

  # banned scores has the maximum score of each parent configureation where 
  # some nodes are banned. The plus1 numbering is also structured according 
  # to some rule that i dont know.
  # it is not the inverse of scoretables since we use max (or sum).


  return(ret)
}
