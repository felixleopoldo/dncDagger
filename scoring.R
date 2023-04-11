rm(list = ls())
library("Rcpp")
library("Jmisc")
sourceAll(path = "R", echo=FALSE, verbose=FALSE)
sourceCpp("src/cppfns.cpp")
source("export_to_gobnilp_score_tables.R")

get_scores <- function(filename,  scoretype = c("bge", "bde", "bdecat"),
                      bgepar = list(am = 1, aw = NULL),
                      bdepar = list(chi = 0.5, edgepf = 2),
                      seed=1) {
  
  set.seed(seed)

  
  
  MAP <- TRUE

  data <- read.csv(filename, check.names = FALSE)
  if (scoretype =="bge") {
      myscore <- scoreparameters(scoretype = scoretype, data, bgepar = bgepar)
  } else if (scoretype == "bde") {
      #data <- dplyr::mutate_if(data, is.factor, ~ as.numeric(.x)) - 1
      myscore <- scoreparameters(scoretype = scoretype, data[-1, ], bdepar = bdepar)
  }

  res <- iterativeMCMC(myscore, chainout = TRUE, scoreout = TRUE, MAP = MAP, verbose=FALSE) # , startspace = startspace)


  # this changes the score tables for each plus1 iteration.
  #print("iterativeMCMC max score")
  #print(res$result$maxorder)
  print(res$result$score)

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
  # tables[[1]] has the scores of parent combinations with out any plus1 parent
  # tables[[j]], j>1 has the scores of parent combinations with
  # a plus1 parent that Jack knows.

  # banned scores has the maximum score of each parent configureation where 
  # some nodes are banned. The plus1 numbering is also structured according 
  # to some rule that i dont know.
  # it is not the inverse of scoretables since we use max (or sum).


  return(ret)
}
