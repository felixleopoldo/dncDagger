rm(list = ls())
library("Rcpp")
library("Jmisc")
sourceAll(path = "R")
sourceCpp("src/cppfns.cpp")


get_scores <- function(filename) {
  # Discrete data

  set.seed(1)
  data2 <- read.delim(filename, sep = ",", as.is = FALSE)
  data <- dplyr::mutate_if(data2, is.factor, ~ as.numeric(.x)) - 1
  colnames(data) <- seq(0, ncol(data) - 1)

  myscore <- scoreparameters(scoretype = "bde", data, bdepar = list(chi = 1, edgepf = 1))

  # data <- read.csv(filename, check.names = FALSE)
  # myscore <- scoreparameters(scoretype = "bge", data, bgepar = list(am = 0.1))

  MAP <- TRUE

  print("iterativeMCMC scoretable for node 1 (0)")
  res <- iterativeMCMC(myscore, chainout = TRUE, plus1it = 4, scoreout = TRUE, MAP = MAP) # , startspace = startspace)

  # this changes the score tables for each plus1 iteration.
  print("iterativeMCMC max score")
  print(res$result$maxorder)
  print(res$result$score)
  # print(res$result$scoretable)
  print(res)

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

  return(ret)
}
