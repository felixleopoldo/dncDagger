rm(list = ls())
library("Rcpp")
library("Jmisc")
sourceAll(path = "R")
sourceCpp("src/cppfns.cpp")


get_scores <- function(filename){

  # Discrete data
  # filename <- "data/myasiandata.csv"
  # data <- read.csv(filename, check.names = FALSE)[-1,]
  # myscore <- scoreparameters(scoretype = "bde", data, bdepar = list(chi = 0.5, edgepf = 2))

  data <- read.csv(filename, check.names = FALSE)
  myscore <- scoreparameters(scoretype = "bge", data, bgepar = list(am = 0.1))

  MAP=TRUE
  startorder <- seq(dim(data)[2])
  startspace <- definestartspace(alpha = NULL, myscore, cpdag = TRUE, algo = "pc")
  scoretable <- getScoreTable(myscore, scoreout = TRUE, MAP = MAP, startspace = startspace)

  res <- orderMCMCFelix(myscore, scoreout = TRUE, MAP = MAP, startorder = startorder, scoretable = scoretable, startspace = startspace)

  ret <- list()
  ret$parenttable <- lapply(res$ptab$parenttable, function(a) {
    df <- data.frame(a)
    df[is.na(df)] <- 0
    df <- df - 1
    m <- as.matrix(df)
    m = as.matrix(df, row.names = 0, col.names = 0)
    rownames(m) <- colnames(m) <- NULL
    return(m)
  })

  ret$aliases <- lapply(res$ptab$aliases, function(a) a - 1)
  ret$numparents <- res$ptab$numparents
  ret$rowmaps_backwards <- lapply(res$rowmaps, function(a) a$backwards - 1)
  ret$plus1listsparents <- lapply(res$plus1lists$parents, function(a) a - 1)
  ret$scoretable <- scoretable$table
  ret$bannedscore <- res$bannedscore
  ret$MAP <- MAP
  return(ret)
}