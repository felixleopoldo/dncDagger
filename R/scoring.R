rm(list = ls())
library("Rcpp")
library("Jmisc")
wd <- getwd()
setwd(paste(wd,"BiDAG", sep="/"))
sourceAll(path = "R", echo=FALSE, verbose=FALSE)
sourceCpp("src/cppfns.cpp")
setwd(wd)

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
      myscore <- scoreparameters(scoretype = scoretype, data[-1, ], bdepar = bdepar)
  }

  ret <- get_plus1_score_essentials_for_cpp(myscore, seed=seed)

  return(ret)
}

get_plus1_score_essentials_for_cpp <- function(myscore, plus1it=NULL, seed=1) {
  set.seed(seed)
  MAP <- TRUE

  res <- iterativeMCMC(myscore, chainout = TRUE, scoreout = TRUE, MAP = MAP, plus1it=plus1it, verbose=TRUE) #this is bidag version 2.0.0

  print("bidag score")
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
  ret$space <- res$result$endspace

  return(ret)
}