rm(list = ls())
library("Rcpp")

library("Jmisc")
sourceAll(path = "R")
sourceCpp("src/cppfns.cpp")

sumsc <- function(scores) {
  return(log(sum(exp(scores - max(scores)))) + max(scores))
}

# filename <- "data/myasiandata.csv"
# data <- read.csv(filename, check.names = FALSE)[-1,]
# myscore <- scoreparameters(scoretype = "bde", data, bdepar = list(chi = 0.5, edgepf = 2))

filename <- "data/avneigs8p30n300.csv"
#filename <- "data/p20n300gaussdata.csv"
#filename <- "data/myvstructdata.csv"
# #filename <- "data/p50n300gaussdata.csv"
# #filename <- "data/jackdata.csv"
data <- read.csv(filename, check.names = FALSE)
myscore <- scoreparameters(scoretype = "bge", data, bgepar = list(am = 0.1))

MAP=TRUE
startorder <- seq(dim(data)[2])
startspace <- definestartspace(alpha = NULL, myscore, cpdag = TRUE, algo = "pc")
scoretable <- getScoreTable(myscore, scoreout = TRUE, MAP = MAP, startspace = startspace)

res <- orderMCMCFelix(myscore, scoreout = TRUE, MAP = MAP, startorder = startorder, scoretable = scoretable, startspace = startspace)

set.seed(3)

order <- startorder #sample(startorder, replace = FALSE)

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

saveRDS(ret, file = paste(filename, "rds", sep = "."))

#Check scores for some ordering. And compare to C++ node scores.
#order <- 1+c(19,14,18,12,17,8,16,7,9,15,10,11,5,1,6,3,13,4,2,0)
if(MAP==FALSE){
  scores <- orderscorePlus1(
                            myscore$n,
                            order,
                            seq(myscore$nsmall),
                            res$ptab$parenttable,
                            res$ptab$aliases,
                            res$ptab$numparents,
                            res$rowmaps,
                            res$plus1lists,
                            scoretable$table,
                            res$bannedscore,
                            order)
} else {
tmp <-list()
tmp$maxmatrix <- res$bannedscore
scores <- orderscorePlus1max(
                          myscore$n,
                          order,
                          seq(myscore$nsmall),
                          res$ptab$parenttable,
                          res$ptab$aliases,
                          res$ptab$numparents,
                          res$plus1lists, # OBS! Tho order is changedhere compared to orderscorePlus1
                          res$rowmaps, # OBS! Tho order is changedhere compared to orderscorePlus1
                          scoretable$table,
                          tmp,
                          order)
}

#-436.407,-412.631,-377.402,-252.436,-426.504,-411.064,-306.037,-427.981,-226.323,-408.531,-331.536,-224.603,-373.228,-299.423,-234.915,-322.707,-300.573,-222.472,-251.861,-167.25,
#ret
#res <- orderMCMCFelix(myscore, scoreout = TRUE, MAP = FALSE, startorder = startorder, scoretable=NULL)

#res <- orderMCMC(myscore, scoreout = TRUE, MAP = FALSE)
#score_tables <- res$scoretable$tables
#res$scoretable$tables
#res$scoretable$therow

# order score
# log scores by node, need to add
# each node or row 

# store in log spcace but add in real space

# score tables are already summed over foppibe edge sets.

# look up roe, look ut list
# allowed 1,2 ,3, 4,8
#- 76.887573
#scores <- c(-77.13198, -78.90716, -81.00421, -80.23902, -80.30517)

# allowed and excluded are complement. lenthg -


#print(sumsc(scores))
# allowed<- 1,3,4
#scondnoce <- c(-30.74563, -26.91450, -4.35449)

#sumsc(scondnoce)

# the plus one tables ar the extra 
# what we did by hand could be done by masking.
# allowedlist has tthe pliu one nodes
# wichout the plus one, would the first table only be needed
# bottleneck might be to check the allowed columns

# try to avoid recomputation of whole oerde score.
