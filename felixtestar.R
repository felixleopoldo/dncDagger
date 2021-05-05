rm(list = ls())
library("Rcpp")
library("Jmisc")
sourceAll(path = "R")
sourceCpp("src/cppfns.cpp")

sumsc <- function(scores) {
  return(log(sum(exp(scores - max(scores)))) + max(scores))
}

data <- read.csv("myasiandata.csv")[-1,]

myscore <- scoreparameters(scoretype = "bde", data, bdepar = list(chi = 0.5, edgepf = 2))

startorder <- c(1, 4, 7, 2, 6, 3, 8, 5)

startspace <- definestartspace(alpha = NULL, myscore, cpdag = TRUE, algo = "pc")
startspace[8, 3] <- 1
startspace
scoretable <- getScoreTable(myscore, scoreout = TRUE, MAP = FALSE, startspace = startspace)
scoretable

res <- orderMCMCFelix(myscore, scoreout = TRUE, MAP = FALSE, startorder = startorder, scoretable = scoretable, startspace = startspace)
res

length(res$ptab) 

set.seed(3)

order <- sample(startorder, replace = FALSE)

scores <- orderscorePlus1Felix(
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
print(order)
print(scores$totscores)
ret <- list()
# These and scoretable should be read into c++ in some way.

ret$order <- order - 1 # fix indexing 
ret$parenttable <- lapply(res$ptab$parenttable, function(a) {
                                                    df <- data.frame(a)
                                                    df[is.na(df)] <- 0
                                                    df <- df-1
                                                    m <- as.matrix(df)
                                                    m = as.matrix(df,row.names=0,col.names=0)
                                                    rownames(m) <- colnames(m) <- NULL
                                                    return(m)
                                                    })

ret$aliases <- lapply(res$ptab$aliases, function(a) a-1)
ret$numparents <- res$ptab$numparents
ret$rowmaps_backwards <- lapply(res$rowmaps, function(a) a$backwards -1)
ret$plus1listsparents <- lapply(res$plus1lists$parents, function(a) a-1) 
ret$scoretable <- scoretable$table
ret$bannedscore <- res$bannedscore 

ret
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