rm(list = ls())
library("Rcpp")
library("Jmisc")
sourceAll(path = "R")
sourceCpp("src/cppfns.cpp")

sumsc <- function(scores) {
  return(log(sum(exp(scores - max(scores)))) + max(scores))
}

#filename <- "myasiandata.csv"
filename <- "jackdata.csv"
data <- read.csv("jackdata.csv", check.names =FALSE)
#data <- read.csv(filename, check.names = FALSE)[-1,]
#data <- read.csv("myhepar2data2000.csv")[-1,]
vars = names(data)
#myscore <- scoreparameters(scoretype = "bdecat", data, bdecatpar = list(chi = 0.5, edgepf = 2))
#myscore <- scoreparameters(scoretype = "bde", data, bdecatpar = list(chi = 0.5, edgepf = 2))
myscore <- scoreparameters(scoretype = "bge", data, bdecatpar = list(am = 0.1))

startorder <- seq(dim(data)[2])

startspace <- definestartspace(alpha = NULL, myscore, cpdag = TRUE, algo = "pc")

scoretable <- getScoreTable(myscore, scoreout = TRUE, MAP = FALSE, startspace = startspace)
scoretable
vars
omcmcres <- orderMCMC(myscore, scoreout = TRUE, chainout = TRUE, plus1 = TRUE, MAP = FALSE, startorder = startorder, scoretable = scoretable, startspace = startspace)
omcmcres$maxorder
#omcmcres$score
#plot(omcmcres$trace, type="l")

max(omcmcres$traceadd$orderscores)
plot(omcmcres$traceadd$orderscores)
#print(as.integer(omcmcres$maxorder) - 1)
#print(omcmcres$score)


omcmc_orderindex <- c()
for (v in omcmcres$maxorder) {
  omcmc_orderindex <- c(omcmc_orderindex, which(vars == v)[[1]] - 1)
}
omcmc_orderindex
# smcorder <- c(32, 48, 37, 50, 8, 7, 27, 55, 52, 28, 66, 53, 25, 33, 11, 56, 12, 69, 31, 13, 49, 62, 6, 38, 30, 21, 46, 34, 24, 58, 2, 67, 20, 65, 63, 60, 40, 57, 36, 61, 64, 5, 18, 44, 35, 14, 9, 1, 4, 0, 39, 29, 17, 59, 3, 22, 54, 43, 45, 68, 10, 26, 47, 42, 16, 19, 51, 15, 41, 23)+1

# vars[smcorder]




res <- orderMCMCFelix(myscore, scoreout = TRUE, MAP = FALSE, startorder = startorder, scoretable = scoretable, startspace = startspace)
res

length(res$ptab)

set.seed(3)

order <- startorder #sample(startorder, replace = FALSE)

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

ret <- list()
# These and scoretable should be read into c++ in some way.

ret$order <- order - 1 # fix indexing 
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
print(res$bannedscore)

saveRDS(ret, file = paste(filename, "rds", sep = "."))

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
