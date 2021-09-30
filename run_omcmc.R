#library(BiDAG)
library("gRbase")
library("Rcpp")
library("Jmisc")
sourceAll(path = "R")
sourceCpp("src/cppfns.cpp")

#adjmat <- as.matrix(read.csv(file.path("data", "adjmat.csv"), check.names=FALSE))
#order_true <- as.integer(topo_sort(adjmat))
#print(order_true)

# filename <- "data/myasiandata.csv"
# data <- read.csv(filename, check.names = FALSE)[-1,]
# myscore <- scoreparameters(scoretype = "bde", data, bdecatpar = list(chi = 0.5, edgepf = 2))

filename <- "data/p20n300gaussdata.csv"
data <- read.csv(filename, check.names = FALSE)
myscore <- scoreparameters(scoretype = "bge", data, bgepar = list(am = 0.1))

startspace <- definestartspace(alpha = NULL, myscore, cpdag = TRUE, algo = "pc")
#scoretable <- getScoreTable(myscore, scoreout = TRUE, MAP = TRUE, startspace = startspace)

set.seed(1)
# This has a different start space probably...
#startorder <- seq(dim(data)[2])
startorder <- c(19,14,18,12,17,8,16,7,9,15,10,11,5,1,6,3,13,4,2,0)+1
omcmcres <- orderMCMC(myscore, scoreout = TRUE, chainout = TRUE, plus1 = TRUE, MAP = TRUE,
                      startorder = startorder, 
#                      scoretable= scoretable,
                        startspace = startspace,
                        iterations = 1000)
# 18 19 17 15 10 14  8 11  5 12 16  7  9  1  6  3  2 13  4  0
# -6416.047
# after 10 min

oscores <- omcmcres$traceadd$orderscores
orders <-  omcmcres$traceadd$orders

print(oscores)
print(as.integer(unlist(omcmcres$maxorder))-1)
max(oscores)

# orderscorePlus1(n, scorenodes, scorepositions, parenttable, aliases, numparents,
#                 rowmaps, plus1lists, scoretable, scoresmatrices, permy)


saveRDS(omcmcres, "omcmc_output.rds")
