library(BiDAG)
library("gRbase")
library("Rcpp")
library("Jmisc")
sourceAll(path = "R")
sourceCpp("src/cppfns.cpp")

#adjmat <- as.matrix(read.csv(file.path("data", "adjmat.csv"), check.names=FALSE))
#order_true <- as.integer(topo_sort(adjmat))
#print(order_true)

filename <- "data/myasiandata.csv"
data <- read.csv(filename, check.names = FALSE)[-1,]
myscore <- scoreparameters(scoretype = "bde", data, bdecatpar = list(chi = 0.5, edgepf = 2))

#filename <- "data/jackdata.csv"
#data <- read.csv(filename, check.names = FALSE)
#myscore <- scoreparameters(scoretype = "bge", data, bdecatpar = list(am = 0.1))

startorder <- seq(dim(data)[2])
#startspace <- definestartspace(alpha = NULL, myscore, cpdag = TRUE, algo = "pc")
#scoretable <- getScoreTable(myscore, scoreout = TRUE, MAP = FALSE, startspace = startspace)

set.seed(1)
#omcmcres <- orderMCMC(myscore, scoreout = TRUE, chainout = TRUE, plus1 = TRUE, MAP = FALSE,
#                      startorder = startorder, stepsave = 10,
#                      startspace = startspace, iterations = 10000000000000000)


omcmcres <- orderMCMC(myscore, scoreout = TRUE, chainout = TRUE, plus1 = TRUE, MAP = FALSE,
                      startorder = startorder, stepsave = 10, iterations = 10000000000000000)

saveRDS(omcmcres, "omcmc_output.rds")
