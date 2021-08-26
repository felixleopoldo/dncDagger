library("Rcpp")
library("Jmisc")
sourceAll(path = "R")

filename <- "jackdata.csv"
data <- read.csv("jackdata.csv", check.names = FALSE)

myscore <- scoreparameters(scoretype = "bge", data, bdecatpar = list(am = 0.1))
startorder <- seq(dim(data)[2])
startspace <- definestartspace(alpha = NULL, myscore, cpdag = TRUE, algo = "pc")
scoretable <- getScoreTable(myscore, scoreout = TRUE, MAP = FALSE, startspace = startspace)

set.seed(1)
omcmcres <- orderMCMC(myscore, scoreout = TRUE, chainout = TRUE, plus1 = TRUE, MAP = FALSE,
                      startorder = startorder, scoretable = scoretable,
                      startspace = startspace, iterations = 10000000)

saveRDS(omcmcres, "omcmc_output.rds")
