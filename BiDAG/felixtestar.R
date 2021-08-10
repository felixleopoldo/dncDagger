rm(list = ls())
library("Rcpp")
library("Jmisc")
sourceAll(path = "R")
sourceCpp("src/cppfns.cpp")

set.seed(1)

data <- read.csv("../../git/benchpress/resources/data/mydatasets/myasiandata.csv")[-1,]

myscore <- scoreparameters(scoretype="bde", data, bdepar = list(chi = 0.5, edgepf = 2))
res <- orderMCMCFelix(myscore, scoreout = TRUE, MAP = FALSE)
#res <- orderMCMC(myscore, scoreout = TRUE, MAP = FALSE)
#score_tables <- res$scoretable$tables
#res$scoretable$tables

# order score
# log scores by node, need to add
# each node or row 

# store in log spcace but add in real space


# score tables are already summed over foppibe edge sets.

# look up roe, look ut list
# allowed 1,2 ,3, 4,8
-76.887573 
scores <- c(-77.13198, -78.90716, -81.00421,  -80.23902, -80.30517 )

# allowed and excluded are complement. lenthg -

sumsc <- function(scores){
    return(log(sum(exp(scores - max(scores)))) + max(scores))
}

print(sumsc(scores))
# allowed<- 1,3,4
scondnoce <- c(-30.74563, -26.91450, -4.35449 )

print(sumsc(scondnoce))

# the plus one tables ar the extra 
# what we did by hand could be done by masking.
# allowedlist has tthe pliu one nodes
# wichout the lus one, it wyould be the first table only needed
# bottleneck might be to check the allowed columns

# try to avoid recomputation of whole oerde score.