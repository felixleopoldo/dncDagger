
source("helper_functions.R")

sumsc <- function(scores) {
  return(log(sum(exp(scores - max(scores)))) + max(scores))
}

# filename <- "data/myasiandata.csv"
# data <- read.csv(filename, check.names = FALSE)[-1,]
# myscore <- scoreparameters(scoretype = "bde", data, bdepar = list(chi = 0.5, edgepf = 2))

#filename <- "data/avneigs8p30n300.csv"
filename <- "data/p20n300gaussdata.csv"
#filename <- "data/myvstructdata.csv"
# #filename <- "data/p50n300gaussdata.csv"
# #filename <- "data/jackdata.csv"


ret <- get_scores(filename)
print(ret$scoretable)
#print(ret$bannedscore)
saveRDS(ret, file = paste(filename, "rds", sep = "."))

#order <- startorder #sample(startorder, replace = FALSE)
#Check scores for some ordering. And compare to C++ node scores.
#order <- 1+c(19,14,18,12,17,8,16,7,9,15,10,11,5,1,6,3,13,4,2,0)
# if(MAP==FALSE){
#   scores <- orderscorePlus1(
#                             myscore$n,
#                             order,
#                             seq(myscore$nsmall),
#                             res$ptab$parenttable,
#                             res$ptab$aliases,
#                             res$ptab$numparents,
#                             res$rowmaps,
#                             res$plus1lists,
#                             scoretable$table,
#                             res$bannedscore,
#                             order)
# } else {
# tmp <-list()
# tmp$maxmatrix <- res$bannedscore
# scores <- orderscorePlus1max(
#                           myscore$n,
#                           order,
#                           seq(myscore$nsmall),
#                           res$ptab$parenttable,
#                           res$ptab$aliases,
#                           res$ptab$numparents,
#                           res$plus1lists, # OBS! Tho order is changedhere compared to orderscorePlus1
#                           res$rowmaps, # OBS! Tho order is changedhere compared to orderscorePlus1
#                           scoretable$table,
#                           tmp,
#                           order)
# }

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
