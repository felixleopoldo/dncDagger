source("R/scoring.R")
source("R/opruner.r") 

#filename <- "data/p20n300gaussdata.csv"
filename <- "data/asiadata.csv"
data <- read.csv(filename, check.names = FALSE)

scoretype <- "bge"
bgepar <- list(am=1, aw=NULL)

if (scoretype =="bge") {
    bidag_scores <- BiDAG::scoreparameters(scoretype = scoretype, data, bgepar = bgepar)
} else if (scoretype == "bde") {
    bidag_scores <- BiDAG::scoreparameters(scoretype = scoretype, data[-1, ], bdepar = bdepar)
}
set.seed(1)
print("getting cpp friendoy bidag scores")
cpp_friendly_scores <- get_plus1_score_essentials_for_cpp(bidag_scores, plus1it=2, iterations=NULL) # from iterativeMCMC


#initial_suborder <- seq(1, floor(ncol(data)/2)) - 1
#initial_suborder <- list(3,1,3)
initial_suborder <- list()

print("initial suborder:")
print(initial_suborder)

opr <- optimal_order(cpp_friendly_scores, initial_suborder)
adjmat <- optimal_dag(bidag_scores, cpp_friendly_scores$space, opr$order)
colnames(adjmat) <- colnames(data)
print(opr)