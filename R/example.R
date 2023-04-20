
source("R/opruner.r") 
library("BiDAG")

filename <- "data/p20n300gaussdata.csv"
data <- read.csv(filename, check.names = FALSE)

scoretype <- "bge"
bgepar <- list(am=1, aw=NULL)

if (scoretype =="bge") {
    bidag_scores <- BiDAG::scoreparameters(scoretype = scoretype, data, bgepar = bgepar)
} else if (scoretype == "bde") {
    bidag_scores <- BiDAG::scoreparameters(scoretype = scoretype, data[-1, ], bdepar = bdepar)
}

cpp_friendly_scores <- get_plus1_score_essentials_for_cpp(bidag_scores) # from iterativeMCMC

opr <- opruner(cpp_friendly_scores)
opr

omcmc = BiDAG::orderMCMC(bidag_scores, plus1=TRUE, MAP=TRUE, iterations=1, 
                         stepsave=1, startorder=opr$order + 1)

omcmc$score