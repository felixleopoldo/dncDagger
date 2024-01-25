library("Rcpp")
# library("Jmisc")
# library("pcalg")
# library("ggplot2")
# library("testit")
library(argparser)
library("BiDAG")
Sys.setenv("PKG_CXXFLAGS" = "-Wall -pipe -Wno-unused -pedantic -Wall -L /usr/lib/R/lib -l R -L /usr/local/lib/R/site-library/RInside/lib/ -l RInside -Wl,-rpath,/usr/local/lib/R/site-library/RInside/lib  -I /usr/local/lib/R/site-library/RInside/include/ -I /usr/local/lib/R/site-library/Rcpp/include/ -I /usr/share/R/include/ -std=c++17 -O3")
sourceCpp("include/opruner_right.cpp",  verbose = TRUE)


optimal_order <- function(cpp_friendly_scores, initial_sub_order){
    ret <- r_opruner_right(cpp_friendly_scores, initial_sub_order)
    ret$suborder_cond_score <- 0
    ret$order <- ret$order + 1 # since c++ enumerate nodes from 0    
    p <- length(ret$node_scores)
    suborder <- ret$order[1:(p-length(initial_sub_order))]
    for (node in suborder) {
        ret$suborder_cond_score <- ret$suborder_cond_score + ret$node_scores[node]
    } 
    ret$suborder <- suborder
    return(ret)
}

optimal_dag <- function(bidag_scores, space, order) {    
  omcmc = BiDAG::orderMCMC(bidag_scores, plus1=TRUE, MAP=TRUE, iterations=1, 
                           stepsave=1, startorder=order, startspace=space)

  return(omcmc$DAG)
}
