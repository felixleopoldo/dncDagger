library("Rcpp")
# library("Jmisc")
# library("pcalg")
# library("ggplot2")
# library("testit")
library(argparser)
# source("R/scoring.R")
library("BiDAG")
Sys.setenv("PKG_CXXFLAGS" = "-Wall -pipe -Wno-unused -pedantic -Wall -L /usr/lib/R/lib -l R -L /usr/local/lib/R/site-library/RInside/lib/ -l RInside -Wl,-rpath,/usr/local/lib/R/site-library/RInside/lib  -I /usr/local/lib/R/site-library/RInside/include/ -I /usr/local/lib/R/site-library/Rcpp/include/ -I /usr/share/R/include/ -std=c++17 -O3")
sourceCpp("include/opruner_right.cpp",  verbose = TRUE)

optimal_order <- function(cpp_friendly_scores){
  return(r_opruner_right(cpp_friendly_scores))
}

optimal_dag <- function(bidag_scores, space, order) {
  print(space)
  omcmc = BiDAG::orderMCMC(bidag_scores, plus1=TRUE, MAP=TRUE, iterations=1, 
                           stepsave=1, startorder=order + 1, startspace=space)
  print("DAG score")
  print(omcmc$score)
  return(omcmc$DAG)
}
