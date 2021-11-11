library("Rcpp")
library("Jmisc")
library("pcalg")
#sourceAll(path = "R")
source("helper_functions.R")
# This function generates Gaussian data from a DAG
# following the topological order.

rmvDAG <- function(trueDAGedges, N) {
  trueDAG <- 1 * (trueDAGedges != 0) # the edge presence in the DAG
  n <- ncol(trueDAG) # number of variables
  data <- matrix(0, nrow = N, ncol = n) # to store the simulated data
  top_order <- rev(BiDAG:::DAGtopartition(n, trueDAG)$permy) # go down order
  for (jj in top_order) {
    parents <- which(trueDAG[, jj] == 1) # find parents
    lp <- length(parents) # number of parents
    if (lp == 0) {
      # no parents
      data[, jj] <- 0
    } else if (lp == 1) {
      # one parent
      data[, jj] <- data[, parents] * trueDAGedges[parents, jj]
    } else {
      # more than one parent
      data[, jj] <- colSums(t(data[, parents]) * trueDAGedges[parents, jj])
    }
    # add random noise
    data[, jj] <- data[, jj] + rnorm(N)
  }
  data
}

wFUN <- function(m, lb, ub) {
  # This function gives edges weights between the bounds
  # with both positive and negative signs
  runif(m, lb, ub) * sample(c(-1, 1), m, replace = TRUE)
}

Sys.setenv("PKG_CXXFLAGS"="-Wall -pipe -Wno-unused -pedantic -Wall -L /usr/lib/R/lib -l R  -pthread -lpcre -llzma -lbz2 -lz -lrt -ldl -lm -L /usr/local/lib/R/site-library/RInside/lib/ -l RInside -Wl,-rpath,/usr/local/lib/R/site-library/RInside/lib  -llapack -lblas  -I /usr/local/lib/R/site-library/RInside/include/ -I /usr/local/lib/R/site-library/Rcpp/include/ -I /usr/share/R/include/ -I thread-pool-master/ -std=c++17 -O3")

sourceCpp("seq_opt.cpp", verbose=TRUE)

seed <- 1
set.seed(seed)
seeds <- rep(4)

reps <- 5
n <- 20
ns <- c(10, 20, 30)
d <- 2
ds <- c(2, 4, 8)
lb <- 0.1
ub <- 1
N <- 100
timing <- list()

for (n in ns){
    timing[[n]] <- list()
    for (d in ds){
        timing[[n]][d] <- list(c())
        adjmat <- 1 * (as(pcalg::randDAG(n, d = d, method = "er"), "matrix") != 0)
        weight_mat <- adjmat
        colnames(weight_mat) <- seq(n)
        n_edges <- sum(adjmat)
        weight_mat[which(weight_mat == 1)] <- wFUN(n_edges, lb = lb, ub = ub)
        trueDAGedges <- weight_mat
        for(i in seq(reps)){
            print(paste("n: ",n, "d: ", d, "i: ",i))
            data <- rmvDAG(trueDAGedges, N)
            colnames(data) <- seq(n)
            write.table(data, file = paste("data/",i,".csv", sep=""), row.names = FALSE, quote = FALSE, col.names = TRUE, sep = ",")
            start <- proc.time()[1]
            ret <- get_scores(paste("data/",i,".csv", sep=""))
            res <- r_sequential_opt(ret)
            totaltime <- proc.time()[1] - start
            timing[[n]][[d]] <- c(timing[[n]][[d]], as.numeric(totaltime))
        }
    }
}

for (n in ns){
    for (d in ds){
        print(paste("n: ",n, "d: ", d))
        print(timing[[n]][[d]])
        print(paste("median: ", median(timing[[n]][[d]]), "mean: ", mean(timing[[n]][[d]]),"std: ", sd(timing[[n]][[d]])))
    }
}

saveRDS(timings, file="timings.rds")