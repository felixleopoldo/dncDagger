library("Rcpp")
library("Jmisc")
library("pcalg")
library("ggplot2")
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

reps <- 20
ns <- seq(15, 25)
ds <- c(1, 1.5, 2)
lb <- 0.25
ub <- 1
N <- 300

timing <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("N","lb","ub","n", "d", "seed", "totaltime", "max_particles", "tot_particles")

colnames(timing) <- x

for (n in ns){
    for (d in ds){
        for(i in seq(reps)){
            set.seed(i)
            adjmat <- 1 * (as(pcalg::randDAG(n, d = d, method = "er"), "matrix") != 0)
            weight_mat <- adjmat
            colnames(weight_mat) <- seq(n)
            n_edges <- sum(adjmat)
            set.seed(i)
            weight_mat[which(weight_mat == 1)] <- wFUN(n_edges, lb = lb, ub = ub)
            trueDAGedges <- weight_mat
            print(paste("n: ",n, "d: ", d, "seed: ",i))
            set.seed(i)
            data <- rmvDAG(trueDAGedges, N)
            colnames(data) <- seq(n)
            #filename <- paste("data/n=",n,"d=",d,"i=",i,",seed=",seed,".csv", sep="")
            filename <- paste("data/simdata.csv", sep="")
            write.table(data, file = filename, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = ",")
            results_filename <- paste("results/n=",n,"_d=",d,"_seed=",i,"_lb=",lb,"_ub=",ub,"_N=",N,".csv", sep="")
            if(file.exists(results_filename)) {
              print(paste(results_filename,"already exists"))
            } else {                  
              start <- proc.time()[1]
              ret <- get_scores(filename)
              res <- r_sequential_opt(ret)
              totaltime <- proc.time()[1] - start
              df <- data.frame(N=c(N), ub=c(ub),lb=c(lb),n=c(n),d=c(d),seed=c(i), totaltime=c(as.numeric(totaltime)), max_particles=c(res$max_n_particles), tot_particles=c(res$tot_n_particles)) 
              write.csv(df, file = results_filename, row.names = FALSE)     
            }
            df <- read.csv(results_filename)
            timing <- rbind(timing, df)
            write.csv(timing, file = "results/timesparticles.csv", row.names = FALSE)
        }
    }
}
print(timing)
