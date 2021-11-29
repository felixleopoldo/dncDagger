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

seed <- 2
set.seed(seed)

reps <- 100
ns <- seq(15, 25)
ds <- c(1, 1.5, 2)
lb <- 0.25
ub <- 1
N <- 300

timing <- data.frame(matrix(ncol = 5, nrow = 0))
x <- c("n", "d", "rep", "totaltime", "max_particles")
colnames(timing) <- x

for (n in ns){
    for (d in ds){
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
            #filename <- paste("data/n=",n,"d=",d,"i=",i,".csv", sep="")
            filename <- paste("data/simdata.csv", sep="")
            write.table(data, file = filename, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = ",")
            start <- proc.time()[1]
            ret <- get_scores(filename)
            res <- r_sequential_opt(ret)
            totaltime <- proc.time()[1] - start            
            df <- data.frame(n=c(n),d=c(d),rep=c(i), totaltime=c(as.numeric(totaltime)), max_particles=c(res$max_n_particles))
            timing <- rbind(timing, df)
            #write.csv(timing, file = "data/timesparticles.csv", row.names = FALSE)            
        }
    }
}
print(timing)


ggplot(timing, aes(x=as.factor(n), y=totaltime, col=as.factor(d))) + geom_boxplot()
#write.table(timing, file = "timings.csv", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = ",")

