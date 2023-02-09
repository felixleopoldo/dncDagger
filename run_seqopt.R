library("Rcpp")
library("Jmisc")
library("pcalg")
library("ggplot2")
library("testit")
library(argparser)
#library("foreach")
#library("doParallel")
# sourceAll(path = "R")
source("helper_functions.R")

# This function generates Gaussian data from a DAG
# following the topological order.

#setup parallel backend to use many processors
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)

# print(cores)
# print(cl)
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

Sys.setenv("PKG_CXXFLAGS" = "-Wall -pipe -Wno-unused -pedantic -Wall -L /usr/lib/R/lib -l R  -pthread -lpcre -llzma -lbz2 -lz -lrt -ldl -lm -L /usr/local/lib/R/site-library/RInside/lib/ -l RInside -Wl,-rpath,/usr/local/lib/R/site-library/RInside/lib  -llapack -lblas  -I /usr/local/lib/R/site-library/RInside/include/ -I /usr/local/lib/R/site-library/Rcpp/include/ -I /usr/share/R/include/ -I thread-pool-master/ -std=c++17 -O3")

sourceCpp("seq_opt.cpp", verbose = TRUE)

p <- arg_parser("Order pruning")

p <- add_argument(p, "--output_dir", help = "output dir", default = "results")
p <- add_argument(p, "--filename", help = "Output filename")
p <- add_argument(p, "--seeds_from", help = "Seeds from")
p <- add_argument(p, "--seeds_to", help = "Seeds to")
argv <- parse_args(p)
print(as.integer(argv$seeds_from))

reps <- seq(as.integer(argv$seeds_from), as.integer(argv$seeds_to))

ns <- seq(15, 26)
ds <- seq(0, 2, 0.1) #c(0, 1.5, 2)
#ds <- c(0.9)
lb <- 0.25
ub <- 1
N <- 300

timing <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("N", "lb", "ub", "n", "d", "seed", "totaltime", "max_particles", "tot_particles")

skipseeds <- c(212, 314)
#colnames(timing) <- x
#.GlobalEnv$gaussCItest <- gaussCItest # makes gaussCItest global so that it can be reached inside foreach

results <- list.files("results")

dir.create("results")
dir.create("data")
for (n in ns) {
  print(paste("n:", n))
  for (d in ds) {
    print(paste("d:", d))
    for (i in reps) {
      if(i %in% skipseeds){
        next
      }
    #foreach(i=seq(201,reps),.packages=c('pcalg', 'Rcpp')) %dopar% {
     

      results_filename <- paste("results/n=", n, "_d=", d, "_seed=", i, "_lb=", lb, "_ub=", ub, "_N=", N, ".csv", sep = "")
      #if (file.exists(results_filename)) {
      if (basename(results_filename) %in% results) {
         #print(paste(results_filename,"already exists"))
      } else {
        print(paste("n: ", n, "d: ", d, "seed: ", i))

        set.seed(i)
        adjmat <- 1 * (as(pcalg::randDAG(n, d = d, method = "er"), "matrix") != 0)
        weight_mat <- adjmat
        colnames(weight_mat) <- seq(n)
        n_edges <- sum(adjmat)
        set.seed(i)
        weight_mat[which(weight_mat == 1)] <- wFUN(n_edges, lb = lb, ub = ub)
        trueDAGedges <- weight_mat
        #print(paste("n: ", n, "d: ", d, "seed: ", i))
        set.seed(i)
        data <- rmvDAG(trueDAGedges, N)
        #print(data)
        colnames(data) <- seq(n)
        filename <- paste("data/n=",n,"d=",d,"_seed=", i, "_lb=", lb, "_ub=", ub, "_N=", N, ".csv", sep="")
        #filename <- paste("data/simdata.csv", sep = "")
        write.table(data, file = filename, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = ",")

        ret <- get_scores(filename)
        start <- proc.time()[1]
        res <- r_sequential_opt(ret)
        totaltime <- proc.time()[1] - start

        # We run the itsearch here too for comparison of the scores.
        # itsearch is not guraranteed to be optimal but foor these small scale
        # graphs it seems to usually be.
        set.seed(1) # This is just for iterative MCMC and will be overwritten
        data <- read.csv(filename, check.names = FALSE)
        myscore <- scoreparameters(scoretype = "bge", data, bgepar = list(am = 0.1))
        #itres <- iterativeMCMC(myscore, chainout = TRUE, scoreout = TRUE, MAP = TRUE) 
        #itscore <- itres$result$score
        
        
        #print("Score from iterative MCMC")
        #print(itscore)

        print("Score from order opt")
        print(res$log_score)
        #assert("same scores", abs(itscore - res$log_score) < 1e-5)
        df <- data.frame(N = c(N), ub = c(ub), lb = c(lb), n = c(n), d = c(d), seed = c(i), totaltime = c(as.numeric(totaltime)), max_particles = c(res$max_n_particles), tot_particles = c(res$tot_n_particles))
        write.csv(df, file = results_filename, row.names = FALSE)
      }
      #df <- read.csv(results_filename)
      #timing <- rbind(timing, df)
      #write.csv(timing, file = "results/timesparticles.csv", row.names = FALSE)
    }
  }
}

print("Joining csvs")
for (n in ns) {
  print(n)
  for (d in ds) {
    #for (i in seq(reps)) {
    for (i in reps) {
      if(i %in% skipseeds){
        next
      }

      results_filename <- paste("results/n=", n, "_d=", d, "_seed=", i, "_lb=", lb, "_ub=", ub, "_N=", N, ".csv", sep = "")
      df <- read.csv(results_filename)
      timing <- rbind(timing, df)
      #write.csv(timing, file = "results/timesparticles.csv", row.names = FALSE)
    }
  }
}

write.csv(timing, file = "results/timesparticles.csv", row.names = FALSE)
print(timing)
