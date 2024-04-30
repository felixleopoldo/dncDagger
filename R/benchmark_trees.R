library("Rcpp")
library("Jmisc")
library("pcalg")
library(argparser)
source("R/scoring.R")
source("R/opruner.r")
library(igraph)

# Example usage:
# $ Rscript R/benchmark_opruner.R  --filename joinedresults.csv --seeds_from 1 --seeds_to 1

p <- arg_parser("Order pruning")
p <- add_argument(p, "--output_dir", help = "output dir", default = "results_trees")
p <- add_argument(p, "--filename", help = "Output filename") # the filename of the joined results
p <- add_argument(p, "--seeds_from", help = "Seeds from")
p <- add_argument(p, "--seeds_to", help = "Seeds to")
argv <- parse_args(p)
print(as.integer(argv$seeds_from))

reps <- seq(as.integer(argv$seeds_from), as.integer(argv$seeds_to))

ns <- seq(20, 25) # Number of nodes 
ds <- seq(1) # graph density (avg indegree)
lb <- 0.25 # SEM parameters lower bound
ub <- 1 # SEM parameters upper bound
N <- 1000 # number of samples
scoretype = "bge"
# score parameters:
am <- 0.1
aw <- NULL
chi <- 0.5
edgepf <- 2

timing <- data.frame(matrix(ncol = 14, nrow = 0))
x <- c("alg", "N", "lb", "ub", "n", "d", "seed", "totaltime", "max_particles", "tot_particles")

skipseeds <- c(333) # some problem

is_matching <- function(adjmat){
  n <- nrow(adjmat)
  for(i in 1:n){
    if(sum(adjmat[i,]) + sum(adjmat[,i]) != 1){
      return(FALSE)
    }
  }
  # also check if mathcing then n is odd

  return(TRUE)
}

results <- list.files(argv$output_dir)

dir.create(argv$output_dir)

for (n in ns) {
    print(paste("n:", n))
    for (d in ds) {
      print(paste("d:", d))
      for (i in reps) {
        if(i %in% skipseeds){
          next
        }
        datastr <- paste("n=",n,"d=",d,"_seed=", i, "_lb=", lb, "_ub=", ub, "_N=", N, sep="")
        if(scoretype == "bge"){
          name <- paste(datastr, "_scoretype=", scoretype, "_am=", am, "_aw=", format(aw), ".csv", sep = "")
        } else if (scoretype == "bde") {
          name <- paste(datastr, "_scoretype=", scoretype, "_chi=", chi, "_edgepf=", edgepf, ".csv", sep = "")
        }
        results_filename <- paste(argv$output_dir,"/", name , sep = "")

        if (basename(results_filename) %in% results) {
          
        } else {
          print(paste("n: ", n, "d: ", d, "seed: ", i))
          set.seed(i)          
          #adjmat <- as.matrix(igraph::as_adj(igraph::sample_tree(n, directed=TRUE)))
          adjmat <- sample_matching(n)
          #print(adjmat)
          
          weight_mat <- adjmat
          colnames(weight_mat) <- seq(n)
          n_edges <- sum(adjmat)
          set.seed(i)
          weight_mat[which(weight_mat == 1)] <- wFUN(n_edges, lb = lb, ub = ub)
          trueDAGedges <- weight_mat
          #print(paste("n: ", n, "d: ", d, "seed: ", i))
          set.seed(i)
          data <- rmvDAG(trueDAGedges, N)
          colnames(data) <- seq(n)
          filename <- paste("data/tree_", datastr, ".csv", sep="")
          print(filename)
          write.table(data, file = filename, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = ",")
          set.seed(1)
          ret <- get_scores(filename, scoretype=scoretype, 
                            bgepar=list(am=am, aw=aw), 
                            bdepar=list(chi=chi, edgepf=edgepf), # one of these should be ignored
                            plus1it=2) 
          
          print("Running optimal order pruning")         
          start <- proc.time()[1] 
          res <- optimal_order(ret, c())
          op_totaltime <- as.numeric(proc.time()[1] - start)
          print("Total time orderpruner")
          print(op_totaltime)
          print(res)
        
          dag <- r_order_to_dag(ret, res$order)

          if (is_matching(dag)) {
             print("Estimated matching")
          } else {
            print("No matching")
          }

        #   igraph::is_tree(g, mode="total")

          
          print("Running D&C")         
          start <- proc.time()[1] 
          res_dnc <- r_dnc(ret)
          dnc_totaltime <- as.numeric(proc.time()[1] - start)
          print("Total time D&C")
          print(dnc_totaltime)
          print(res_dnc)

          # We run the itsearch here too for comparison of the scores.
          # itsearch is not guraranteed to be optimal but foor these small scale
          # graphs it seems to usually be.
          #set.seed(1) # This is just for iterative MCMC and will be overwritten
            #   data <- read.csv(filename, check.names = FALSE)
            #   myscore <- scoreparameters(scoretype = scoretype, data, bgepar = list(am = am, aw=aw))
            #   itres <- iterativeMCMC(myscore, chainout = TRUE, scoreout = TRUE, MAP = TRUE, plus1it = NULL, iterations = NULL) 
            #   itscore <- itres$result$score
            #   print("Score from iterative MCMC")
            #   print(itscore)

          #print("Score from order pruner")
          print("Score from D&C order pruner")
          print(res$log_score)
          df_dnc <- data.frame(alg=c("dnc"), N = c(N), ub = c(ub), lb = c(lb), n = c(n), d = c(d), seed = c(i), 
                          totaltime = c(dnc_totaltime), 
                          max_particles = c(res_dnc$max_n_particles), tot_particles = c(res_dnc$tot_n_particles), 
                          scoretype=c(scoretype), am=c(as.numeric(am)), aw=c(format(aw)), chi=c(chi), 
                          edgepf=c(edgepf))

          df_op <- data.frame(alg=c("orderpruner"),N = c(N), ub = c(ub), lb = c(lb), n = c(n), d = c(d), seed = c(i), 
                          totaltime = c(op_totaltime), 
                          max_particles = c(res$max_n_particles), tot_particles = c(res$tot_n_particles), 
                          scoretype=c(scoretype), am=c(as.numeric(am)), aw=c(format(aw)), chi=c(chi), 
                          edgepf=c(edgepf))
          
          

            df <- rbind(df_dnc, df_op)
          write.csv(df, file = results_filename, row.names = FALSE)

        }
      }
    }
}


print("Joining csvs")
for (n in ns) {
    print(n)
    for (d in ds) {
      for (i in reps) {
        if(i %in% skipseeds){
          next
        }

        datastr <- paste("n=",n,"d=",d,"_seed=", i, "_lb=", lb, "_ub=", ub, "_N=", N, sep="")
        if(scoretype == "bge"){
          name <- paste(datastr, "_scoretype=", scoretype, "_am=", format(am), "_aw=", format(aw), ".csv", sep = "")
        } else if (scoretype == "bde") {
          name <- paste(datastr, "_scoretype=", scoretype, "_chi=", format(chi), "_edgepf=", format(edgepf), ".csv", sep = "")
        }
        results_filename <- paste(argv$output_dir,"/", name, sep = "")
        df <- read.csv(results_filename)
        timing <- rbind(timing, df)
      }
    }
}

write.csv(timing, file = argv$filename, row.names = FALSE)
print(timing)
