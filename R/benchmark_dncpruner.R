library("Rcpp")
library("Jmisc")
library("pcalg")
library(argparser)
source("R/scoring.R")
source("R/opruner.R")

# Example usage:
# $ Rscript R/benchmark_opruner.R  --filename joined_results.csv --seeds_from 1 --seeds_to 3

p <- arg_parser("Order pruning")
p <- add_argument(p, "--output_dir", help = "output dir", default = "results")
p <- add_argument(p, "--filename", help = "Output filename", default ="joined_results.csv") # the filename of the joined results
p <- add_argument(p, "--seeds_from", help = "Seeds from")
p <- add_argument(p, "--seeds_to", help = "Seeds to")
argv <- parse_args(p)
print(as.integer(argv$seeds_from))

reps <- seq(as.integer(argv$seeds_from), as.integer(argv$seeds_to))

ns <- c(15, 16, 17) # Number of nodes 
ds <- seq(0, 2, 0.1) # graph density (avg indegree)

lb <- 0.25 # SEM parameters lower bound
ub <- 1 # SEM parameters upper bound
N <- 300 # number of samples
scoretype = "bge"
# score parameters:
am <- 0.1
aw <- NULL
chi <- 0.5
edgepf <- 2

algs <- c("dnc")

gobnilp <- TRUE

timing <- data.frame(matrix(ncol = 14, nrow = 0))
x <- c("alg", "N", "lb", "ub", "n", "d", "seed", "totaltime", "max_particles", "tot_particles")

skipseeds <- c() # some problem

results <- list.files(argv$output_dir)

dir.create(argv$output_dir)
dir.create(paste0(argv$output_dir,"/gobnilp_scores"))
#dir.create("results/gobnilp_scores")


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
        gobnilp_scores_filename  <- NULL
        if ("gobnilp" %in% algs) {
            gobnilp_scores_filename  <- paste(argv$output_dir,"/gobnilp_scores/", name , sep = "")
        } 

        if (basename(results_filename) %in% results) {
          # Do nothing
        } else {
          print(paste("n: ", n, "d: ", d, "seed: ", i))
          
          # create empty data frame to be filled up with results
            df <- data.frame(matrix(ncol = 15, nrow = 0))
            colnames(df) <- c("alg", "N", "ub", "lb", "n", "d", "seed", "totaltime", "max_particles", "tot_particles", "scoretype", "am", "aw", "chi", "edgepf")

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
          colnames(data) <- seq(n)
          filename <- paste("data/", datastr, ".csv", sep="")
          print(filename)
          write.table(data, file = filename, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = ",")
          set.seed(1)
          ret <- get_scores(filename, scoretype=scoretype, 
                            bgepar=list(am=am, aw=aw), 
                            bdepar=list(chi=chi, edgepf=edgepf), # one of these should be ignored
                            plus1it=2,
                            gobnilp_scores_filename=gobnilp_scores_filename
                            )
          
          if ("dnc" %in% algs) {

            print("Running D&C")         
            start <- proc.time()[1] 
            res_dnc <- r_dnc(ret)
            dnc_totaltime <- as.numeric(proc.time()[1] - start)
            print("Total time D&C")
            print(dnc_totaltime)

            print("Score from D&C")
            print(res_dnc$log_score)
            df_dnc <- data.frame(alg=c("dnc"), N = c(N), ub = c(ub), lb = c(lb), n = c(n), d = c(d), seed = c(i), 
                            totaltime = c(dnc_totaltime), 
                            max_particles = c(res_dnc$max_n_particles), tot_particles = c(res_dnc$tot_n_particles), 
                            scoretype=c(scoretype), am=c(as.numeric(am)), aw=c(format(aw)), chi=c(chi), 
                            edgepf=c(edgepf))
            df <- rbind(df, df_dnc)

          }
          if ("gobnilp" %in% algs) {
            # Compose the gibnilp.set file
            # Name the gobnilp.set file
            dir.create("gobnilp/") # weher the settings files are stored
            dir.create(paste0(argv$output_dir,"/gobnilp/"))
            gobnilp_conf_name <- paste0("gobnilp/", basename(name),  ".set")
            gobnilp_time_name <- paste0(argv$output_dir,"/gobnilp/", basename(name),  ".txt")
            setstr <- paste0('gobnilp/outputfile/scoreandtime = "', gobnilp_time_name, '"')

            write(setstr, gobnilp_conf_name)
            # Read the time  from the time file as a csv file
            
            # Using the singularity version of 
            output <- system(paste0("singularity exec gobnilp_4347c64.sif bash -c '/myappdir/gobnilp/bin/gobnilp ", "-g=", gobnilp_conf_name, " " , gobnilp_scores_filename, "'"), intern=TRUE)

            gobnilp_time <- read.csv(gobnilp_time_name, sep="\t", header = FALSE)

            df_gobnilp <- data.frame(alg=c("gobnilp"), N=c(N), ub = c(ub), lb = c(lb), n = c(n), d = c(d), seed = c(i), 
                            totaltime = c(gobnilp_time[["V2"]]), 
                            max_particles = c("NULL"), tot_particles = c("NULL"), 
                            scoretype=c(scoretype), am=c(as.numeric(am)), aw=c(format(aw)), chi=c(chi), 
                            edgepf=c(edgepf))
            
            df <- rbind(df, df_gobnilp)

        }            
          
          if ("orderpruner" %in% algs) {
            print("Running optimal order pruning")         
            start <- proc.time()[1] 
            res <- optimal_order(ret, c())
            op_totaltime <- as.numeric(proc.time()[1] - start)
            print("Total time orderpruner")
            print(op_totaltime)
            print(res)
            df_op <- data.frame(alg=c("orderpruner"),N = c(N), ub = c(ub), lb = c(lb), n = c(n), d = c(d), seed = c(i), 
                            totaltime = c(op_totaltime), 
                            max_particles = c(res$max_n_particles), tot_particles = c(res$tot_n_particles), 
                            scoretype=c(scoretype), am=c(as.numeric(am)), aw=c(format(aw)), chi=c(chi), 
                            edgepf=c(edgepf))
            df <- rbind(df_dnc, df_op)
            
          }          

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
