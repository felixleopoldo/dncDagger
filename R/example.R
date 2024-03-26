source("R/scoring.R")
source("R/opruner.r") 

set.seed(2)
# Generate data
N = 50
ndim <- 8
#dag <- randDAG(ndim, 3, method ="interEr", par1=5, par2=0.01, DAG = TRUE, weighted = FALSE, wFUN = list(runif, min=0.1, max=1))
dag <- randDAG(ndim, 3, method ="interEr", par1=2, par2=0.01, DAG = TRUE, weighted = FALSE, wFUN = list(runif, min=0.1, max=1))
#dag <- randDAG(ndim, 2, method ="er", par1=4, par2=0.01, DAG = TRUE, weighted = FALSE, wFUN = list(runif, min=0.1, max=1))
adjmat <- 1 * t(as(dag, "matrix") ) # transpose?

colnames(adjmat) <- paste0("X", 1:ndim)

G_true <- igraph::graph_from_adjacency_matrix(adjmat, mode="directed")

print("G_true:")
print(G_true)
print("components:")
print(igraph::components(G_true))

# Also ge the optimal estiamted DAG to check the components.
weight_mat <- adjmat
n_edges <- sum(adjmat)
weight_mat[which(weight_mat == 1)] <- wFUN(n_edges, lb = 0.25, ub = 1)
dataset <- data.frame(rmvDAG(weight_mat, N))
write.csv(dataset, "data/testing.csv", row.names=FALSE)
#filename <- "data/p20n300gaussdata.csv"
#filename <- "data/asiadata.csv"
#filename <- "data/asiadata_double.csv"
#filename <- "data/n=20d=0_seed=1_lb=0.25_ub=1_N=300.csv"
filename <- "data/testing.csv"
data <- read.csv(filename, check.names = FALSE)
ndim <- ncol(dataset)

set.seed(1)
print("getting cpp friendly bidag scores")
# cpp_friendly_scores <- get_plus1_score_essentials_for_cpp(bidag_scores, plus1it=2, iterations=NULL) # from iterativeMCMC

cpp_friendly_scores <<- get_scores(filename, scoretype="bge", 
                                  bgepar=list(am=0.1, aw=NULL), 
                                  bdepar=list(chi=0.5, edgepf=2), # one of these should be ignored
                                  plus1it=2)

# print("running optimal order pruning")
# initial_suborder <- c()
# start <- proc.time()[1]
# opr <- optimal_order(cpp_friendly_scores, initial_suborder)
# print("optimal order:")
# print(opr)
# adjmat <- optimal_dag(cpp_friendly_scores$bidag_scores, cpp_friendly_scores$space, opr$order)
# G_opt <- igraph::graph_from_adjacency_matrix(adjmat, mode="directed")
# totaltime <- as.numeric(proc.time()[1] - start)
# print("opruner: Total time")
# print(totaltime)

# print("G_opt:")
# print(G_opt)
# # colnames(adjmat) <- colnames(data)
# print(opr)
Rprof()
start <- proc.time()[1]  
print("running dnc")
bidag_scores <<- cpp_friendly_scores$bidag_scores
res <- dnc2()#cpp_friendly_scores, cpp_friendly_scores$bidag_scores)
print("dnc: Total time")
totaltime <- as.numeric(proc.time()[1] - start)
print(totaltime)

print("Time spent on finding the optimal DAG")
print(res$tot_order_to_dag_time)

print("Total time after subtraction")
print(totaltime - res$tot_order_to_dag_time)

print("Tot get order time")
print(res$tot_order_time)

dag <- igraph::graph_from_adjacency_matrix(res$adjmat, mode="directed")
png("dnc.png")
plot(dag)
dev.off() 