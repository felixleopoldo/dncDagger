source("R/scoring.R")
source("R/opruner.r") 

set.seed(2)
# Generate data
N = 100
ndim <- 32
dag <- randDAG(ndim, 3, method ="interEr", par1=8, par2=0.01, DAG = TRUE, weighted = FALSE, wFUN = list(runif, min=0.1, max=1))
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
data <- data.frame(rmvDAG(weight_mat, N))

#filename <- "data/p20n300gaussdata.csv"
#filename <- "data/asiadata.csv"
#filename <- "data/asiadata_double.csv"
#data <- read.csv(filename, check.names = FALSE)
ndim <- ncol(data)

scoretype <- "bge"
bgepar <- list(am=1, aw=NULL)

if (scoretype =="bge") {
    bidag_scores <- BiDAG::scoreparameters(scoretype = scoretype, data, bgepar = bgepar)
} else if (scoretype == "bde") {
    bidag_scores <- BiDAG::scoreparameters(scoretype = scoretype, data[-1, ], bdepar = bdepar)
}

set.seed(1)
print("getting cpp friendly bidag scores")
cpp_friendly_scores <- get_plus1_score_essentials_for_cpp(bidag_scores, plus1it=2, iterations=NULL) # from iterativeMCMC

# initial_suborder <- c()
# opr <- optimal_order(cpp_friendly_scores, initial_suborder)
# print("optimal order:")
# print(opr)
# adjmat <- optimal_dag(bidag_scores, cpp_friendly_scores$space, opr$order)
# G_opt <- igraph::graph_from_adjacency_matrix(adjmat, mode="directed")
# print("G_opt:")
# print(G_opt)
# # colnames(adjmat) <- colnames(data)
# print(opr)

res <- dnc(cpp_friendly_scores, bidag_scores)

print(res)
dag <- igraph::graph_from_adjacency_matrix(res$adjmat, mode="directed")
png("dnc.png")
plot(dag)
dev.off() 