source("R/scoring.R")
source("R/opruner.r") 

wFUN <- function(m, lb, ub) {
  # This function gives edges weights between the bounds
  # with both positive and negative signs
  runif(m, lb, ub) * sample(c(-1, 1), m, replace = TRUE)
}

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



set.seed(1)
# Generate data
N = 1000
p <- 40
#dag <- randDAG(p, 4, method ="interEr", par1=4, par2=0.01, DAG = TRUE, weighted = FALSE, wFUN = list(runif, min=0.1, max=1))
dag <- randDAG(p, 2, method ="er", par1=4, par2=0.01, DAG = TRUE, weighted = FALSE, wFUN = list(runif, min=0.1, max=1))
adjmat <- 1 * t(as(dag, "matrix") ) # transpose?

print("adjmat:")
print(adjmat)

G_true <- igraph::graph_from_adjacency_matrix(adjmat, mode="directed")

print("G_true:")
print(G_true)
print("components:")
print(igraph::components(G_true))

weight_mat <- adjmat
n_edges <- sum(adjmat)
weight_mat[which(weight_mat == 1)] <- wFUN(n_edges, lb = 0.25, ub = 1)
data <- data.frame(rmvDAG(weight_mat, N))

# filename <- "data/p20n300gaussdata.csv"
# # filename <- "data/asiadata.csv"
# # #filename <- "data/asiadata_double.csv"
# data <- read.csv(filename, check.names = FALSE)


p <- ncol(data)

print("p:") 
print(p)

#print("data:")
#print(head(data.frame(data)))
#print(data)


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

# This is now working if aliases is empty
aliases <- lapply(cpp_friendly_scores$aliases, function(a) a + 1)

print("Creating diff matrices")
diff_matrices <- get_diff_matrices(cpp_friendly_scores$rowmaps, cpp_friendly_scores$scoretable, aliases, labels(data)[[2]])
H_min <- diff_matrices$H_min
H_max <- diff_matrices$H_max

#initial_suborder <- list(5,1,3)
initial_suborder <- list()

# print("initial suborder:")
# print(initial_suborder)

# opr <- optimal_order(cpp_friendly_scores, initial_suborder)
# ## adjmat <- optimal_dag(bidag_scores, cpp_friendly_scores$space, opr$order)
# # colnames(adjmat) <- colnames(data)
# print(opr)

H_min_adj <- (H_min> 0)*1
H_max_adj <- (H_max> 0)*1

#G_opt <- igraph::graph_from_adjacency_matrix(adjmat, mode="directed")
G_H_min <- igraph::graph_from_adjacency_matrix(H_min_adj, mode="directed")
G_H_max <- igraph::graph_from_adjacency_matrix(H_max_adj, mode="directed")

# print("H_min:")
# print(H_min_adj)

# print("H_max:")
# print(H_max_adj)

# Compontents of H_min. Possible component in G.
print("Components of H_min:")
print(igraph::components(G_H_min))

# Compontents of H_max. Ome componetn of this may contain several componets of G.
print("Components of H_max:")
print(igraph::components(G_H_max))

# Should be able to run order opt on suborder.

# Vector of component numbers
membership <- igraph::components(G_H_min)$membership
n_components <- igraph::components(G_H_min)$no

#print("membership:")
#print(membership)


# Run order search on each component from H_min and see which components in Hmin
# That the parents belong to

# Fact: If a components exists in H_min then it is a component or a sub component in H_max.
optimal_components <- list()


# Go through the components of H_min
#for (component_id1 in seq(n_components)) {
component_id1 <- 1
while (component_id1 <= n_components) {
    
    component1 <- seq(1, p)[membership == component_id1]
    initial_suborder <- seq(1, p)[-component1]

    # Run order search for the nodes in component_id1, with al the orher nodes as initial suborder, i.e
    # possible parents. Then check to which component the parents belong to.
    tmp <- optimal_order(cpp_friendly_scores, initial_suborder)

    print("optimal suborder:")
    print(tmp$suborder)

    copmonent1_adjmat <- optimal_dag(bidag_scores, cpp_friendly_scores$space, tmp$order)
    G_opt <- igraph::graph_from_adjacency_matrix(copmonent1_adjmat, mode="directed")
    vertices1 <- igraph::V(G_opt)[membership == component_id1]$name
    
    # Now check to which component the parents belong to
    for (component_id2 in seq(n_components)) {    
        if (component_id1 == component_id2) next

        # Get the edges between the two components
        vertices2 <- igraph::V(G_opt)[membership == component_id2]$name
        between_edges <- igraph::E(G_opt)[vertices1 %<-% vertices2]

        if (length(between_edges) > 0) {
            print("edges between groups")
            print(igraph::as_ids(between_edges))
            print("merging groups")
            #n_components <- n_components - 1
            #membership[membership == component_id2] <- component_id1
        } else {
            print("no edges between these groups")
        }
    }
    component_id1 <- component_id1 + 1

}


# print(igraph::as_ids(igraph::E(G_opt)[1:6 %<-% 7:16]))

# print(igraph::E(G_opt))

## For the components of H_max run the following:

## For the components of H_min run the following:
## 1. Get the components of H_min
## 1.1. For each component, run the order search and see if the nodes in the optimal
## order has parents outside the component. If so, add them to the component. 
## 1.2. Rerun the order search and see if the nodes in the optimal.