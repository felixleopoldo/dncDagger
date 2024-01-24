source("R/scoring.R")
source("R/opruner.r") 

#filename <- "data/p20n300gaussdata.csv"
#filename <- "data/asiadata.csv"
filename <- "data/asiadata_double.csv"
data <- read.csv(filename, check.names = FALSE)
p <- ncol(data)
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


aliases <- lapply(cpp_friendly_scores$aliases, function(a) a + 1)
diff_matrices <- get_diff_matrices(cpp_friendly_scores$rowmaps, cpp_friendly_scores$scoretable, aliases, labels(data)[[2]])
H_min <- diff_matrices$H_min
H_max <- diff_matrices$H_max

initial_suborder <- list(5,1,3)

print("initial suborder:")
print(initial_suborder)

opr <- optimal_order(cpp_friendly_scores, initial_suborder)
adjmat <- optimal_dag(bidag_scores, cpp_friendly_scores$space, opr$order)
colnames(adjmat) <- colnames(data)
print(opr)

H_min_adj <- (H_min> 0)*1
H_max_adj <- (H_max> 0)*1

G_opt <- igraph::graph_from_adjacency_matrix(adjmat, mode="directed")
G_H_min <- igraph::graph_from_adjacency_matrix(H_min_adj, mode="directed")
G_H_max <- igraph::graph_from_adjacency_matrix(H_max_adj, mode="directed")

print("H_min:")
print(H_min_adj)

print("H_max:")
print(H_max_adj)

# Compontents of H_min. Possible component in G.
print("Components of H_min:")
print(igraph::components(G_H_min)$membership)

# Compontents of H_max. Ome componetn of this may contain several componets of G.
print("Components of H_max:")
print(igraph::components(G_H_max)$membership)

# Should be able to run order opt on suborder.

# Vector of component numbers
membership <- igraph::components(G_H_min)$membership
n_components <- igraph::components(G_H_min)$no

# Run order search on each component from H_min and see which components in Hmin
# That the parents belong to

for (component_id1 in seq(n_components)) {
    component1 <- seq(1, p)[membership == component_id1]
    vertices1 <- igraph::V(G_opt)[membership == component_id1]$name

    initial_suborder <- seq(1, p)[-component1]-1

    print("initial suborder:")
    print(initial_suborder)
    # Run order search for the nodes in component_id1
    tmp <- optimal_order(cpp_friendly_scores, initial_suborder)
    copmonent1_adjmat <- optimal_dag(bidag_scores, cpp_friendly_scores$space, tmp$order)

    # Now check to which component the parents belong to
    for (component_id2 in seq(n_components)) {    
        if (component_id1 == component_id2) next

        # Run order search for the nodes in component_id1 having
        # parents in component_id2 as possible parents.
        # This should only run once actually.

        vertices2 <- igraph::V(G_opt)[membership == component_id2]$name
        between_edges <- igraph::E(G_opt)[vertices1 %<-% vertices2]

        if (length(between_edges) > 0) {
            print("edges between groups")
            print(between_edges)
            print(igraph::as_ids(between_edges))
            print("merging groups")
            #membership[membership == component_id2] <- component_id1
        }
    }
}

print("edges between groups")
print(igraph::E(G_opt)[1:6 %->% 7:16])
print(igraph::E(G_opt)[1:6 %<-% 7:16])

print(igraph::as_ids(igraph::E(G_opt)[1:6 %<-% 7:16]))
print(igraph::as_ids(igraph::E(G_opt)[1:6 %--% 7:16]))
print(igraph::as_ids(igraph::E(G_opt)[1:6 %->% 7:16]))

print(igraph::E(G_opt))

## For the components of H_max run the following:

## For the components of H_min run the following:
## 1. Get the components of H_min
## 1.1. For each component, run the order search and see if the nodes in the optimal
## order has parents outside the component. If so, add them to the component. 
## 1.2. Rerun the order search and see if the nodes in the optimal.