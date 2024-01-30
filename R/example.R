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
N = 100
p <- 32
dag <- randDAG(p, 3, method ="interEr", par1=4, par2=0.01, DAG = TRUE, weighted = FALSE, wFUN = list(runif, min=0.1, max=1))
#dag <- randDAG(p, 2, method ="er", par1=4, par2=0.01, DAG = TRUE, weighted = FALSE, wFUN = list(runif, min=0.1, max=1))
adjmat <- 1 * t(as(dag, "matrix") ) # transpose?

colnames(adjmat) <- paste0("X", 1:p)

G_true <- igraph::graph_from_adjacency_matrix(adjmat, mode="directed")

print("G_true:")
print(G_true)
print("components:")
print(igraph::components(G_true))

# Also ge the optimal estiamted DAG to check the componets

weight_mat <- adjmat
n_edges <- sum(adjmat)
weight_mat[which(weight_mat == 1)] <- wFUN(n_edges, lb = 0.25, ub = 1)
data <- data.frame(rmvDAG(weight_mat, N))

# filename <- "data/p20n300gaussdata.csv"
# # filename <- "data/asiadata.csv"
# # #filename <- "data/asiadata_double.csv"
# data <- read.csv(filename, check.names = FALSE)

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

# This is now working if aliases is empty
aliases <- lapply(cpp_friendly_scores$aliases, function(a) a + 1)

print("Creating diff matrices")
diff_matrices <- get_diff_matrices(cpp_friendly_scores$rowmaps, cpp_friendly_scores$scoretable, aliases, labels(data)[[2]])
H_min <- diff_matrices$H_min
H_max <- diff_matrices$H_max

#initial_suborder <- list(5,1,3)
#initial_suborder <- list()
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



#print("membership:")
#print(membership)


# Run order search on each component from H_min and see which components in Hmin
# That the parents belong to

# Fact: If a components exists in H_min then it is a component or a sub component in H_max.
optimal_components <- list()

# TODO: This should take the membership vector, not the graph
component_dependence <- function(membership, bidag_scores, cpp_friendly_scores) {
    # Vector of component numbers
    n_components <- max(membership)

    # Go through the components of H_min
    #for (component_id1 in seq(n_components)) {
    component_id1 <- 1

    component_score_sum <- 0    

    # Adjacentcy matrix for the component dependence
    comp_dep <- matrix(0, n_components, n_components)

    while (component_id1 <= n_components) {

        component1 <- seq(1, p)[membership == component_id1]
        initial_suborder <- seq(1, p)[-component1]

        # Run order search for the nodes in component_id1, with al the orher nodes as initial suborder, i.e
        # possible parents. Then check to which component the parents belong to.
        print("calculating optimal order for suborder:")
        print(component1)
        tmp <- optimal_order(cpp_friendly_scores, initial_suborder)

        print("optimal suborder and score:")
        print(tmp$suborder)
        # print suborder score
        print(tmp$suborder_cond_score)

        copmonent1_adjmat <- optimal_dag(bidag_scores, cpp_friendly_scores$space, tmp$order)

        G_opt <- igraph::graph_from_adjacency_matrix(copmonent1_adjmat, mode="directed")
        vertices1 <- igraph::V(G_opt)[membership == component_id1]$name
        
        # Now check to which component the parents belong to
        merged_components <- c()
        for (component_id2 in seq(n_components)) {    
            if (component_id1 == component_id2) next
            component2 <- seq(1, p)[membership == component_id2]
            # Get the edges between the two components
            vertices2 <- igraph::V(G_opt)[membership == component_id2]$name
            between_edges <- igraph::E(G_opt)[vertices1 %<-% vertices2]

            if (length(between_edges) > 0) {
                comp_dep[component_id2, component_id1] <- 1 # i.e. it is at least one edge/parent from comp 2 to 1
                # There are edges from component2 to component1
                print("components:")
                print(paste(component_id1,":"))
                print(component1)
                print(paste(component_id2, ":"))
                print(component2)
                
                print("The edges between the two components are:")
                print(between_edges)
            } 
        }
        component_id1 <- component_id1 + 1        
        
    }

    return(comp_dep)
    

}

get_cycles <- function(g) {
    g <- G_compdep
    Cycles = NULL
    for(v1 in igraph::V(g)) {
        for(v2 in igraph::neighbors(g, v1, mode="out")) {
            Cycles = c(Cycles, 
                lapply(igraph::all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p)))
        }
    }
    return(Cycles)
}


#print("Sum of the individual components scores:")
#print(component_score_sum)
membership <- igraph::components(G_H_min)$membership
#n_components <- igraph::components(G_H_min)$no

adjmat_compdep <- component_dependence(membership, bidag_scores, cpp_friendly_scores)

# merge components that are dependent on each other

print("Component dependence graph:")

G_compdep <- igraph::graph_from_adjacency_matrix(adjmat_compdep, mode="directed")
print(igraph::E(G_compdep))

Cycles <- get_cycles(G_compdep)
print("Cycles:")
print(Cycles)

# Make the cycle components bidirected. And take out those components.
n_components <- max(membership)
adjmat_compdep2 <- matrix(0, n_components, n_components)
print("G_compdep2:")
print(adjmat_compdep2)
for (cycle in Cycles) {
    print("cycle:")
    print(cycle)
    # Create edges in both directions for elements in the cycle
    for (i in seq(1, length(cycle))) {
        for (j in seq(1, length(cycle))) {
            if (i == j) next
            #print(paste(i, j))
            #print(paste(cycle[i], cycle[j]))
            adjmat_compdep2[cycle[i], cycle[j]] <- 1
            adjmat_compdep2[cycle[j], cycle[i]] <- 1
        }
    }
}
print("adjmat_compdep2:")
print(adjmat_compdep2)
# This is just to get the new merged components, where the cycles are merged.
G_compdep2 <- igraph::graph_from_adjacency_matrix(adjmat_compdep2, mode="undirected")
membership_comp <- igraph::components(G_compdep2)$membership
print("membership:")
print(membership_comp)
n_merged_components <- max(membership_comp)


## DAG for the merged components.
merged_components_adjmat <- matrix(0, n_merged_components, n_merged_components)

# Now fill in with the non cycle edges in adjmat_compdep2
# Go through edges in the G_compdep.

for (i in seq(1, n_components)) {
    for (j in seq(1, n_components)) {
        #print(paste(i, j))
        # check if the edge i->j or i<-j exists in G_compdep
        if ((adjmat_compdep[i, j] == 0) && (adjmat_compdep[j,i]==0)) next
        if (i == j) next
        merged_i <- membership_comp[i]
        merged_j <- membership_comp[j]    
        #print("merged components:")
        #print(paste(merged_i, merged_j))
        if (merged_i == merged_j) next # same merged component, ie no edge

        # Here we know that i and j are not merged
        if ((adjmat_compdep[i, j] == 1) && (adjmat_compdep[j,i]==0)) {
            merged_components_adjmat[merged_i, merged_j] <- 1
        }

        if((adjmat_compdep[i, j] == 0) && (adjmat_compdep[j,i]==1)) {
            merged_components_adjmat[merged_j, merged_i] <- 1
        }
        # if not both i and j are in the same component, 
        # check the new components of i and j and add the 
    }
}

print("merged_components_adjmat:")
print(merged_components_adjmat)
merged_components_graph <- igraph::graph_from_adjacency_matrix(merged_components_adjmat, mode="directed")

sorted_compnents <- igraph::topo_sort(merged_components_graph, mode= "in")
print("sorted_compnents:")
print(sorted_compnents)

# TODO: Still need the directed part for the order!
# check for the edges that are not part of a cycle and add these too

# Now, go trough the nodes in the original graph, and check their components
# and cchek in the merged components which components they belong to.

ultimate_membership <- rep(0, p)
for (i in seq(1, p)) {    
    #print(paste("node:", i))    
    original_component <- membership[i]
    merged_component <- membership_comp[original_component]
    ultimate_membership[i] <- merged_component
}
print("ultimate_membership:")
print(ultimate_membership)


# Remove the remaining directed edges

# Get the components of the new graph


## The current components

## The new components
# Create these by merging the current components guided by the


# print(igraph::as_ids(igraph::E(G_opt)[1:6 %<-% 7:16]))

# print(igraph::E(G_opt))

## For the components of H_max run the following:

## For the components of H_min run the following:
## 1. Get the components of H_min
## 1.1. For each component, run the order search and see if the nodes in the optimal
## order has parents outside the component. If so, add them to the component. 
## 1.2. Rerun the order search and see if the nodes in the optimal.