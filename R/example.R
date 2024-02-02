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

# TODO: This should take the membership vector, not the graph
component_dependence <- function(membership, bidag_scores, cpp_friendly_scores, input_opt_suborders=NULL) {
    # Vector of component numbers
    n_components <- max(membership)
    #print(paste("n_components:", n_components))
    #print(input_opt_suborders)
    # Go through the components of H_min
    #for (component_id1 in seq(n_components)) {
    component_id1 <- 1

    component_score_sum <- 0    

    # Adjacentcy matrix for the component dependence
    comp_dep <- matrix(0, n_components, n_components)
    # list of order for each component
    opt_sub_orders <- list()
    opt_sub_orders[[1]] <- 1234
    while (component_id1 <= n_components) {

        component1 <- seq(1, p)[membership == component_id1]
        initial_suborder <- seq(1, p)[-component1]

        # Run order search for the nodes in component_id1, with al the orher nodes as initial suborder, i.e
        # possible parents. Then check to which component the parents belong to.
        print(paste("******* calculating optimal order for suborder:", component_id1, "********"))
        print(component1)

        # Check if we already have calculations for this components
        component1_exists <- FALSE
        if (!is.null(input_opt_suborders)) {
            # Go through all the components in input_opt_suborders and check if they are component1
            #print("--------- Checking if component exists in input_opt_suborders")            
            for (comp_id in seq(1, length(input_opt_suborders))) {
               #print(paste("input_opt_suborders[[comp_id]]$suborder:"))
               # print(sort(input_opt_suborders[[comp_id]]$suborder))                
                
                if (setequal(input_opt_suborders[[comp_id]]$suborder, component1)) {
                    #print("Found component in input_opt_suborders")
                    # print(input_opt_suborders[[comp_id]]$suborder)
                    # print(input_opt_suborders[[comp_id]]$score)
                    # print(opt_sub_orders)
                    # print(paste("component_id1:", component_id1))
                    

                    # print(opt_sub_orders)
                    # print("trying to assign")
                    opt_sub_orders[[component_id1]] <- input_opt_suborders[[comp_id]]

                    component1_exists <- TRUE
                    break
                }
            }
        }
        # If we have a new component, calculate the optimal order for it
        if (component1_exists==FALSE) {
            tmp <- optimal_order(cpp_friendly_scores, initial_suborder)        
            
            copmonent1_adjmat <- optimal_dag(bidag_scores, cpp_friendly_scores$space, tmp$order) #suborder?
            G_opt <- igraph::graph_from_adjacency_matrix(copmonent1_adjmat, mode="directed")            
            #print("component_id1:")
            #print(component_id1)
            opt_sub_orders[[component_id1]] <- list("suborder"=tmp$suborder, 
                                                    "score"=tmp$suborder_cond_score,
                                                    "dag"=G_opt)
        
            #print("optimal suborder and score:")
            #print(tmp$suborder)
            # print suborder score
            #print(tmp$suborder_cond_score)
        }

        # In any case we have to evaluate where the potential parents come from, 
        # as the component ids might have been renamed.
        vertices1 <- igraph::V(opt_sub_orders[[component_id1]]$dag)[membership == component_id1]$name
        
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
                print(paste("Edges from component:", component_id2))
                print(component2)
                print(between_edges)
            } 
        }
        component_id1 <- component_id1 + 1        
        
    }

    # note that opt_sub_orders will only be valied after the las round when all the merging is done.
    return(list("comp_dep"=comp_dep, "opt_sub_orders"=opt_sub_orders)) 

}

get_cycles <- function(g) {
    #g <- G_compdep
    Cycles = NULL
    for(v1 in igraph::V(g)) {
        for(v2 in igraph::neighbors(g, v1, mode="out")) {
            Cycles = c(Cycles, 
                lapply(igraph::all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p)))
        }
    }
    return(Cycles)
}

merged_neig_cycles <- function(adjmat_compdep2){
    G_compdep <- igraph::graph_from_adjacency_matrix(adjmat_compdep2, mode="directed")
    Cycles <- get_cycles(G_compdep)
    #print("Cycles:")
    #print(Cycles)


    # Make the cycle components bidirected. And take out those components.
    n_components <- ncol(adjmat_compdep2)
    #n_components <- max(membership) change this to the above.
    adjmat_compdep2 <- matrix(0, n_components, n_components)

    for (cycle in Cycles) {
        #print("cycle:")
        #print(cycle)
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
    #print("adjmat_compdep2:")
    #print(adjmat_compdep2)
    # This is just to get the new merged components, where the cycles are merged.
    G_compdep2 <- igraph::graph_from_adjacency_matrix(adjmat_compdep2, mode="undirected")
    membership_comp <- igraph::components(G_compdep2)$membership

    return(membership_comp)
}

merged_component_dependencies <- function(adjmat_compdep, membership_comp){
    n_merged_components <- max(membership_comp)
    n_components <- length(membership_comp)

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

    return(merged_components_graph)
}

merged_components_membership <- function(membership, membership_comp){
    # Translates the membership vector of the original graph to the new merged components
    print("Determin new node memberships in the merged components:")
    ultimate_membership <- rep(0, p)
    for (i in seq(1, p)) {    
        #print(paste("node:", i))
        original_component <- membership[i]
        #print(paste("original_component:", original_component))
        merged_component <- membership_comp[original_component]
        #print(paste("merged_component:", merged_component))
        ultimate_membership[i] <- merged_component
        #print(paste("ultimate_membership:", ultimate_membership[i]))
    }
    return(ultimate_membership)
}


set.seed(1)
# Generate data
N = 100
p <- 16
dag <- randDAG(p, 2, method ="interEr", par1=4, par2=0.01, DAG = TRUE, weighted = FALSE, wFUN = list(runif, min=0.1, max=1))
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
# filename <- "data/asiadata.csv"
#filename <- "data/asiadata_double.csv"
#data <- read.csv(filename, check.names = FALSE)

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

# initial_suborder <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,26,27,28,29,30,31) +1 
#initial_suborder <- c()
#opr <- optimal_order(cpp_friendly_scores, initial_suborder)
#print("optimal order:")
#print(opr)


# ## adjmat <- optimal_dag(bidag_scores, cpp_friendly_scores$space, opr$order)
#G_opt <- igraph::graph_from_adjacency_matrix(adjmat, mode="directed")
# # colnames(adjmat) <- colnames(data)
# print(opr)

H_min_adj <- (H_min> 0)*1
H_max_adj <- (H_max> 0)*1
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


# Run order search on each component from H_min and see which components in Hmin
# That the parents belong to

#print("Sum of the individual components scores:")
#print(component_score_sum)

membership <- igraph::components(G_H_min)$membership

div_n_conquer <- function(membership, bidag_scores, cpp_friendly_scores, input_opt_suborders=NULL) {
    
    # Get the component dependence graph    
    ret <- component_dependence(membership, bidag_scores, cpp_friendly_scores, input_opt_suborders=input_opt_suborders)
    #print(ret)
    adjmat_compdep <- ret$comp_dep
    G_compdep <- igraph::graph_from_adjacency_matrix(adjmat_compdep, mode="directed")
    print("Component dependence graph before possible merging:")
    print(igraph::E(G_compdep))

    # Merge components that are dependent on each other in a cycle and neigboring cycles.
    # This is for the graph of comopnents. So its the new components membership after merging the cycles.
    print("Component dependence graph:")
    membership_comp <- merged_neig_cycles(adjmat_compdep)
    print("membership after merging:")
    print(membership_comp)
    
    # Get dependency graph for the merged components
    merged_components_graph <- merged_component_dependencies(adjmat_compdep, membership_comp)

    # Now, go trough the nodes in the original graph, and check their components
    # and check in the merged components which components they belong to.
    ultimate_membership <- merged_components_membership(membership, membership_comp)
    print("ultimate_membership:")
    print(ultimate_membership)
    
    # Could be done aboe but like this we shold be able to verify that the merging is correct.
    # Here we can chag if there where no merged components, i.e. no cycles.
    # By checking if the number of components is the same as the number of merged components.
    if (max(membership_comp) == length(membership_comp)) {
        print("No cycles, returning the original order")
        return(list("membership"=membership, 
                    "merged_components_graph"= G_compdep, 
                    "no_cycles"=TRUE,
                    "component_order"=ret,
                    "opt_sub_orders"=ret$opt_sub_orders))
    }
    
    return(list("membership"=ultimate_membership,
                "merged_components_graph"= merged_components_graph,
                "no_cycles"=FALSE,
                "component_order"=NULL,
                "opt_sub_orders"=ret$opt_sub_orders))
}

concat_suborders <- function(ultimate_membership, merged_components_graph, component_order){
    ## Paste the pieaces of suborders acccording to the topological order of merged components.    
    sorted_components <- igraph::topo_sort(merged_components_graph, mode= "in")
    print("sorted_components:")
    print(sorted_components)

    full_order <- c()
    total_score <- 0
    for (comp_id in sorted_components) {
        suborder <- component_order$opt_sub_orders[[comp_id]]$suborder
        print(paste("suborder:", comp_id))
        print(suborder)
        print("score:")
        print(component_order$opt_sub_orders[[comp_id]]$score)
        full_order <- c(full_order, suborder) #maybe the other way around
        total_score <- total_score + component_order$opt_sub_orders[[comp_id]]$score
    }
    return(list("order"=full_order, "score"=total_score))
}


round <- 1
opt_sub_orders <- NULL
while (TRUE) {   
    print(paste("############# Round ",round," #################"))
    ret <- div_n_conquer(membership, bidag_scores, cpp_friendly_scores, opt_sub_orders)
    ultimate_membership <- ret$membership
    merged_components_graph <- ret$merged_components_graph
    if (ret$no_cycles) { # no cycles to merge.
        res <- concat_suborders(ultimate_membership, merged_components_graph, ret$component_order)
        print("optorder:")
        print(res$order)
        print(res$score)
        break
    }
    opt_sub_orders <- ret$opt_sub_orders
    membership <- ultimate_membership # Update the membership to the new merged components
    round <- round + 1
}
