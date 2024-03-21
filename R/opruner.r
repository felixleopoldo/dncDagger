library("Rcpp")
# library("Jmisc")
# library("pcalg")
# library("ggplot2")
# library("testit")
library(argparser)
library("BiDAG")
Sys.setenv("PKG_CXXFLAGS" = "-Wall -pipe -Wno-unused -pedantic -Wall -L /usr/lib/R/lib -l R -L /usr/local/lib/R/site-library/RInside/lib/ -l RInside -Wl,-rpath,/usr/local/lib/R/site-library/RInside/lib  -I /usr/local/lib/R/site-library/RInside/include/ -I /usr/local/lib/R/site-library/Rcpp/include/ -I /usr/share/R/include/ -std=c++17 -O3")
sourceCpp("include/opruner_right.cpp",  verbose = TRUE)


optimal_order <- function(cpp_friendly_scores, initial_sub_order){

    # Check if the initial suborder is empty
    print("Running optimal order")
    ret <- r_opruner_right(cpp_friendly_scores, initial_sub_order)
    ret$suborder_cond_score <- 0
    ret$order <- ret$order + 1 # since c++ enumerate nodes from 0
    p <- length(ret$node_scores)
    suborder <- ret$order[1:(p-length(initial_sub_order))] # The suborder is the firs part of the order
    # print(paste("optimal_order: suborder:"))
    # print(paste(suborder))
    # print(paste("optimal_order: initial_sub_order:"))
    # print(paste(initial_sub_order))
    # print("optimal order")
    # print(ret$order)
    # #initial_sub_order

    for (node in suborder) {
        ret$suborder_cond_score <- ret$suborder_cond_score + ret$node_scores[node]
    }
    ret$suborder <- suborder
    return(ret)
}

optimal_dag <- function(bidag_scores, space, order) {
  omcmc = BiDAG::orderMCMC(bidag_scores, plus1=TRUE, MAP=TRUE, iterations=1,
                           stepsave=1, startorder=order, startspace=space)

  return(omcmc$DAG)
}

structure_components <- function(isolated_comp_membership, separable_comp_membership){
    isocomps <- list()
    membership <- separable_comp_membership
    for (i in seq(1, max(isolated_comp_membership))) {
        nodes <- which(isolated_comp_membership == i)
        isocomps[[i]] <- list()
        isocomps[[i]]$nodes <- nodes
        isocomps[[i]]$subcompids <- c()
        isocomps[[i]]$subcomps <- list()
        isocomps[[i]]$connected <- FALSE # This means that the components Could actually be connected
    }

    # Find the sub components in H_min
    #print("Finding super components for:")
    for (j in seq(1, max(membership))) {

        subnodes <- which(membership == j)

        if(length(unique(isolated_comp_membership[subnodes])) == 1) {
            isocompid <- isolated_comp_membership[subnodes][1]
            isocomps[[isocompid]]$subcompids <- c(j, isocomps[[isocompid]]$subcompids)

            k <- length(isocomps[[isocompid]]$subcomps)

            isocomps[[isocompid]]$subcomps[[k+1]] <- list("nodes" = subnodes,
                                                          "score" = NULL,
                                                          "opt_adjmat" = NULL,
                                                          "subadjmat" = NULL,
                                                          "suborder" = NULL)
            isocomps[[isocompid]]$max_n_particles <- 0
            isocomps[[isocompid]]$tot_n_particles <- 0
            isocomps[[isocompid]]$tot_order_to_dag_time <- 0
            isocomps[[isocompid]]$tot_order_time <- 0

            # Check if it is the same, ie the component is isolated and connected.
            if(length(subnodes) == length(isocomps[[isocompid]]$nodes)){
                isocomps[[isocompid]]$connected <- TRUE
            }
        }
    }
    

    return(isocomps)
}

dnc <- function(cpp_friendly_scores, bidag_scores) {
    #bidag_scores <- cpp_friendly_scores$bidag_scores
    # This is now working if aliases is empty
    aliases <- lapply(cpp_friendly_scores$aliases, function(a) a + 1)
    var_labels <- cpp_friendly_scores$labels
    print("Creating diff matrices")
    diff_matrices <- get_diff_matrices(cpp_friendly_scores$rowmaps, cpp_friendly_scores$scoretable, aliases, var_labels)
    H_min <- diff_matrices$H_min
    H_max <- diff_matrices$H_max
    p <- nrow(H_min)
    H_min_adj <- (H_min> 0)*1
    H_max_adj <- (H_max> 0)*1
    G_H_min <- igraph::graph_from_adjacency_matrix(H_min_adj, mode="directed")
    G_H_max <- igraph::graph_from_adjacency_matrix(H_max_adj, mode="directed")
    # print("H_min:")
    # print(H_min_adj)
    # print("H_max:")
    # print(H_max_adj)

    # Components of H_min. Possible component in G.
    print("Components of H_min:")
    print(igraph::components(G_H_min))

    # Components of H_max.
    # One componeny of this may contain several componets of G.
    # We call them isolated components.
    print("Components of H_max:")
    print(igraph::components(G_H_max))

    # Run order search on each component from H_min and see which components in Hmin
    # That the parents belong to

    membership <- igraph::components(G_H_min)$membership
    isolated_comp_membership <- igraph::components(G_H_max)$membership

    round <- 1
    opt_sub_orders <- NULL
    max_n_particles <- 0
    tot_n_particles <- 0
    tot_order_to_dag_time <- 0

    # maybe this loop or the div_n_conquer should be for each isolated component instead?
    while (TRUE) {
        print(paste("############# Round ",round," #################"))


        ret <- restructure_components(membership, isolated_comp_membership, bidag_scores, cpp_friendly_scores, opt_sub_orders)
        ultimate_membership <- ret$membership
        merged_components_graph <- ret$merged_components_graph
        max_n_particles <- max(max_n_particles, ret$max_n_particles)
        tot_n_particles <- tot_n_particles + ret$tot_n_particles
        #print(ret$tot_order_to_dag_time)
        tot_order_to_dag_time <- tot_order_to_dag_time + ret$tot_order_to_dag_time

        if (ret$no_cycles) { # no cycles to merge.
            # Concatenate subboorder according to the DAG of component dependencies
            res <- concat_suborders(ultimate_membership, merged_components_graph, ret$component_order)
            # join the sub matrices
            full_adjmat <- matrix(0, p, p)
            for (suborder in ret$opt_sub_orders) {
                full_adjmat <- full_adjmat + suborder$subadjmat
            }
            #fulldag <- igraph::graph_from_adjacency_matrix(full_adjmat, mode="directed")

            res$adjmat <- full_adjmat
            #res$dag <- fulldag
            res$max_n_particles <- max_n_particles
            res$tot_n_particles <- tot_n_particles
            res$tot_order_to_dag_time <- tot_order_to_dag_time
            return(res)
        }

        opt_sub_orders <- ret$opt_sub_orders
        membership <- ultimate_membership # Update the membership to the new merged components
        round <- round + 1
    }
}

dnc2 <- function(cpp_friendly_scores, bidag_scores) {
    #bidag_scores <- cpp_friendly_scores$bidag_scores
    # This is now working if aliases is empty
    aliases <- lapply(cpp_friendly_scores$aliases, function(a) a + 1)
    var_labels <- cpp_friendly_scores$labels
    print("Creating diff matrices")
    diff_matrices <- get_diff_matrices(cpp_friendly_scores$rowmaps, cpp_friendly_scores$scoretable, aliases, var_labels)
    H_min <- diff_matrices$H_min
    H_max <- diff_matrices$H_max
    p <- nrow(H_min)
    H_min_adj <- (H_min> 0)*1
    H_max_adj <- (H_max> 0)*1
    G_H_min <- igraph::graph_from_adjacency_matrix(H_min_adj, mode="directed")
    G_H_max <- igraph::graph_from_adjacency_matrix(H_max_adj, mode="directed")
    # print("H_min:")
    # print(H_min_adj)
    # print("H_max:")
    # print(H_max_adj)

    # Components of H_min. Possible component in G.
    print("Components of H_min:")
    print(igraph::components(G_H_min))

    # Components of H_max.
    # One componeny of this may contain several componets of G.
    # We call them isolated components.
    print("Components of H_max:")
    print(igraph::components(G_H_max))

    # Run order search on each component from H_min and see which components in Hmin
    # That the parents belong to

    membership <- igraph::components(G_H_min)$membership
    isolated_comp_membership <- igraph::components(G_H_max)$membership

    round <- 1
    opt_sub_orders <- NULL
    max_n_particles <- 0
    tot_n_particles <- 0
    tot_order_to_dag_time <- 0

    # Order the components in H_min after the isolated components (induced by H_max).
    start <- proc.time()[1]
    isocomps <- structure_components(isolated_comp_membership, membership)
    totaltime <- as.numeric(proc.time()[1] - start)
    print("Total time for finding super components:")
    print(totaltime)

    # We initialize the process by computing the component dependencies once for each
    # isolated component.
    start <- proc.time()[1]
    #todagtime <- 0
    #tot_order_time <- 0
    for (i in seq(1, max(isolated_comp_membership))) {
        print(paste("Component:", i))
        isocomps[[i]] <- component_dependence2(isocomps[[i]], bidag_scores, cpp_friendly_scores)
        # Update isocomps costs(time and particle counts)   
        isocomps$max_n_particles <- max(isocomps$max_n_particles, isocomps[[i]]$max_n_particles)
        isocomps$tot_n_particles <- isocomps$tot_n_particles + isocomps[[i]]$tot_n_particles
        isocomps$tot_order_to_dag_time <- isocomps$tot_order_to_dag_time + isocomps[[i]]$tot_order_to_dag_time

        #todagtime <- todagtime + isocomps[[i]]$tot_order_to_dag_time
        #tot_order_time <- tot_order_time + isocomps[[i]]$tot_order_time
    }
    totaltime <- as.numeric(proc.time()[1] - start)

    #print(isocomps)
    # print("Total time for component dependence:")
    # print(totaltime)
    # print("Total time for order to dag:")
    # print(todagtime)
    # print("subtracted time")
    # print(totaltime - todagtime)
    # print("tot get order time")
    # print(tot_order_time)
    # print("returned")
    # print(xx)
    # as isloated components doesnt have dependencies outside the component, they should ideally be treated
    # without the need for the other variables.

    # The divide and conquer procedure should continue until all components in
    # isocomp are connected. (If we decide to split up if they are not connected.)

    # maybe this loop or the div_n_conquer should be for each isolated component instead?
    print("Starting divide and conquer")
    print("number of isocomps")
    print(max(isolated_comp_membership))
    for (i in seq(1, max(isolated_comp_membership))) {
        
        print(paste("***********treating isolated component", i))
        #print(isocomps[[i]])
        # print nodes in the subcomponents 
        for (sc in isocomps[[i]]$subcomps) print(sc$nodes)

        if (isocomps[[i]]$connected) {
            print("Its already connected, i.e. has only one components Hmin=Hmax. So no dependencies to be handled, just use the optimal order.")
            next
        }

        round <- 1
        while (TRUE) {
            print(paste("############# Round ",round," #################"))
            print("restructuring component")
            isocomps[[i]] <- restructure_components2(isocomps[[i]], bidag_scores, cpp_friendly_scores, opt_sub_orders)
            
            print("restructuring component done")
            if (isocomps[[i]]$no_cycles) break # no cycles to merge.
            round <- round + 1
        }
    }

    # Print and concat the whole isocomp
    print("Concatenating the suborders of each isolated component")
    isocomp_score <- 0
    for (i in seq(1, max(isolated_comp_membership))) {
        # Concatenate subboorder according to the DAG of component dependencies        
        #print(isocomps[[i]])
        print(paste("Concatenating isocomp:",i))
        res <- concat_suborders2(isocomps[[i]]) # returns joined isocomp order and score
        print("concat suborders done")
        # join the sub matrices
        full_adjmat <- matrix(0, p, p)
        for (sc in isocomps[[i]]$subcomp) {
            full_adjmat <- full_adjmat + sc$subadjmat
        }
        #fulldag <- igraph::graph_from_adjacency_matrix(full_adjmat, mode="directed")
        print(res$order)
        print(res$score)
        res$adjmat <- full_adjmat
        isocomp_score <- isocomp_score + res$score

    }
    print("tot score")
    print(isocomp_score)
    return(isocomps)
}

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

component_dependence2 <- function(isocomp, bidag_scores, cpp_friendly_scores) {

    p <- bidag_scores$n
    components <- isocomp$subcomps
    n_components <- length(components)
    component_id1 <- 1
    component_score_sum <- 0

    # Adjacentcy matrix for the component dependence
    comp_dep <- matrix(0, n_components, n_components)
    # list of order for each component
    max_n_particles <- 0
    tot_n_particles <- 0
    tot_order_to_dag_time <- 0
    tot_order_time <- 0

    for (component_id1 in seq(1, n_components)) {
        component1 <- components[[component_id1]] #copy
        initial_suborder <- seq(1, p)[-component1$nodes] #not in optimal order

        # Run order search for the nodes in component_id1, with al the orher nodes as initial suborder, i.e
        # possible parents. Then check to which component the parents belong to.
        #print(paste("******* calculating optimal order for suborder:", component_id1, "********"))
        #print(component1)

        # If we have a new component, calculate the optimal order for it
        if (is.null(component1$score)) { #Dont forget to restor to NULL when components are merged!
            #print("Calculating optimal order for new component")
            start <- proc.time()[1]
            tmp <- optimal_order(cpp_friendly_scores, initial_suborder) # TODO: It seems like i should make this faster, perhaps by now scoring the initial nodes.
            totaltime <- as.numeric(proc.time()[1] - start)
            isocomp$tot_order_time <- isocomp$tot_order_time + totaltime
            isocomp$max_n_particles <- max(isocomp$max_n_particles, max_n_particles)
            isocomp$tot_n_particles <- isocomp$tot_n_particles + tot_n_particles

            start <- proc.time()[1]
            component1_adjmat <- optimal_dag(bidag_scores, cpp_friendly_scores$space, tmp$order)
            #component1_adjmat <- matrix(0, nrow=p, ncol=p)
            totaltime <- as.numeric(proc.time()[1] - start)
            isocomp$tot_order_to_dag_time <- isocomp$tot_order_to_dag_time + totaltime

            # Now, remove iuncoming nodes to all nodes except for those in the component.
            # Then we can paste the all together in the end.
            subadjmat <- component1_adjmat
            # remove all edges not leadning to compnent one
            # This also means that edges from the outsi into the component is alllowed.
            subadjmat[, -component1$nodes] <- 0
            components[[component_id1]]$score <- tmp$suborder_cond_score
            components[[component_id1]]$opt_adjmat <- component1_adjmat
            components[[component_id1]]$subadjmat <- subadjmat
            components[[component_id1]]$suborder <- tmp$suborder

                        
            print("optimal suborder and score:")
            print(components[[component_id1]]$suborder)
            print(components[[component_id1]]$score)
        }

        # In any case we have to evaluate where the potential parents come from,
        # as the component ids might have been renamed.
        # Check if connected component somewhere.
        print("Checking for component dependence")
        for (component_id2 in seq(n_components)) {
            if (isocomp$connected) {
                print("connected component")
                break # Connected component, no need to check for dependencies
            }
            if (component_id1 == component_id2) next
            # Get the edges between the two components
            if(any(component1$opt_adjmat[components[[component_id2]]$nodes, component1$nodes]) > 0) {
                comp_dep[component_id2, component_id1] <- 1 # i.e. it is at least one edge/parent from comp 2 to 1
            }
        }
        print("Checking for component dependence. Done")

    }

    isocomp$comp_dep <- comp_dep
    isocomp$subcomps <- components
    
    isocomp$max_n_particles <- max(isocomp$max_n_particles, max_n_particles)
    isocomp$tot_n_particles <- isocomp$tot_n_particles + tot_n_particles
    isocomp$tot_order_to_dag_time <- isocomp$tot_order_to_dag_time + tot_order_to_dag_time
    isocomp$tot_order_time <- isocomp$tot_order_time + tot_order_time
   
    return(isocomp)
}

# TODO: This should take the membership vector, not the graph
component_dependence <- function(membership, bidag_scores, cpp_friendly_scores, input_opt_suborders=NULL) {
    # Vector of component numbers
    p <- length(membership)

    n_components <- max(membership)
    #print(paste("n_components:", n_components))
    #print(input_opt_suborders)
    # Go through the components of H_min

    component_id1 <- 1
    component_score_sum <- 0

    # Adjacentcy matrix for the component dependence
    comp_dep <- matrix(0, n_components, n_components)
    # list of order for each component
    opt_sub_orders <- list()

    max_n_particles <- 0
    tot_n_particles <- 0
    tot_order_to_dag_time <- 0


    while (component_id1 <= n_components) {

        component1 <- seq(1, p)[membership == component_id1]
        initial_suborder <- seq(1, p)[-component1]

        # Run order search for the nodes in component_id1, with al the orher nodes as initial suborder, i.e
        # possible parents. Then check to which component the parents belong to.
        #print(paste("******* calculating optimal order for suborder:", component_id1, "********"))
        #print(component1)

        # Check if we already have calculations for this components
        component1_exists <- FALSE
        if (!is.null(input_opt_suborders)) {
            # Go through all the components in input_opt_suborders and check if they are component1
            #print("--------- Checking if component exists in input_opt_suborders ----------")
            for (comp_id in seq(1, length(input_opt_suborders))) {
                if (setequal(input_opt_suborders[[comp_id]]$suborder, component1)) {
                    #print("Found component in input_opt_suborders")
                    opt_sub_orders[[component_id1]] <- input_opt_suborders[[comp_id]]
                    component1_exists <- TRUE
                    break
                }
            }
        }
        # If we have a new component, calculate the optimal order for it
        if (component1_exists==FALSE) {
            #print("Calculating optimal order for new component")
            start <- proc.time()[1]
            tmp <- optimal_order(cpp_friendly_scores, initial_suborder)
            totaltime <- as.numeric(proc.time()[1] - start)
            print("Total time to get order")
            print(totaltime)

            #print("done")
            max_n_particles <- max(max_n_particles, tmp$max_n_particles)
            tot_n_particles <- tot_n_particles + tmp$tot_n_particles

            start <- proc.time()[1]
            print("Calculating optimal DAG for new component")
            component1_adjmat <- optimal_dag(bidag_scores, cpp_friendly_scores$space, tmp$order)
            print("done")
            totaltime <- as.numeric(proc.time()[1] - start)
            #print("Total order to dag time //////")
            #print(totaltime)
            tot_order_to_dag_time <- tot_order_to_dag_time + totaltime

            G_opt <- igraph::graph_from_adjacency_matrix(component1_adjmat, mode="directed")
            # Now, remove iuncoming nodes to all nodes except for those in the component.
            # Then we can paste the all together in the end.
            for (i in seq(1, p)) {
                if (membership[i] != component_id1) {
                    component1_adjmat[, i] <- 0
                }
            }
            #print("component_id1:")
            #print(component_id1)
            opt_sub_orders[[component_id1]] <- list("suborder"=tmp$suborder,
                                                    "score"=tmp$suborder_cond_score,
                                                    "dag"=G_opt,
                                                    "subadjmat"=component1_adjmat)

            #print("optimal suborder and score:")
            #print(tmp$suborder)
            #print(tmp$suborder_cond_score)
        } else {
            #print("Component already exists")
        }

        # In any case we have to evaluate where the potential parents come from,
        # as the component ids might have been renamed.
        G_opt <- opt_sub_orders[[component_id1]]$dag # Since it may already be calculated
        vertices1 <- igraph::V(G_opt)[membership == component_id1]$name

        # Now check to which component the parents belong to
        merged_components <- c()

        # Here I only have to look at the components in the same isolated component.
        for (component_id2 in seq(n_components)) {
            if (component_id1 == component_id2) next
            component2 <- seq(1, p)[membership == component_id2]
            # Get the edges between the two components

            # This may be slow
            vertices2 <- igraph::V(G_opt)[membership == component_id2]$name
            between_edges <- igraph::E(G_opt)[vertices1 %<-% vertices2]

            if (length(between_edges) > 0) {
                comp_dep[component_id2, component_id1] <- 1 # i.e. it is at least one edge/parent from comp 2 to 1
                # There are edges from component2 to component1
                #print(paste("Edges from component:", component_id2))
                #print(component2)
                #print(between_edges)
            }
        }
        component_id1 <- component_id1 + 1

    }

    # note that opt_sub_orders will only be valied after the las round when all the merging is done.
    return(list("comp_dep"=comp_dep,
                "opt_sub_orders"=opt_sub_orders,
                "max_n_particles"=max_n_particles,
                "tot_n_particles"=tot_n_particles,
                "tot_order_to_dag_time"=tot_order_to_dag_time))

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

    # Create undirected edges between all nodes in each cycle.
    # Then neigboring ones will be connected, as intended
    # So we can then just take the components.
    n_components <- ncol(adjmat_compdep2)
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

    # go trough and create new subcomps

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

    #print("merged_components_adjmat:")
    #print(merged_components_adjmat)
    #merged_components_graph <- igraph::graph_from_adjacency_matrix(merged_components_adjmat, mode="directed")
    #return(merged_components_graph)
    return(merged_components_adjmat)
}

merged_components_membership <- function(membership, membership_comp){
    # Translates the membership vector of the original graph to the new merged components
    #print("Determin new node memberships in the merged components:")
    p <- length(membership)
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

# TODO: This does not use info aobout isolated components (i.e. H_max). But it should.
# This algorithm works something like:
# 1. Get the dependency graph for the components (which in turn consists of nodes)
# 2. Fin the cycles in the components dependence and merge them into new components.
# 3. Find the dependence between the new, merged components.
# 4. Now we go back to considering the nodes, and update which component they now belong to
# 5. Return all info, and if there was any cycles, if not this component is done.
restructure_components <- function(membership, isolated_comp_membership, bidag_scores, cpp_friendly_scores, input_opt_suborders=NULL) {

    # 1.
    # Get the component dependence graph
    print("Calculating component dependence")
    ret <- component_dependence(membership, bidag_scores, cpp_friendly_scores, input_opt_suborders=input_opt_suborders)
    print("done")
    #print(ret)

    adjmat_compdep <- ret$comp_dep
    G_compdep <- igraph::graph_from_adjacency_matrix(adjmat_compdep, mode="directed")
    print("Component dependence graph before possible merging:")
    print(igraph::E(G_compdep))
    print("and the adjmat")
    print(adjmat_compdep)

    # 2.
    # Merge components that are dependent on each other in a cycle and neigboring cycles.
    # This is for the graph of comopnents. So its the new components membership after merging the cycles.
    #print("Component dependence graph:")
    membership_comp <- merged_neig_cycles(adjmat_compdep)
    print("membership after merging:")
    print(membership_comp)

    # 3.
    # Get dependency graph for the merged components
    merged_components_adjmat <- merged_component_dependencies(adjmat_compdep, membership_comp)
    merged_components_graph <-  igraph::graph_from_adjacency_matrix(merged_components_adjmat, mode="directed")
    # 4.
    # Now, go trough the nodes in the original graph, and check their components
    # and check in the merged components which components they belong to.
    ultimate_membership <- merged_components_membership(membership, membership_comp)
    #print("ultimate_membership:")
    #print(ultimate_membership)

    # Could be done aboe but like this we should be able to verify that the merging is correct.
    # Here we can chag if there where no merged components, i.e. no cycles.
    # By checking if the number of components is the same as the number of merged components.
    if (max(membership_comp) == length(membership_comp)) {
        #print("No cycles, returning the original order")
        return(list("membership"=membership,
                    "merged_components_graph"= G_compdep,
                    "no_cycles"=TRUE,
                    "component_order"=ret,
                    "opt_sub_orders"=ret$opt_sub_orders,
                    "max_n_particles"=ret$max_n_particles,
                    "tot_n_particles"=ret$tot_n_particles,
                    "tot_order_to_dag_time"=ret$tot_order_to_dag_time))
    }

    return(list("membership"=ultimate_membership,
                "merged_components_graph"= merged_components_graph,
                "no_cycles"=FALSE,
                "component_order"=NULL,
                "opt_sub_orders"=ret$opt_sub_orders,
                "max_n_particles"=ret$max_n_particles,
                "tot_n_particles"=ret$tot_n_particles,
                "tot_order_to_dag_time"=ret$tot_order_to_dag_time))
}

# TODO: This does not use info aobout isolated components (i.e. H_max). But it should.
# This algorithm works something like:
# 1. Get the dependency graph for the components (which in turn consists of nodes)
# 2. Fin the cycles in the components dependence and merge them into new components.
# 3. Find the dependence between the new, merged components.
# 4. Now we go back to considering the nodes, and update which component they now belong to
# 5. Return all info, and if there was any cycles, if not this component is done.
restructure_components2 <- function(isocomp, bidag_scores, cpp_friendly_scores, 
                                    input_opt_suborders=NULL) {

    # 1.
    # Get the component dependence graph
    print("Calculating component dependence")
    isocomp <- component_dependence2(isocomp, bidag_scores, cpp_friendly_scores)
    print("done")
    #print(ret)

    adjmat_compdep <- isocomp$comp_dep
    print("compdep adjmat before merging")
    print(adjmat_compdep)
    G_compdep <- igraph::graph_from_adjacency_matrix(adjmat_compdep, mode="directed")
    print("Component dependence graph before possible merging:")
    print(igraph::E(G_compdep))

    # 2.
    # Merge components that are dependent on each other in a cycle and neigboring cycles.
    # This is for the graph of components. So its the new components membership after merging the cycles.
    #print("Component dependence graph:")
    membership_comp <- merged_neig_cycles(adjmat_compdep)
    print("membership after merging:")
    print(membership_comp)

    # go through and create the new components
    # merg some of them anc keep some of them
    n_newcomps <- max(membership_comp)
    #print("printing entire isolate component")
    #print(ret)
    new_subcomps <- list()
    for (i in seq(1, n_newcomps)) {
        print(paste("Component:", i))
        # here we get all the components that below to a new component
        comps_to_merge <- which(membership_comp == i)
        # if only one, dont merge, use the same as before (a copy)
        if (length(comps_to_merge) == 1) {
            print("No need to merge, using: ")
            print(isocomp$subcomps[[comps_to_merge]]$suborder)
            #print(isocomp$subcomps[[comps_to_merge]])
            new_subcomps[[i]] <- isocomp$subcomps[[comps_to_merge]]
        } else { # merge them if more then one
            new_subcomp <- list()
            print("Merging components:")
            print(comps_to_merge)
            new_subcomp_nodes <- c()
            for (comp_no in comps_to_merge){
                new_subcomp_nodes <- c(new_subcomp_nodes, isocomp$subcomps[[comp_no]]$nodes)
            }
            new_subcomp$nodes <- new_subcomp_nodes
            new_subcomp$score <- NULL
            new_subcomp$opt_adjmat <- NULL
            new_subcomp$subadjmat <- NULL
            new_subcomp$suborder <- NULL
            new_subcomps[[i]] <- new_subcomp
        }
    }

    # 3.
    # Get dependency graph for the merged components
    merged_components_adjmat <- merged_component_dependencies(adjmat_compdep, membership_comp)
    merged_components_graph <- igraph::graph_from_adjacency_matrix(merged_components_adjmat, mode="directed")
    # print merged_components_graph
    print("merged_components_graph:")
    print(igraph::E(merged_components_graph))
    print("merged_components_graph done")
    
    # 4. Update the isolate component
    #isocomp <- ret$isocomp
    isocomp$subcomps <- new_subcomps
    #isocomp$comp_dep <- merged_components_graph # should probably be an adjmat
    isocomp$comp_dep <- merged_components_adjmat

  
    if (max(membership_comp) == length(membership_comp)) {
        #print("No cycles, returning the original order")
        isocomp$no_cycles <- TRUE
        
    } else {
        isocomp$no_cycles <- FALSE
    }
    return(isocomp)
}


concat_suborders <- function(ultimate_membership, merged_components_graph, component_order){
    ## Paste the pieaces of suborders acccording to the topological order of merged components.
    sorted_components <- igraph::topo_sort(merged_components_graph, mode= "in")
    #print("sorted_components:")
    #print(sorted_components)

    full_order <- c()
    total_score <- 0
    for (comp_id in sorted_components) {
        suborder <- component_order$opt_sub_orders[[comp_id]]$suborder
        #print(paste("suborder:", comp_id))
        #print(suborder)
        #print("score:")
        #print(component_order$opt_sub_orders[[comp_id]]$score)
        full_order <- c(full_order, suborder) #maybe the other way around
        total_score <- total_score + component_order$opt_sub_orders[[comp_id]]$score
    }
    return(list("order"=full_order, "score"=total_score))
}

concat_suborders2 <- function(isocomp){
    print("in concat2")
    print(isocomp$comp_dep)
    #print(igraph::V(isocomp$comp_dep))
    ## Paste the pieces of suborders acccording to the topological order of merged components.
    G_compdep <- igraph::graph_from_adjacency_matrix(isocomp$comp_dep, mode="directed")
    print(G_compdep)
    sorted_components <- igraph::topo_sort(G_compdep, mode= "in")
    print("sorted_components:")
    print(length(sorted_components))
    print(sorted_components)

    full_order <- c()
    total_score <- 0
    for (comp_id in sorted_components) {
        suborder <- isocomp$subcomps[[comp_id]]$suborder
        print(paste("suborder:", comp_id))
        print(suborder)
        print("score:")
        print(isocomp$subcomps[[comp_id]]$score)
        full_order <- c(full_order, suborder) #maybe the other way around
        total_score <- total_score + isocomp$subcomps[[comp_id]]$score
    }

    return(list("order"=full_order, "score"=total_score))
}
