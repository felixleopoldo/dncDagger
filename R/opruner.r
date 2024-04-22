library("Rcpp")
# library("Jmisc")
# library("pcalg")
# library("ggplot2")
# library("testit")
library(argparser)
library("BiDAG")
Sys.setenv("PKG_CXXFLAGS" = "-fconcepts -Wall -pipe -Wno-unused -pedantic -Wall -L /usr/lib/R/lib -l R -L /usr/local/lib/R/site-library/RInside/lib/ -l RInside -Wl,-rpath,/usr/local/lib/R/site-library/RInside/lib  -I /usr/local/lib/R/site-library/RInside/include/ -I /usr/local/lib/R/site-library/Rcpp/include/ -I /usr/share/R/include/ -std=c++17 -O3")
sourceCpp("include/opruner_right.cpp",  verbose = TRUE)
sourceCpp("include/dnc.cpp",  verbose = TRUE)


optimal_order <- function(cpp_friendly_scores, initial_sub_order){

    # Check if the initial suborder is empty
    #print("Running optimal order")
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
    #isocomps <- list()
    membership <- separable_comp_membership
    for (i in seq(1, max(isolated_comp_membership))) {
        nodes <- which(isolated_comp_membership == i)
        isocomps[[i]] <<- list()
        isocomps[[i]]$nodes <<- nodes
        isocomps[[i]]$subcompids <<- c()
        isocomps[[i]]$subcomps <<- list()
        isocomps[[i]]$connected <<- FALSE # This means that the components Could actually be connected
    }

    # Find the sub components in H_min
    #print("Finding super components for:")
    for (j in seq(1, max(membership))) {
         subnodes <- which(membership == j)
 
        if(length(unique(isolated_comp_membership[subnodes])) == 1) {
            isocompid <- isolated_comp_membership[subnodes][1]
            isocomps[[isocompid]]$subcompids <<- c(j, isocomps[[isocompid]]$subcompids)

            k <- length(isocomps[[isocompid]]$subcomps)

            isocomps[[isocompid]]$subcomps[[k+1]] <<- list("nodes" = subnodes,
                                                          "score" = NULL,
                                                          "opt_adjmat" = NULL,
                                                          "subadjmat" = NULL,
                                                          "suborder" = NULL)
            isocomps[[isocompid]]$max_n_particles <<- 0
            isocomps[[isocompid]]$tot_n_particles <<- 0
            isocomps[[isocompid]]$tot_order_to_dag_time <<- 0
            isocomps[[isocompid]]$tot_order_time <<- 0

            # Check if it is the same, ie the component is isolated and connected.
            if(length(subnodes) == length(isocomps[[isocompid]]$nodes)){
                isocomps[[isocompid]]$connected <<- TRUE
            }
        }
    }
    # print number of subcomponents of each isolated component
    for (i in seq(1, max(isolated_comp_membership))) {
        print(paste("Isolated component:", i, " has this number of subcomponents"))
        print(length(isocomps[[i]]$subcompids))
        # print(isocomps[[i]]$subcompids)
        # # print the nodes in each subcomponent
        # for (j in seq(1, length(isocomps[[i]]$subcomps))) {
        #     print(paste("Subcomponent:", j))
        #     print(isocomps[[i]]$subcomps[[j]]$nodes)
        # }

    }


    isocomps$max_n_particles <<- 0
    isocomps$tot_n_particles <<- 0
    isocomps$tot_order_to_dag_time <<- 0
    isocomps$tot_order_time <<- 0

    #return(isocomps)
}


dnc2 <- function() {
    #bidag_scores <- cpp_friendly_scores$bidag_scores
    # This is now working if aliases is empty
    starttot <- proc.time()
    
    start <- proc.time()
    aliases <- lapply(cpp_friendly_scores$aliases, function(a) a + 1)
    var_labels <- cpp_friendly_scores$labels
    #print("Creating diff matrices")
    #print("Scoretables")
    #print(cpp_friendly_scores$scoretable)
    diff_matrices <- get_diff_matrices(cpp_friendly_scores$rowmaps, cpp_friendly_scores$scoretable, aliases, var_labels)
    H_min <- diff_matrices$H_min
    H_max <- diff_matrices$H_max
    p <- nrow(H_min)
    H_min_adj <- (H_min> 0)*1
    H_max_adj <- (H_max> 0)*1
    G_H_min <- igraph::graph_from_adjacency_matrix(H_min_adj, mode="directed")
    G_H_max <- igraph::graph_from_adjacency_matrix(H_max_adj, mode="directed")
  
    # Components of H_min. Possible component in G.
    print("Components of H_min:")
    print(igraph::components(G_H_min, mode="weak"))

    # Components of H_max.
    # One component of this may contain several components of G.
    # We call them isolated components.
    print("Components of H_max:")
    print(igraph::components(G_H_max, mode="weak"))

    # Run order search on each component from H_min and see which components in Hmin
    # That the parents belong to

    membership <- igraph::components(G_H_min)$membership
    isolated_comp_membership <- igraph::components(G_H_max)$membership
    
    totaltime <- proc.time() - start
    print("Total time for creating diff matrices:")
    print(totaltime)
    round <- 1
    opt_sub_orders <- NULL

    # Order the components in H_min after the isolated components (induced by H_max).
    start <- proc.time()
    # THIS IS A GLOBAL VARIABLE
    #isocomps <<- list()
    structure_components(isolated_comp_membership, membership)
    #isocomps <<- structure_components(isolated_comp_membership, membership)
    totaltime <- proc.time() - start
    
    #print("Total time for finding super components:")
    #print(totaltime)

    # We initialize the process by computing the component dependencies once for each
    # isolated component.

    print("Going through the isolated components and get the sub component dependencies.")
    start <- proc.time()
    for (i in seq(1, max(isolated_comp_membership))) {
        #print(paste("Isolated component:", i))
        component_dependence2(i)#, bidag_scores, cpp_friendly_scores)
        #print("Isolated component done")
        # Update isocomps costs(time and particle counts)
        isocomps$max_n_particles <<- max(isocomps$max_n_particles, isocomps[[i]]$max_n_particles)
        isocomps$tot_n_particles <<- isocomps$tot_n_particles + isocomps[[i]]$tot_n_particles
        isocomps$tot_order_to_dag_time <<- isocomps$tot_order_to_dag_time + isocomps[[i]]$tot_order_to_dag_time
        isocomps$tot_order_time <<- isocomps$tot_order_time + isocomps[[i]]$tot_order_time
    }
    totaltime <- proc.time() - start
    #print("--------------------  The isocomps after initial round")
    #print(isocomps)
    print("*********** First round setup of Isolated components done.")
    # The divide and conquer procedure should, operate individually ion eas isolated component.
    print("************** Starting divide and conquer. Going through the isolated components and get the internal component dependencies.")
    #print("number of isocomps")
    #print(max(isolated_comp_membership))
    for (i in seq(1, max(isolated_comp_membership))) {

        #print(paste("*********** Treating isolated component", i))
        # print nodes in the subcomponents
        start <- proc.time() - isocomps[[i]]$tot_order_to_dag_time
        if (isocomps[[i]]$connected) {
            #print("Its already connected, i.e. has only one components Hmin=Hmax. So no dependencies to be handled, just use the optimal order.")
            next
        }

        round <- 1
        while (TRUE) {
            #print(paste("############# Merging round ",round," #################"))
            #print("restructuring component")
            #isocomps[[i]] <- restructure_components2(isocomps[[i]], bidag_scores, cpp_friendly_scores, opt_sub_orders)
            #restructure_components2(i, bidag_scores, cpp_friendly_scores, opt_sub_orders)
            start <- proc.time()
            restructure_components2(i)
            # print time minus the time for order to dag
            #print("Total time for restructure component:")
            #print(proc.time() - start - isocomps[[i]]$tot_order_to_dag_time)
            
            #print("restructuring component done")
            if (isocomps[[i]]$no_cycles) break # no cycles to merge.
            round <- round + 1
        }
        #print(">>>> Time for isolated component - ordertodag")
        #print(proc.time() - start - isocomps[[i]]$tot_order_to_dag_time)
    }

    # Print and concat the whole isocomp
    #print("Concatenating the suborders of each isolated component")
    isocomps$order <<- c()
    isocomp_score <- 0
    for (i in seq(1, max(isolated_comp_membership))) {
        # Concatenate subboorder according to the DAG of component dependencies
        ##print(paste("Concatenating isocomp:",i))
        res <- concat_suborders2(isocomps[[i]]) # returns joined isocomp order and score
        #print("concat suborders done")
        # join the sub matrices
        full_adjmat <- matrix(0, p, p)
        for (sc in isocomps[[i]]$subcomp) {
            full_adjmat <- full_adjmat + sc$subadjmat
            # Calculating the costs for all the components together
            isocomps$tot_n_particles <<- isocomps$tot_n_particles + sc$tot_n_particles
            isocomps$tot_order_to_dag_time <<- isocomps$tot_order_to_dag_time + sc$tot_order_to_dag_time
            isocomps$max_n_particles <<- max(isocomps$max_n_particles, sc$max_n_particles)
        }
        
        res$adjmat <- full_adjmat
        isocomp_score <- isocomp_score + res$score
        isocomps$order <<- c(isocomps$order, res$order)
       
    }

    #print("tot score")
    #print(isocomp_score)

    isocomps$score <<- isocomp_score
    # print("Total time for dnc:")
    # print(proc.time() - starttot)
    # print("Total time for order to dag:")
    # print(isocomps$tot_order_to_dag_time)
    # print("subtracted time")
    # print(proc.time() - starttot - isocomps$tot_order_to_dag_time)
    # #return(isocomps)
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

component_dependence2 <- function(isocomp_id) {
#component_dependence2 <- function(isocomp_id, bidag_scores, cpp_friendly_scores) {
#component_dependence2 <- function(isocomp, bidag_scores, cpp_friendly_scores) {

    p <- bidag_scores$n
    #components <- isocomp$subcomps
    n_components <- length(isocomps[[isocomp_id]]$subcomps)
    #print(paste("Number of sub components:", n_components))
    component_id1 <- 1
    component_score_sum <- 0

    # Adjacency matrix for the component dependence.
    comp_dep <- matrix(0, n_components, n_components)

    for (component_id1 in seq(1, n_components)) {
        # Run order search for the nodes in component_id1, with all the other 
        # nodes as initial suborder, i.e possible parents. Then check to which 
        # component the parents belong to.         
        # If we have a new component, calculate the optimal order for it.
        if (is.null(isocomps[[isocomp_id]]$subcomps[[component_id1]]$score)) { # Dont forget to restor to NULL when components are merged!
            #print("Calculating optimal order for new component")
            
            # TODO: It seems like i should make this faster, perhaps by now scoring the initial nodes.
            initial_suborder <- seq(1, p)[-isocomps[[isocomp_id]]$subcomps[[component_id1]]$nodes] #not in optimal order
            start1 <- proc.time()
            tmp <- optimal_order(cpp_friendly_scores, initial_suborder) 
            totaltime1 <- proc.time() - start1

            start2 <- proc.time()
            #component1_adjmat <- optimal_dag(bidag_scores, cpp_friendly_scores$space, tmp$order) # component1_adjmat <- matrix(0, nrow=p, ncol=p)
            component1_adjmat <- r_order_to_dag(cpp_friendly_scores, tmp$order)
            totaltime2 <- proc.time() - start2

            # Store the "cost" for finding the suborder and DAG
            isocomps[[isocomp_id]]$tot_order_time <<- isocomps[[isocomp_id]]$tot_order_time + totaltime1
            isocomps[[isocomp_id]]$max_n_particles <<- max(isocomps[[isocomp_id]]$max_n_particles, tmp$max_n_particles)
            isocomps[[isocomp_id]]$tot_n_particles <<- isocomps[[isocomp_id]]$tot_n_particles + tmp$tot_n_particles
            isocomps[[isocomp_id]]$tot_order_to_dag_time <<- isocomps[[isocomp_id]]$tot_order_to_dag_time + totaltime2

        
            # Remove all edges not leadning to compnent one,
            # this also means that edges from the outside into the component is alllowed.
            isocomps[[isocomp_id]]$subcomps[[component_id1]]$score <<- tmp$suborder_cond_score
            isocomps[[isocomp_id]]$subcomps[[component_id1]]$opt_adjmat <<- component1_adjmat # O(p^2)?            
            isocomps[[isocomp_id]]$subcomps[[component_id1]]$subadjmat <<- component1_adjmat
            isocomps[[isocomp_id]]$subcomps[[component_id1]]$subadjmat[, -isocomps[[isocomp_id]]$subcomps[[component_id1]]$nodes] <<- 0
            isocomps[[isocomp_id]]$subcomps[[component_id1]]$suborder <<- tmp$suborder
        }

        # In any case we have to evaluate where the potential parents come from,
        # as the component ids might have been renamed.
        # Check if connected component somewhere.
        #print("Checking for component dependence")
        for (component_id2 in seq(n_components)) {
            if (isocomps[[isocomp_id]]$connected) {
                #print("connected component")
                break # Connected component, no need to check for dependencies
            }
            if (component_id1 == component_id2) next
            # Get the edges between the two components. This might take some time.
            if(any(isocomps[[isocomp_id]]$subcomps[[component_id1]]$opt_adjmat[isocomps[[isocomp_id]]$subcomps[[component_id2]]$nodes, 
                                                                               isocomps[[isocomp_id]]$subcomps[[component_id1]]$nodes]) > 0) {
                comp_dep[component_id2, component_id1] <- 1 # i.e. it is at least one edge/parent from comp 2 to 1
            }
        }
        #print("Checking for component dependence. Done")
        
    }

    isocomps[[isocomp_id]]$comp_dep <<- comp_dep
    # print the component dependence
    print("Component dependence:")
    print(comp_dep)
    #isocomp$subcomps <- components
    #return(isocomp)
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
            start <- proc.time()
            tmp <- optimal_order(cpp_friendly_scores, initial_suborder)
            totaltime <- proc.time() - start
            print("Total time to get order")
            print(totaltime)

            #print("done")
            max_n_particles <- max(max_n_particles, tmp$max_n_particles)
            tot_n_particles <- tot_n_particles + tmp$tot_n_particles

            start <- proc.time()
            print("Calculating optimal DAG for new component")
            component1_adjmat <- optimal_dag(bidag_scores, cpp_friendly_scores$space, tmp$order)
            print("done")
            totaltime <- proc.time() - start
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


# This algorithm works something like:
# 1. Get the dependency graph for the components (which in turn consists of nodes)
# 2. Find the cycles in the components dependence and merge them into new components.
# 3. Find the dependence between the new, merged components.
# 4. Now we go back to considering the nodes, and update which component they now belong to
# 5. Return all info, and if there was any cycles, if not this component is done.
# restructure_components2 <- function(isocomp, bidag_scores, cpp_friendly_scores,
#                                     input_opt_suborders=NULL) {
restructure_components2 <- function(isocomp_id#, #bidag_scores, cpp_friendly_scores,
                                    ) {
    # 1.
    # Get the component dependence graph
    #isocomp <- component_dependence2(isocomp, bidag_scores, cpp_friendly_scores)
    start <- proc.time()
    tot_order_to_dag_time <- isocomps[[isocomp_id]]$tot_order_to_dag_time
    component_dependence2(isocomp_id)#, bidag_scores, cpp_friendly_scores)
    tot_order_to_dag_time <- isocomps[[isocomp_id]]$tot_order_to_dag_time - tot_order_to_dag_time 
    adjmat_compdep <- isocomps[[isocomp_id]]$comp_dep
    #print("Time for finding component dependence")
    #print(proc.time() - start - tot_order_to_dag_time)

    # 2.
    # Merge components that are dependent on each other in a cycle and neigboring cycles.
    # This is for the graph of components. So its the new components membership after merging the cycles.
    
    start <- proc.time()
    membership_comp <- merged_neig_cycles(adjmat_compdep) # slow
    #print("Time for finding new components by merging cycles")
    #print(proc.time() - start)

    start <- proc.time()
    # Go through and create the new components
    # merge some of them and keep some of them.
    n_newcomps <- max(membership_comp)
    new_subcomps <- list()
    for (i in seq(1, n_newcomps)) {
        #print(paste("Internal component:", i))
        # Here we get all the components that below to a new component
        comps_to_merge <- which(membership_comp == i)
        # If only one, dont merge, use the same as before (a copy)
        if (length(comps_to_merge) == 1) {
            #print("No need to merge, using: ")
            new_subcomps[[i]] <- isocomps[[isocomp_id]]$subcomps[[comps_to_merge]] # this may me slow
        } else { 
            # Merge them if more then one
            new_subcomp <- list()
            #print("Merging components:")
            new_subcomp_nodes <- c()
            for (comp_no in comps_to_merge){                
                new_subcomp_nodes <- c(new_subcomp_nodes, isocomps[[isocomp_id]]$subcomps[[comp_no]]$nodes)
            }
            new_subcomp$nodes <- new_subcomp_nodes
            new_subcomp$score <- NULL
            new_subcomp$opt_adjmat <- NULL
            new_subcomp$subadjmat <- NULL
            new_subcomp$suborder <- NULL
            new_subcomps[[i]] <- new_subcomp
        }
    }
    # print("Time for maybe creating new sub components")
    # print(proc.time() - start)

    # 3.
    # Get dependency graph for the merged components
    start <- proc.time()
    merged_components_adjmat <- merged_component_dependencies(adjmat_compdep, membership_comp)
    # print("Time for finding new component dependencies")
    # print(proc.time() - start)

    # 4. Update the isolate component
    start <- proc.time()
    isocomps[[isocomp_id]]$subcomps <<- new_subcomps
    isocomps[[isocomp_id]]$comp_dep <<- merged_components_adjmat

    if (max(membership_comp) == length(membership_comp)) {
        #print("No cycles, returning the original order")
        isocomps[[isocomp_id]]$no_cycles <<- TRUE
    } else {
        isocomps[[isocomp_id]]$no_cycles <<- FALSE
    }


    # totaltime <- proc.time() - start
    # print("Time for updating isolated component:")
    # print(totaltime)

    #return(isocomp)
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
    ## Paste the pieces of suborders acccording to the topological order of merged components.
    G_compdep <- igraph::graph_from_adjacency_matrix(isocomp$comp_dep, mode="directed")
    sorted_components <- igraph::topo_sort(G_compdep, mode= "in")

    full_order <- c()
    total_score <- 0
    for (comp_id in sorted_components) {
        suborder <- isocomp$subcomps[[comp_id]]$suborder
        full_order <- c(full_order, suborder)
        total_score <- total_score + isocomp$subcomps[[comp_id]]$score
    }

    return(list("order"=full_order, "score"=total_score))
}
