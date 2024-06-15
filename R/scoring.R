rm(list = ls())
library("Rcpp")
library("Jmisc")

insertSource("R/iterativeMCMC.R", package = "BiDAG")

source("R/export_to_gobnilp_score_tables.R")

get_scores <- function(filename,  scoretype = c("bge", "bde", "bdecat"),
                      bgepar = list(am = 1, aw = NULL),
                      bdepar = list(chi = 0.5, edgepf = 2), 
                      plus1it=NULL, 
                      iterations=NULL,
                      gobnilp_scores_filename=NULL
                      ) {
    MAP <- TRUE

    data <- read.csv(filename, check.names = FALSE)
    if (scoretype =="bge") {
        myscore <- scoreparameters(scoretype = scoretype, data, bgepar = bgepar)
    } else if (scoretype == "bde") {
        myscore <- scoreparameters(scoretype = scoretype, data[-1, ], bdepar = bdepar)
    }
    ret <- get_plus1_score_essentials_for_cpp(myscore, 
                                              plus1it=plus1it, 
                                              iterations=iterations, 
                                              gobnilp_scores_filename=gobnilp_scores_filename)

    aliases <- lapply(ret$aliases, function(a) a + 1)
    diff_matrices <- get_diff_matrices(ret$rowmaps, ret$scoretable, aliases, labels(data)[[2]])

    ret$H_min <- diff_matrices$H_min
    ret$H_max <- diff_matrices$H_max
    ret$H_min_adj <- (ret$H_min> 0)*1
    ret$H_max_adj <- (ret$H_max> 0)*1
    ret$bidag_scores <- myscore
    ret$labels <- labels(data)[[2]]
    return(ret)
}

get_diff_matrices <- function(rowmaps, scoretable, aliases, var_labels){

    nvars <- length(rowmaps)
    H_min = matrix(, nrow = nvars, ncol = nvars)
    H_max = matrix(, nrow = nvars, ncol = nvars)
    colnames(H_min) <- var_labels
    colnames(H_max) <- var_labels
    rownames(H_min) <- var_labels
    rownames(H_max) <- var_labels

    for (i in seq(nvars)) {
        var <- rowmaps[[i]]
        # total number of possible parents
        n_pos_parents <- length(aliases[[i]])# sqrt(length(var$forward))

        # For the possible parents, i.e. not plus1 parents
        # We exclude each parent in turn and see how the score changes
        for (parent_ind in seq(n_pos_parents)){
            # if no possible parents, skip
            if(n_pos_parents == 0) next # since seq is weird
            parent <- aliases[[i]][[parent_ind]] # parent to exclude
            for (hash in var$forward){

                # Check if the parent is in the hash
                check <- (hash-1) %% 2^(parent_ind)
                if (check < 2^(parent_ind-1)){
                    # has parent parent
                    #print(paste0("Has parent ind: ", parent_ind))

                    # Compute the hash with the parent excluded
                    hash_with_parent <- hash + 2^(parent_ind-1)
                    # Using the no plus1 score table, i.e. index 1
                    score_diff <- scoretable[[i]][[1]][var$backwards[[hash_with_parent]]] - scoretable[[i]][[1]][var$backwards[[hash]]]
                    # Update the H matrices
                    if (is.na(H_max[i, parent])){
                        H_max[i,parent] <- score_diff
                    } else {
                        H_max[i, parent] <- max(H_max[i, parent], score_diff)
                    }

                    if (is.na(H_min[i, parent])){
                        H_min[i,parent] <- score_diff
                    } else {
                        H_min[i, parent] <- min(H_min[i, parent], score_diff)
                    }
                 }
            }
        }

        ## For the plus1 parents
        # plus1parents are those parents that are not in the aliases
        plus1parents <- c()
        plus1parent_inds <- c()

        j <- 1
        for (label in var_labels){
            if (j == i) {  # skip itself
                j <- j + 1
                next
            }
            if (!(label %in% labels(aliases[[i]]))){
                # label not in aliases so it is a plus1 parent
                plus1parents <- c(plus1parents, label)
                plus1parent_inds <- c(plus1parent_inds, j)
            }
            j <- j + 1
        }

        for(j in seq(1, length(plus1parents))){
            score_diffs <- scoretable[[i]][[j+1]] - scoretable[[i]][[1]] # the first one is the no plus1 score table. Subtracting all at once.
            H_max[i, plus1parent_inds[j]] <- max(score_diffs)
            H_min[i, plus1parent_inds[j]] <- min(score_diffs)
        }
    }
    # print("H_max:")
    # print((H_max > 0) * 1)
    # print("H_min:")
    # print((H_min > 0) * 1)
    return(list(H_min = H_min, H_max = H_max))
}

get_plus1_score_essentials_for_cpp <- function(myscore, plus1it=NULL, iterations=NULL, gobnilp_scores_filename=NULL) {
  
  res <- iterativeMCMC(myscore, chainout = TRUE, scoreout = TRUE, MAP = TRUE,
                       plus1it=plus1it, iterations=iterations, verbose=TRUE) #this is bidag version 2.0.0


  print("score from MCMC:")
  print(res$result$score)

    if (!is.null(gobnilp_scores_filename)) {

        #print("Writing gobnilp scoretables")
        gobnscores<-write_gobnilp_scores(res$result$scoretable$tables, 
                                res$result$scoretable$adjacency, 
                                gobnilp_scores_filename)
        #cat(gobnscores)
        write(gobnscores, file = gobnilp_scores_filename)
        #print("done")
    }

  ret <- list()
  ret$parenttable <- lapply(res$ptab$parenttable, function(a) {
    df <- data.frame(a)
    df[is.na(df)] <- 0
    df <- df - 1
    m <- as.matrix(df)
    m <- as.matrix(df, row.names = 0, col.names = 0)
    rownames(m) <- colnames(m) <- NULL
    return(m)
  })

  ret$aliases <- lapply(res$ptab$aliases, function(a) a - 1)
  ret$numparents <- res$ptab$numparents
  ret$rowmaps_backwards <- lapply(res$rowmaps, function(a) a$backwards - 1)
  ret$rowmaps_forward <- lapply(res$rowmaps, function(a) a$forward - 1)

  ret$plus1listsparents <- lapply(res$plus1lists$parents, function(a) a - 1)
  ret$scoretable <- res$result$scoretable$table
  ret$bannedscore <- res$bannedscore
  ret$MAP <- TRUE
  ret$space <- res$result$endspace
  ret$rowmaps <- res$rowmaps
  ret$maxmatrices <- res$maxmatrices
  ret$maxmatrix <- res$maxmatrices$maxmatrix
  ret$maxrow <- lapply(res$maxmatrices$maxrow, function(a) a - 1)

  #ret$bidag_scores # put it in here too

  #print("max matrices/banned scores:")
  #print(ret$bannedscore)
  return(ret)
}

