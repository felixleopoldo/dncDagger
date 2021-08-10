#'This function implements the order MCMC algorithm for the structure learning of Bayesian networks. This function can be used
#'for MAP discovery and for sampling from the posterior distribution of DAGs given the data.
#'Due to the superexponential size of the search space as the number of nodes increases, the 
#'MCMC search is performed on a reduced search space.
#'By default the search space is limited to the skeleton found through the PC algorithm by means of conditional independence tests 
#'(using the functions \code{\link[pcalg]{skeleton}} and \code{\link[pcalg]{pc}} from the `pcalg' package [Kalisch et al, 2012]).
#'It is also possible to define an arbitrary search space by inputting an adjacency matrix, for example estimated by partial correlations or other network algorithms.
#'Also implemented is the possibility to expand the default or input search space, by allowing each node in the network to have one additional parent.  This offers improvements in the learning and sampling of Bayesian networks. 
#' @param scorepar an object of class \code{scoreparameters}, containing the data and score parameters, see constructor function \code{\link{scoreparameters}}
#' @param MAP logical, if TRUE (default) the search targets the MAP DAG (a DAG with maximum score),
#' if FALSE at each MCMC step a DAG is sampled from the order proportionally to its score
#' @param plus1 logical, if TRUE (default) the search is performed on the extended search space
#' @param chainout logical, if TRUE the saved MCMC steps are returned, TRUE by default
#' @param scoreout logical, if TRUE the search space and score tables are returned, FALSE by default
#' @param moveprobs a numerical vector of 4 values in \code{\{0,1\}} corresponding to the probabilities of the following MCMC moves in the order space
#' \itemize{
#' \item exchanging 2 random nodes in the order
#' \item exchanging 2 adjacent nodes in the order
#' \item placing a single node elsewhere in the order
#' \item staying still
#' }
#' @param iterations integer, the number of MCMC steps, the default value is \eqn{5n^{2}\log{n}}
#' @param stepsave integer, thinning interval for the MCMC chain, indicating the number of steps between two output iterations, the default is \code{iterations/1000}
#' @param alpha numerical significance value in \code{\{0,1\}} for the conditional independence tests at the PC algorithm stage (by default \eqn{0.4} for \eqn{n<50}, \eqn{20/n} for \eqn{n>50})
#' @param gamma tuning parameter which transforms the score by raising it to this power, 1 by default
#' @param cpdag logical, if TRUE the CPDAG returned by the PC algorithm will be used as the search
#'space, if FALSE (default) the full undirected skeleton will be used as the search space
#' @param hardlimit integer, limit on the size of parent sets in the search space; by default 14 when MAP=TRUE and 20 when MAP=FALSE
#' @param verbose logical, if TRUE messages about the algorithm's progress will be printed, FALSE by default
#' @param startspace (optional) a square matrix, of dimensions equal to the number of nodes, which defines the search space for the order MCMC in the form of an adjacency matrix. If NULL, the skeleton obtained from the PC-algorithm will be used. If \code{startspace[i,j]} equals to 1 (0) it means that the edge from node \code{i} to node \code{j} is included (excluded) from the search space. To include an edge in both directions, both \code{startspace[i,j]} and \code{startspace[j,i]} should be 1.
#' @param blacklist (optional) a square matrix, of dimensions equal to the number of nodes, which defines edges to exclude from the search space. If \code{blacklist[i,j]} equals to 1 it means that the edge from node \code{i} to node \code{j} is excluded from the search space.
#' @param scoretable (optional) object of class \code{MCMCspace} containing list of score tables calculated for example by the last iteration of the function \code{iterativeMCMC}. When not NULL, parameter \code{startspace} is ignored.
#' @param startorder (optional) integer vector of length n, which will be used as the starting order in the MCMC algorithm, the default order is random
#' @return Object of class \code{MCMCres}, which contains log-score trace of sampled DAGs as well 
#' as adjacency matrix of the maximum scoring DAG, its score and the order score. See \code{\link{MCMCres}}.
#' The output can optionally include DAGs sampled in MCMC iterations and the score tables. Optional output is regulated by the parameters \code{chainout} and \code{scoreout}.
#'@references Friedman N and Koller D (2003). A Bayesian approach to structure discovery in bayesian networks. Machine Learning 50, 95-125.
#'@references Kalisch M, Maechler M, Colombo D, Maathuis M and Buehlmann P (2012). Causal inference using graphical models with the R package pcalg. Journal of Statistical Software 47, 1-26.
#'@references Geiger D and Heckerman D (2002). Parameter priors for directed acyclic graphical models and the characterization of several probability distributions. The Annals of Statistics 30, 1412-1440.
#'@references Kuipers J, Moffa G and Heckerman D (2014). Addendum on the scoring of Gaussian acyclic graphical models. The Annals of Statistics 42, 1689-1691.
#'@references Spirtes P, Glymour C and Scheines R (2000). Causation, Prediction, and Search, 2nd edition. The MIT Press.
#'@examples
#'\dontrun{
#'#find a MAP DAG with search space defined by PC and plus1 neighbourhood
#'Bostonscore<-scoreparameters("bge",Boston)
#'#estimate MAP DAG
#'orderMAPfit<-orderMCMC(Bostonscore)
#'summary(orderMAPfit)
#'#sample DAGs from the posterior distribution
#'ordersamplefit<-orderMCMC(Bostonscore,MAP=FALSE,chainout=TRUE)
#'plot(ordersamplefit)
#'}
#'@author Polina Suter, Jack Kuipers, the code partly derived from the order MCMC implementation from Kuipers J, Moffa G (2017) <doi:10.1080/01621459.2015.1133426>
#'@export
source("R/other.R")
source("R/usrscorefns.R")
source("R/spacefns.R")
library("pcalg")

orderMCMCFelix <- function(scorepar, MAP = TRUE, plus1 = TRUE, chainout = FALSE, scoreout = FALSE, moveprobs = NULL,
                    iterations = NULL, stepsave = NULL, alpha = 0.05, cpdag = FALSE, gamma = 1,
                    hardlimit = ifelse(plus1, 14, 20), verbose = FALSE,
                    startspace = NULL, blacklist = NULL, startorder = NULL, scoretable = NULL) {
  if (is.null(moveprobs)) {
    prob1 <- 99
    if (scorepar$nsmall > 3) { prob1 <- round(6 * 99 * scorepar$nsmall / (scorepar$nsmall ^ 2 + 10 * scorepar$nsmall - 24)) }
    prob1 <- prob1 / 100
    moveprobs <- c(prob1, 0.99 - prob1, 0.01)
    moveprobs <- moveprobs / sum(moveprobs)
    moveprobs <- c(moveprobs[c(1, 2)], 0, moveprobs[3])
  }
  if (is.null(iterations)) {
    if (scorepar$nsmall < 26) {
      iterations <- 30000
    } else {
      iterations <- (5 * scorepar$nsmall * scorepar$nsmall * log(scorepar$nsmall)) - (5 * scorepar$nsmall * scorepar$nsmall * log(scorepar$nsmall)) %% 1000
    }
  }
  if (is.null(stepsave)) {
    stepsave <- floor(iterations / 1000)
  }

  ordercheck <- checkstartorder(startorder, varnames = scorepar$labels.short, mainnodes = scorepar$mainnodes,
                              bgnodes = scorepar$static, DBN = scorepar$DBN, split = scorepar$split)

  if (ordercheck$errorflag) {
    stop(ordercheck$message)
  } else {
    startorder <- ordercheck$order
  }

  result <- orderMCMCmainFelix(param = scorepar, iterations, stepsave, startorder = startorder,
                    moveprobs = moveprobs, alpha = alpha, cpdag = cpdag, scoretable = scoretable,
                    plus1 = plus1, MAP = MAP, chainout = chainout, scoreout = scoreout,
                    startspace = startspace, blacklist = blacklist, gamma = gamma, verbose = verbose,
                    hardlimit = hardlimit)



  if (plus1) {
    result$info$algo <- "plus1 order MCMC"
  } else {
    result$info$algo <- "base order MCMC"
  }

  if (is.null(startspace)) {
    result$info$spacealgo <- "PC"
  } else {
    result$info$spacealgo <- "user defined matrix"
  }
  result$info$iterations <- iterations
  result$info$samplesteps <- length(result$trace)
  if (MAP) {
    result$info$sampletype <- "MAP"
  } else {
    result$info$sampletype <- "sample"
  }

  attr(result, "class") <- "MCMCres"

  return(result)

}

orderMCMCmainFelix <- function(param, iterations, stepsave, MAP = TRUE, posterior = 0.5,
                        startorder = c(1:n), moveprobs, plus1 = FALSE, chainout = TRUE,
                        scoreout = FALSE, startspace = NULL, blacklist = NULL, gamma = 1, verbose = FALSE, alpha = NULL,
                        hardlimit = ifelse(plus1, 15, 22),
                        cpdag = FALSE, addspace = NULL, scoretable = NULL) {
  result <- list()
  maxobj <- list()
  MCMCtraces <- list()

  n <- param$n
  nsmall <- param$nsmall
  matsize <- ifelse(param$DBN, n + nsmall, n)

  #defining startorder and updatenodes 
  if (!param$DBN) {

    if (param$bgn != 0) {
      updatenodes <- c(1:n)[-param$bgnodes]
    } else {
      updatenodes <- c(1:n)
    }

  } else {
    #for DBNs startorder is defined in main.R
    updatenodes <- c(1:nsmall)
  }
  maxorder <- startorder

  #creating blacklist objects
  if (is.null(blacklist)) {
    blacklist <- matrix(0, nrow = matsize, ncol = matsize)
  }
  diag(blacklist) <- 1
  if (!is.null(param$bgnodes)) {
    for (i in param$bgnodes) {
      blacklist[, i] <- 1
    }
  }

  #defining startskel
  if (!is.null(scoretable)) {
    startskel <- scoretable$adjacency
    blacklist <- scoretable$blacklist
    scoretable <- scoretable$tables
  } else {
    if (is.null(startspace)) {
      startspace <- definestartspace(alpha, param, cpdag = cpdag, algo = "pc")
    }
    startskeleton <- 1 * (startspace & !blacklist)
    if (!is.null(addspace)) {
      startskel <- 1 * ((addspace | startskeleton) & !blacklist)
    } else { startskel <- startskeleton }
  }

  blacklistparents <- list()
  for (i in 1:matsize) {
    blacklistparents[[i]] <- which(blacklist[, i] == 1)
  }

  if (verbose) {
    cat(paste("maximum parent set size is", max(apply(startskel, 2, sum))), "\n")
  }
  if (max(apply(startskel, 2, sum)) > hardlimit) {
    stop("the size of maximal parent set is higher that the hardlimit; redifine the search space or increase the hardlimit!")
  }

  tablestart <- Sys.time()
  #computing score tables
  ptab <- listpossibleparents.PC.aliases(startskel, isgraphNEL = FALSE, n, updatenodes)

  if (verbose) {
    cat("skeleton ready \n")
    flush.console()
  }

    parenttable <- ptab$parenttable # basic parenttable without plus1 lists
    aliases <- ptab$aliases #aliases for each node since all nodes in parent tables are done as 1,2,3,4... not real parent names
    numberofparentsvec <- ptab$numberofparentsvec
    numparents <- ptab$numparents
    plus1lists <- PLUS1(matsize, aliases, updatenodes, blacklistparents)
    rowmaps <- parentsmapping(parenttable, numberofparentsvec, n, updatenodes)
    if (is.null(scoretable)) {
      scoretable <- scorepossibleparents.PLUS1(parenttable, plus1lists, n, param, updatenodes,
                                           rowmaps, numparents, numberofparentsvec)
    }
    posetparenttable <- poset(parenttable, numberofparentsvec, rowmaps, n, updatenodes)

    if (MAP == TRUE) {
      maxmatrices <- posetscoremax(posetparenttable, scoretable, numberofparentsvec,
                               rowmaps, n, plus1lists = plus1lists, updatenodes)
    } else {
      bannedscore <- poset.scores(posetparenttable, scoretable, ptab$numberofparentsvec, rowmaps,
                              n, plus1lists = plus1lists, ptab$numparents, updatenodes)
    }

    if (verbose) {
      cat(paste("score tables calculated, MCMC plus1 starts"), "\n")
      flush.console()
    }

    tableend <- Sys.time()

    MCMCresult <- orderMCMCplus1Felix(n, nsmall, startorder, iterations, stepsave, moveprobs, parenttable,
                                 scoretable, aliases, numparents, rowmaps, plus1lists,
                                 bannedscore, numberofparentsvec, gamma = gamma,
                                 bgnodes = param$bgnodes, matsize = matsize)

    mcmcend <- Sys.time()

  #defining result object
  if (chainout) {

      MCMCtraces$incidence <- lapply(MCMCresult$incidence, function(x) assignLabels(x, param$labels))
      MCMCtraces$orders <- lapply(MCMCresult$orders, order2var, varnames = param$labels)
    
    #MCMCtraces$DAGscores<-MCMCresult$DAGscores
    MCMCtraces$orderscores <- MCMCresult$orderscores
  }


  maxobj <- storemaxMCMC(MCMCresult, param)
  maxN <- which.max(MCMCresult$DAGscores)
  #maxobj$reach<-maxN

  if (scoreout) {
    if (chainout) { output <- 4 }
    else { output <- 3 }
  } else {
    if (chainout) { output <- 2 }
    else { output <- 1 }
  }

  result$DAG <- maxobj$DAG
  result$CPDAG <- graph2m(dag2cpdag(m2graph(result$DAG)))
  result$score <- maxobj$score
  result$maxorder <- maxobj$order
  result$info <- list()
  tabletime <- tableend - tablestart
  if (units(tabletime) == "mins") {
    tabletime <- as.numeric(tabletime * 60)
  }
  mcmctime <- mcmcend - tableend
  if (units(mcmctime) == "mins") {
    mcmctime <- as.numeric(mcmctime * 60)
  }
  result$info$times <- c(tabletime, mcmctime)
  result$trace <- MCMCresult$DAGscores
  switch(as.character(output),
         "1" = {
    # return only maximum DAG and score trace
  },
         "2" = {
    # return all MCMC all saved MCMC steps: incidence, DAGscore, orderscore and order and max result
    result$traceadd <- MCMCtraces
  },
         "3" = {
    # return max DAG, order, last search space incidence and all scoretables
    result$scoretable <- list()
    result$scoretable$adjacency <- startskel
    result$scoretable$tables <- scoretable
    result$scoretable$blacklist <- blacklist
    attr(result$scoretable, "class") <- "MCMCscoretab"
  },
         "4" = {
    # return all MCMC all saved MCMC steps,max result,last search space and scoretables
    result$scoretable <- list()
    result$scoretable$adjacency <- startskel
    result$scoretable$tables <- scoretable
    result$scoretable$blacklist <- blacklist
    attr(result$scoretable, "class") <- "MCMCscoretab"
    result$traceadd <- MCMCtraces
  }
  )
  attr(result, "class") <- "MCMCres"
  return(result)

}

#implements order MCMC on a defined search space, sampling version
#partly derived from <doi:10.1080/01621459.2015.1133426>
orderMCMCbaseFelix <- function(n, nsmall, startorder, iterations, stepsave, moveprobs, parenttable, scoretable, aliases, numparents,
                        rowmaps, scoresmatrices, numberofparentsvec, gamma = 1, bgnodes, matsize) {

  #n - number of nodes (background included)
  #nsmall - number of nodes excluding background
  #matsize - number of rows/columns in adjacency matrix 


  if (!is.null(bgnodes)) {
    mainnodes <- c(1:matsize)[-bgnodes]
  } else { mainnodes <- c(1:matsize) }

  currentpermy <- startorder #starting order represented as a permutation
  currentorderscores <- orderscoreBaseFelix(matsize, currentpermy[1:nsmall], c(1:nsmall), parenttable, aliases, numparents,
                                            rowmaps, scoretable, scoresmatrices, currentpermy) #starting score


    print(currentorderscores)
  currenttotallogscore <- sum(currentorderscores$totscores[mainnodes]) #log total score of all DAGs in the starting order
  currentDAG <- samplescoreplus1(matsize, mainnodes, currentorderscores, plus1lists = NULL, scoretable,
                               scoresmatrices, parenttable, numberofparentsvec, aliases) #score of a single DAG sampled from the starting order

  L1 <- list() # stores the adjacency matrix of a DAG sampled from the orders
  L2 <- vector() # stores its log BGe score
  L3 <- vector() # stores the log BGe score of the entire order
  L4 <- list() # stores the orders as permutations

  zlimit <- floor(iterations / stepsave) + 1 # number of outer iterations
  length(L1) <- zlimit
  length(L2) <- zlimit
  length(L3) <- zlimit
  length(L4) <- zlimit

  L1[[1]] <- currentDAG$incidence #starting DAG adjacency matrix
  L2[1] <- currentDAG$logscore #starting DAG score
  L3[1] <- currenttotallogscore #starting order score
  L4[[1]] <- currentpermy[1:nsmall] #starting order

  moveprobsstart <- moveprobs

  for (z in 2:zlimit) {
    #the MCMC chain loop with 'iteration' steps is in two parts
    for (count in 1:stepsave) {
      #since we only save the results to the lists each 'stepsave'


      chosenmove <- sample.int(4, 1, prob = moveprobs)
      if (chosenmove < 4) {
        # if it is 3 then we stay still

        proposedpermy <- currentpermy #sample a new order by swapping two elements
        switch(as.character(chosenmove),
               "1" = {
          # swap any two elements at random
          sampledelements <- sample.int(nsmall, 2, replace = FALSE) #chosen at random

        },
               "2" = {
          # swap any adjacent elements
          k <- sample.int(nsmall - 1, 1) #chose the smallest at random
          sampledelements <- c(k, k + 1)

        },
               "3" = {
          # swap any adjacent elements
          sampledpos <- sample.int(nsmall, 1)
        }, {
          # if neither is chosen, we have a problem
          print('The move sampling has failed!')
        }) #end switch

        if (chosenmove < 3) {

          proposedpermy[sampledelements] <- currentpermy[rev(sampledelements)] #proposed new order ???
          scorepositions <- c(min(sampledelements):max(sampledelements))
          rescorenodes <- proposedpermy[scorepositions] #we only need to rescore these nodes between the swapped elements to speed up the calculation
          proposedorderrescored <- orderscoreBaseFelix(matsize, rescorenodes, scorepositions, parenttable, aliases, numparents, rowmaps,
                                                scoretable, scoresmatrices, proposedpermy)
          proposedtotallogscore <- currenttotallogscore - sum(currentorderscores$totscores[rescorenodes]) + sum(proposedorderrescored$totscores[rescorenodes]) #and the new log total score by updating only the necessary nodes
          scoreratio <- exp((proposedtotallogscore - currenttotallogscore) * gamma) #acceptance probability

          if (runif(1) < scoreratio) {
            #Move accepted then set the current order and scores to the proposal

            currentpermy <- proposedpermy
            currentorderscores$therow[rescorenodes] <- proposedorderrescored$therow[rescorenodes]
            currentorderscores$totscores[rescorenodes] <- proposedorderrescored$totscores[rescorenodes]
            currenttotallogscore <- proposedtotallogscore

          }
        } else {

          neworder <- positionscorebase(matsize, nsmall, currentorderscores, sampledpos, currentpermy, aliases,
                                           rowmaps, numparents, scoretable, scoresmatrices)
          currentpermy <- neworder$order
          currentorderscores <- neworder$score
          currenttotallogscore <- sum(neworder$totscores)
        }

      }
    }
    currentDAG <- samplescoreplus1(matsize, mainnodes, currentorderscores, plus1lists = NULL, scoretable, scoresmatrices,
                                 parenttable, numberofparentsvec, aliases)
    L1[[z]] <- currentDAG$incidence #store adjacency matrix of a sampled DAG each 'stepsave'
    L2[z] <- currentDAG$logscore #and log score of a sampled DAG
    L3[z] <- currenttotallogscore #and the current order score
    L4[[z]] <- currentpermy[1:nsmall] #and store current order
  }
  result <- list()
  result$incidence <- L1
  result$DAGscores <- L2
  result$orderscores <- L3
  result$orders <- L4

  return(result)
}

#scores a single order base version (base neighbourhood)
orderscoreBaseFelix <- function(n, scorenodes, scorepositions, parenttable, aliases, numparents, rowmaps,
                         scoretable, scoresmatrices, permy) {

  orderscores <- vector("double", n)
  therows <- vector("integer", n)
  k <- 1
  for (i in scorenodes) {
    position <- scorepositions[k]
    if (position == n) {
      #no parents allowed, i.e. only first row, only first list
      orderscores[i] <- scoretable[[i]][1, 1]
      therows[i] <- c(2 ^ numparents[i])
    }
    else {
      bannednodes <- permy[1:position]
      allowednodes <- permy[(position + 1):n]
      bannedpool <- which(aliases[[i]] %in% bannednodes)
      if (numparents[i] == 0 || length(bannedpool) == 0) {
        #all parents allowed or no parents in the parent table
        therows[i] <- c(1)
      }
      else {
        therows[i] <- rowmaps[[i]]$backwards[sum(2 ^ bannedpool) / 2 + 1]
      }
      orderscores[i] <- scoresmatrices[[i]][therows[i], 1]
    }
    k <- k + 1
  }
  scores <- list()
  scores$therow <- therows
  scores$totscores <- orderscores
  return(scores)

}


#implements order MCMC on an extended search space sampling version
#partly derived from <doi:10.1080/01621459.2015.1133426>

orderMCMCplus1Felix <- function(n, nsmall, startorder, iterations, stepsave, moveprobs, parenttable, scoretable, aliases, numparents, rowmaps,
                           plus1lists, scoresmatrices, numberofparentsvec, gamma = 1, bgnodes, matsize) {
  
  #n - number of nodes (background included)
  #nsmall - number of nodes excluding background
  #matsize - number of rows/columns in adjacency matrix 
  print("startorder")
  print(startorder)
  if (!is.null(bgnodes)) {
    mainnodes <- c(1:matsize)[-bgnodes]
  } else { mainnodes <- c(1:matsize) }
  
  currentpermy <- startorder #starting order represented as a permutation
  currentorderscores <- orderscorePlus1Felix(matsize, currentpermy[1:nsmall], c(1:nsmall), parenttable, 
                                        aliases, numparents, rowmaps, plus1lists, scoretable, scoresmatrices, currentpermy) #starting score
  
  # This is the actual order score.
  currenttotallogscore <- sum(currentorderscores$totscores[mainnodes]) #log total score of all DAGs in the starting order
  currentDAG <- samplescoreplus1(matsize, mainnodes, currentorderscores, plus1lists, scoretable, scoresmatrices, parenttable, numberofparentsvec, aliases) #score of a single DAG sampled from the starting order
  
  L1 <- list() # stores the adjacency matrix of a DAG sampled from the orders
  L2 <- vector() # stores its log BGe score of a DAG
  L3 <- vector() # stores the log BGe score of the entire order
  L4 <- list() # stores the orders as permutations
  
  zlimit <- min(floor(iterations / stepsave) + 1, 1000) # number of outer iterations
  
  length(L1) <- zlimit
  length(L2) <- zlimit
  length(L3) <- zlimit
  length(L4) <- zlimit
  
  L1[[1]] <- currentDAG$incidence #starting DAG adjacency matrix
  L2[1] <- currentDAG$logscore #starting DAG score
  L3[1] <- currenttotallogscore #starting order score
  L4[[1]] <- currentpermy[1:nsmall] #starting order
  
  for (z in 2:zlimit) {
    #the MCMC chain loop with 'iteration' steps is in two parts
    for (count in 1:stepsave) {
      #since we only save the results to the lists each 'stepsave'
      
      
      chosenmove <- sample.int(4, 1, prob = moveprobs)
      if (chosenmove < 4) {
        # if it is 3 then we stay still
        
        proposedpermy <- currentpermy #sample a new order by swapping two elements
        switch(as.character(chosenmove),
               "1" = {
                 # swap any two elements at random
                 sampledelements <- sample.int(nsmall, 2, replace = FALSE) #chosen at random
                 
               },
               "2" = {
                 # swap any adjacent elements
                 k <- sample.int(nsmall - 1, 1) #chose the smallest at random
                 sampledelements <- c(k, k + 1)
                 
               },
               "3" = {
                 # swap any adjacent elements
                 sampledpos <- sample.int(nsmall, 1)
               }, {
                 # if neither is chosen, we have a problem
                 cat("The move sampling has failed!  \n")
               }) #end switch
        
        if (chosenmove < 3) {
          scorepositions <- c(min(sampledelements):max(sampledelements))
          proposedpermy[sampledelements] <- currentpermy[rev(sampledelements)] #proposed new order ???
          rescorenodes <- proposedpermy[scorepositions] #we only need to rescore these nodes between the swapped elements to speed up the calculation
          proposedorderrescored <- orderscorePlus1(matsize, rescorenodes, scorepositions, parenttable, aliases, numparents, rowmaps, plus1lists, scoretable, scoresmatrices, proposedpermy)
          proposedtotallogscore <- currenttotallogscore - sum(currentorderscores$totscores[rescorenodes]) + sum(proposedorderrescored$totscores[rescorenodes]) #and the new log total score by updating only the necessary nodes
          scoreratio <- exp((proposedtotallogscore - currenttotallogscore) * gamma) #acceptance probability
          
          if (runif(1) < scoreratio) {
            #Move accepted then set the current order and scores to the proposal
            currentpermy <- proposedpermy
            currenttotallogscore <- proposedtotallogscore
            currentorderscores$therow[rescorenodes] <- proposedorderrescored$therow[rescorenodes]
            currentorderscores$totscores[rescorenodes] <- proposedorderrescored$totscores[rescorenodes]
            currentorderscores$allowedlists[rescorenodes] <- proposedorderrescored$allowedlists[rescorenodes]
          }
        } else if (chosenmove == 3) {
          neworder <- positionscorePlus1(matsize, nsmall, currentorderscores, sampledpos, currentpermy, aliases,
                                         rowmaps, plus1lists, numparents, scoretable, scoresmatrices)
          currentpermy <- neworder$order
          currentorderscores <- neworder$score
          currenttotallogscore <- neworder$tot
        }
      }
    }
    currentDAG <- samplescoreplus1(matsize, mainnodes, currentorderscores, plus1lists, scoretable, scoresmatrices,
                                   parenttable, numberofparentsvec, aliases)
    L1[[z]] <- currentDAG$incidence #store adjacency matrix of a sampled DAG each 'stepsave'
    L2[z] <- currentDAG$logscore #and log score of a sampled DAG
    L3[z] <- currenttotallogscore #and the current order score
    L4[[z]] <- currentpermy[1:nsmall] #and store current order
  }
  result <- list()
  result$incidence <- L1
  result$DAGscores <- L2
  result$orderscores <- L3
  result$orders <- L4
  
  return(result)
}

#scores a single order base version (plus1 neighbourhood)

orderscorePlus1Felix <- function(n, scorenodes, scorepositions, parenttable, aliases, numparents,
                                 rowmaps, plus1lists, scoretable, scoresmatrices, permy) {
  
  print("n")
  print(n)
  print("numparents")
  print(numparents)
  # Scorenodes is the order? /Felix
  print("scorenodes")
  print(scorenodes)
  print("aliases")
  print(aliases)
  print("parenttable")
  print(parenttable)
  print("scoretable")
  print(scoretable)
  
  print("scorepositions")
  print(scorepositions)
  print(as.list(match.call()))
  orderscores <- vector("double", n)
  allowedscorelists <- vector("list", n)
  therows <- vector("integer", n) # what is this? /Felix
  k <- 1

  ## Seems lite there are 8 order scores. What does that mean? / Felix
  for (i in scorenodes) {
    position <- scorepositions[k]
    if (position == n) {
      #no parents allowed, i.e. only first row, only first list
      orderscores[i] <- scoretable[[i]][[1]][1, 1]
      allowedscorelists[[i]] <- c(1)
      therows[i] <- c(2 ^ numparents[i])
    }
    else {
      bannednodes <- permy[1:position]
      allowednodes <- permy[(position + 1):n]
      bannedpool <- which(aliases[[i]] %in% bannednodes)
      if (numparents[i] == 0 || length(bannedpool) == 0) {
        #all parents allowed
        therows[i] <- c(1)
      }
      else {
        therows[i] <- rowmaps[[i]]$backwards[sum(2 ^ bannedpool) / 2 + 1]
      }
      allowedscorelists[[i]] <- c(1, which(plus1lists$parents[[i]] %in% allowednodes) + 1)
      scoresvec <- scoresmatrices[[i]][therows[i], allowedscorelists[[i]]]
      maxallowed <- max(scoresvec)
      orderscores[i] <- maxallowed + log(sum(exp(scoresvec - maxallowed)))
    }
    k <- k + 1
  }
  scores <- list()
  scores$therow <- therows
  scores$allowedlists <- allowedscorelists
  scores$totscores <- orderscores
  
  print("scores")
  print(scores)
  print("therows")
  print(therows)
  return(scores)
  
  
}