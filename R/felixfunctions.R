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


# This is a rip off from orderMCMCmain.R , I think...
getScoreTable <- function(param, iterations, stepsave, MAP = TRUE, posterior = 0.5,
                        startorder = c(1:n), moveprobs, plus1 = FALSE, chainout = TRUE,
                        scoreout = FALSE, startspace = NULL, blacklist = NULL, gamma = 1, verbose = FALSE, alpha = NULL,
                        hardlimit = ifelse(plus1, 15, 22),
                        cpdag = FALSE, addspace = NULL) {
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

  if (is.null(startspace)) {
    startspace <- definestartspace(alpha, param, cpdag = cpdag, algo = "pc")
  }
  startskeleton <- 1 * (startspace & !blacklist)
  if (!is.null(addspace)) {
    startskel <- 1 * ((addspace | startskeleton) & !blacklist)
  } else { startskel <- startskeleton }

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


  # This is where score tables are obtained! /Felix
  # Should also take care of MAP=TRUE.
  scoretable <- scorepossibleparents.PLUS1(parenttable, plus1lists, n, param, updatenodes,
                                           rowmaps, numparents, numberofparentsvec)
  


  res <- list()
  res$tables <- scoretable
  res$adjacency <- startskel
  res$blacklist <- blacklist
  return(res)
}


orderMCMCFelix <- function(scorepar, MAP = TRUE, plus1 = TRUE, chainout = FALSE, scoreout = FALSE, moveprobs = NULL,
                    iterations = NULL, stepsave = NULL, alpha = 0.05, cpdag = FALSE, gamma = 1,
                    hardlimit = ifelse(plus1, 14, 20), verbose = FALSE,
                    startspace = NULL, blacklist = NULL, startorder = NULL, scoretable = NULL) {

  result <- orderMCMCmainFelix(param = scorepar, iterations, stepsave, startorder = startorder,
                                moveprobs = moveprobs, alpha = alpha, cpdag = cpdag, scoretable = scoretable,
                                plus1 = plus1, MAP = MAP, chainout = chainout, scoreout = scoreout,
                                startspace = startspace, blacklist = blacklist, gamma = gamma, verbose = verbose,
                                hardlimit = hardlimit)
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

  updatenodes <- c(1:n)
  maxorder <- startorder

  #defining startskel
 
  startskel <- scoretable$adjacency
  blacklist <- scoretable$blacklist
  scoretable <- scoretable$tables # remove the graph with no extra parent. 

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

  # move this to score tables? / Felix
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


  ########### Adding this ################
  posetparenttable <- poset(parenttable, numberofparentsvec, rowmaps, n, updatenodes)

  if (MAP == TRUE) {
    maxmatrices <- posetscoremax(posetparenttable, scoretable, numberofparentsvec,
                              rowmaps, n, plus1lists = plus1lists, updatenodes)
  } else {
    bannedscore <- poset.scores(posetparenttable, scoretable, ptab$numberofparentsvec, rowmaps,
                            n, plus1lists = plus1lists, ptab$numparents, updatenodes)
  }
  ##############################


  if (verbose) {
    cat(paste("score tables calculated, MCMC plus1 starts"), "\n")
    flush.console()
  }

  ret <- list()
  if (MAP == TRUE) {
   ret$bannedscore <- maxmatrices$maxmatrix
  } else {
    ret$bannedscore <- bannedscore
  }

  ret$plus1lists <- plus1lists
  ret$rowmaps <- rowmaps
  ret$ptab <- ptab
  return(ret)

}

