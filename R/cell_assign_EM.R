## Single cell assignment to phylogenetic clones with EM algorithm for 
## Bernoulli mixtue model and binomial mixture model

#' EM algorithm for assigning cells to clones and estimating parameters in 
#' mixture model.
#' The two Bernoulli components correspond to false positive and false negative.
#' The two binomial components correspond to reads distribution without and with
#' variant.
#' 
#' @param A A matrix of integers. Number of alteration reads in variant i cell j
#' @param D A matrix of integers. Number of reads depth in variant i cell j
#' @param C A matrix of binary values. The clone-variant configuration, whcih 
#' encodes the phylogenetic tree structure. This is the output Z of Canopy
#' @param Psi A vector of float. The fractions of each clone, output P of Canopy
#' @param model A string. The model to use: Bernoulli or binomial
#' @param threshold A value of integer or float. the threshold on count or 
#' fraction of alteration reads when using Bernoulli model
#' @param threshold_type A string. The type of threshold: count or fraction
#' @param max_iter A integer. The maximum number of iterations in EM algorithm. 
#' The real iteration may finish earlier.
#' @param min_iter A integer. The minimum number of iterations in EM algorithm.
#' @param logLik_threshold A float. The threshold of logLikelihood increase for
#' detecting convergence.
#' @param verbose logical(1), should the function output verbose information as
#' it runs? 
#' 
#' @return a list containing \code{alpha}, a float denoting the estimated false 
#' positive rate, \code{beta}, a float denoting the estimated false negative 
#' rate, \code{prob}, the matrix of fitted probabilities of each cell belonging 
#' to each clone, and \code{logLik},
#' the log likelihood of the parameter based on the final cell assignment.
#' 
#' @return a list containing \code{theta}, a vector of two floats denoting the 
#' parameters of the two componets of the base model, i.e., mean of Bernoulli or 
#' binomial model given variant exists or not, \code{prob}, the matrix of 
#' posterior probabilities of each cell belonging to each clone with fitted 
#' parameters, and \code{logLik}, the log likelihood of the final parameters.
#' 
#' @import stats
#' 
#' @export
#' 
#' @examples
cell_assign_EM <- function(A, D, C, Psi=NULL, model="Bernoulli", threshold=1, 
                           threshold_type="count", max_iter=1000, min_iter=10, 
                           logLik_threshold=1e-5, verbose=TRUE) {
  if(is.null(Psi)){Psi <- rep(1/ncol(C), ncol(C))}
  if (dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2] || 
      dim(A)[1] != dim(C)[1] || dim(C)[2] != length(Psi)) {
    stop(paste("A and D must have the same size;\n ",
               "A and C must have the same number of variants;\n",
               "C and Psi must have the same number of clones",sep = ""))
  }
  
  ## preprocessing
  N <- dim(A)[1] # number of variants 
  M <- dim(A)[2] # number of cells
  K <- dim(C)[2] # number of clones
  A[which(D==0)] <- NA
  D[which(D==0)] <- NA
  
  C1 <- C
  C0 <- 1-C
  if (model == "Bernoulli"){
    if (threshold_type == "fraction"){
      A1 <- (A/D >= threshold)
    } else{
      A1 <- (A >= threshold)
    }
    B1 <- 1 - A1
    W_log <- matrix(0, N, M)
  } else{
    A1 <- A                #number of alteration reads
    B1 <- D - A            #number of reference reads
    W_log <- lchoose(D,A)  #log binomial coefficients, need for likelihood
  }
  A1[is.na(A1)] <- 0
  B1[is.na(B1)] <- 0
  W_log[is.na(W_log)] <- 0
  
  S1 <- t(A1) %*% C0  #Alteration reads without variants (C=0): theta[1]
  S2 <- t(B1) %*% C0  #Reference reads without variants  (C=0): 1-theta[1]
  S3 <- t(A1) %*% C1  #Alteration reads with variants (C=1): theta[2]
  S4 <- t(B1) %*% C1  #Reference reads with variants  (C=1): 1-theta[2]
  W0 <- t(W_log) %*% C0 #log product of Binomial cofficient for C=0
  W1 <- t(W_log) %*% C1 #log product of Binomial cofficient for C=1
  
  ## random initialization for EM
  theta <- c(stats::runif(1, 0.001, 0.25), stats::runif(1, 0.3, 0.6))
  logLik <- 0
  logLik_new <- logLik + logLik_threshold * 2
  
  ## EM iterations
  for(it in seq_len(max_iter)){
    # Check convergence
    if((it > min_iter) && ((logLik_new - logLik) < logLik_threshold)) {break}
    logLik <- logLik_new
    
    #E-step
    logLik_mat <- (S1 * log(theta[1]) + S2 * log(1 - theta[1]) + 
                   S3 * log(theta[2]) + S4 * log(1 - theta[2]))
    logLik_mat <- t(t(logLik_mat) + log(Psi))# + W0 + W1
    #logLik_new <- sum(log(rowSums(exp(logLik_mat), na.rm = TRUE)), na.rm = TRUE)
    logLik_vec <- rep(NA, nrow(logLik_mat))
    for (i in seq_len(nrow(logLik_mat))) {
      logLik_vec[i] <- matrixStats::logSumExp(logLik_mat[i,], na.rm = TRUE)
    }
    logLik_new <- sum(logLik_vec, na.rm = TRUE)
    logLik_mat_amplify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
    prob_mat <- exp(logLik_mat_amplify) / rowSums(exp(logLik_mat_amplify))
    
    #M-step
    theta[1] <- sum(prob_mat * S1) / sum(prob_mat * (S1+S2))
    theta[2] <- sum(prob_mat * S3) / sum(prob_mat * (S3+S4))
  }
  if(verbose){print(paste("Total iterations:", it))}
  
  ## return values
  return_list <- list("theta"=theta, "prob"=prob_mat, "logLik"=logLik)
  return_list
}

