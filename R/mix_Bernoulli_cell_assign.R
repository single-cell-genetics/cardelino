## Bernoulli mixtue model for cell assignment

#' EM algorithm for assigning cells to cloens and estimating Bernoulli mixture model.
#' The two Bernoulli components correspond to false positive and false negative.
#' 
#' @param A A matrix of integers. Number of alteration reads in variant i in cell j
#' @param D A matrix of integers. Number of total reads (depth) in variant i in cell j
#' @param C A matrix of binary values. The clone-variant configuration, whcih encodes
#' the phylogenetic tree structure. This is the output Z from Canopy
#' @param Psi A vector of float. The fractional size of clone, output P from Canopy
#' @param threshold A value of integer or float. The threshold on count or fraction of 
#' alteration reads
#' @param threshold_type A string. The type of threshold: count or fraction
#' @param max_iter A integer. The maximum number of iterations in EM algorithm. 
#' The real iteration may finish earlier.
#' 
#' @return a list containing \code{alpha}, a float denoting the estimated false positive 
#' rate, \code{beta}, a float denoting the estimated false negative rate, \code{prob}, the 
#' matrix of fitted probabilities of each cell belonging to each clone, and \code{logLik},
#' the log likelihood of the parameter based on the final cell assignment.
#' 
#' @import stats
#' 
#' @export
#' 
#' @examples
mix_Bernoulli_tree <- function(A, D, C, Psi, threshold=1, threshold_type="count", 
                               max_iter=100){
  if (dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2] || dim(A)[1] != dim(C)[1] || 
      dim(C)[2] != length(Psi)) {
    stop(paste("A and D must have the same size;\n ",
               "A and C must have the same number of variants;\n",
               "C and Psi must have the same number of clones",sep = ""))
  }
  
  ## preprocessing
  N <- dim(A)[1] # number of variants 
  M <- dim(A)[2] # number of cells
  K <- dim(C)[2] # number of clones
  
  if (threshold_type == "fraction"){
    Y <- (A/D >= threshold)
  } else{
    Y <- (A >= threshold)
  }
  
  C1 <- C
  C0 <- 1-C
  Y1 <- Y
  Y0 <- 1-Y
  Y1[is.na(Y1)] <- 0
  Y0[is.na(Y0)] <- 0
  
  S1 <- t(Y1) %*% C0  #False positive (Y=1 and C=0): alpha
  S2 <- t(Y0) %*% C0  #True negative  (Y=0 and C=0): 1-alpha
  S3 <- t(Y1) %*% C1  #True positive  (Y=1 and C=1): 1-beta
  S4 <- t(Y0) %*% C1  #False negative (Y=0 and C=1): beta
  
  ## random initialization for EM
  alpha <- stats::runif(1, 0.001, 0.5)
  beta  <- stats::runif(1, 0.001, 0.5)
  alpha_new <- stats::runif(1, 0.001, 0.5)
  beta_new  <- stats::runif(1, 0.001, 0.5)
  
  lik_mat <- (alpha**S1 * (1-alpha)**S2 * (1-beta)**S3 * beta**S4) * Psi
  prob_mat <- lik_mat / rowSums(lik_mat)
  
  ## EM iterations
  for(t in seq_len(max_iter)){
    #TODO: check convergence by the fold change on logLik
    if (alpha_new == alpha && beta_new == beta){
      print(paste("Total iterations: ", t))
      break
    } else{
      alpha <- alpha_new
      beta  <- beta_new
    }
    
    #E-step
    lik_mat <- (alpha**S1 * (1-alpha)**S2 * (1-beta)**S3 * beta**S4) * Psi
    prob_mat <- lik_mat / rowSums(lik_mat)
    
    #M-step
    alpha_new <- sum(prob_mat * S1) / sum(prob_mat * (S1+S2))
    beta_new  <- sum(prob_mat * S4) / sum(prob_mat * (S3+S4))
  }
  
  logLik <- sum(log(rowSums(lik_mat)))
  
  ## return values
  return_list <- list("alpha"=alpha, "beta"=beta, "prob"=prob_mat, "logLik"=logLik)
  return_list
}

