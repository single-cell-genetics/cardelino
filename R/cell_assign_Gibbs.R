## Single cell assignment to phylogenetic clones with Gibbs sampling for 
## Bernoulli mixtue model and binomial mixture model

#' Gibbs sampling the cell assignments to cloens and the parameters of the
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
#' @param prior0 numeric(2), alpha and beta parameters for the Beta prior 
#' distribution on the inferred false positive rate.
#' @param prior1 numeric(2), alpha and beta parameters for the Beta prior 
#' distribution on the inferred (1 - false negative) rate.
#' @param model A string. The model to use: Bernoulli or binomial
#' @param threshold A value of integer or float. the threshold on count or 
#' fraction of alteration reads when using Bernoulli model
#' @param threshold_type A string. The type of threshold: count or fraction
#' @param min_iter A integer. The minimum number of iterations in the Gibbs 
#' sampling. The real iteration may be longer utile the convergence.
#' @param max_iter A integer. The maximum number of iterations in the Gibbs 
#' sampling, even haven't passed the convergence diagnosis
#' 
#' @return a list containing \code{theta}, a vector of two floats denoting the 
#' parameters of the two componets of the base model, i.e., mean of Bernoulli or 
#' binomial model given variant exists or not, \code{assign}, the matrix of 
#' clonal labels of each cell is assigned, and \code{logLik}, the log likelihood 
#' of the parameters during sampling. All returns have a length of H, as the 
#' Gibbs sampling chain.
#' 
#' @import stats
#' 
#' @export
#' 
#' @examples
cell_assign_Gibbs <- function(A, D, C, Psi, prior0 = c(1, 1), prior1 = c(1, 1),
                              model = "Bernoulli", threshold = 1, 
                              threshold_type = "count", min_iter = 1000, 
                              max_iter = 50000) {
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
  assign_all <- matrix(1, max_iter, M)
  logLik_all <- numeric(max_iter)
  
  theta_all <- matrix(1, max_iter, 2)
  theta_all[1,1] <- stats::runif(1, 0.001, 0.25)
  theta_all[1,2] <- stats::runif(1, 0.3, 0.6)
  
  lik_mat <- (theta_all[1,1]**S1 * (1-theta_all[1,1])**S2 * 
              theta_all[2,1]**S3 * (1-theta_all[2,2])**S4) * Psi
  prob_mat <- lik_mat / rowSums(lik_mat)
  
  logLik_all[1] <- sum(log(rowSums(lik_mat*exp(W0+W1))))
  for (j in seq_len(M)){
    assign_all[1,j] <- sample(seq_len(K), 1, replace=TRUE, prob=prob_mat[j,])
  }
  
  ## Gibbs sampling
  for(t in 2:max_iter){
    ss1 <- 0
    ss2 <- 0
    ss3 <- 0
    ss4 <- 0
    #sample theta
    for (j in seq_len(M)){
      ss1 <- ss1+S1[j,assign_all[t-1,j]]
      ss2 <- ss2+S2[j,assign_all[t-1,j]]
      ss3 <- ss3+S3[j,assign_all[t-1,j]]
      ss4 <- ss4+S4[j,assign_all[t-1,j]]
    }
    theta_all[t, 1] <- rbeta(1, prior0[1] + ss1, prior0[2]+ss2, ncp = 0)
    theta_all[t, 2] <- rbeta(1, prior1[1] + ss3, prior1[2] + ss4, ncp = 0)
    
    #sample assignment
    lik_mat <- (theta_all[t,1]**S1 * (1-theta_all[t,1])**S2 * 
                theta_all[t,2]**S3 *(1-theta_all[t,2])**S4) * Psi
    prob_mat <- lik_mat / rowSums(lik_mat)
    
    logLik_all[t] <- sum(log(rowSums(lik_mat*exp(W0+W1))))
    for (j in seq_len(M)){
      assign_all[t,j] <- sample(seq_len(K), 1, replace=TRUE, prob=prob_mat[j,])
    }
    
    #check convergence
    if (t >= min_iter && t%%100 == 0){
      if (Geweke_Z(theta_all[1:t,1]) <= 2 && Geweke_Z(theta_all[1:t,2]) <= 2){
        break
      }
    }
  }
  
  ## return values
  return_list <- list("theta"=theta_all[1:t,], "assign"=assign_all[1:t,], 
                      "logLik"=logLik_all[1:t])
  return_list
}


#' Geweke diagnostic for MCMC sampling.
#' @param X A vector of MCMC samples. Note, it is one-dimentional.
#' @param first A float between 0 and 1. The initial region of MCMC chain.
#' @param last A float between 0 and 1. The final region of MCMC chain.
#' 
#' @return an absolute value of the Z score.
Geweke_Z <-function(X, first=0.1, last=0.5){
  N = length(X)
  A = X[1:floor(first*N)]
  B = X[ceiling(last*N):N]
  Z = abs(mean(A) - mean(B)) / sqrt(var(A) + var(B))
  
  return(Z)
}

