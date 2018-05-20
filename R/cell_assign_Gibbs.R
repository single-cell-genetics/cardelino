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
#' @param A_germ A matrix of integers. Number of alteration reads in germline 
#' heterozygous site for variant i cell j
#' @param D_germ A matrix of integers. Number of reads depth in germline 
#' heterozygous site for variant i cell j
#' @param prior0 numeric(2), alpha and beta parameters for the Beta prior 
#' distribution on the inferred false positive rate.
#' @param prior1 numeric(2), alpha and beta parameters for the Beta prior 
#' distribution on the inferred (1 - false negative) rate.
#' @param model A string. The model to use: Bernoulli or binomial
#' @param Bern_threshold An integer. the count threshold of alteration reads 
#' when using Bernoulli model
#' @param min_iter A integer. The minimum number of iterations in the Gibbs 
#' sampling. The real iteration may be longer utile the convergence.
#' @param max_iter A integer. The maximum number of iterations in the Gibbs 
#' sampling, even haven't passed the convergence diagnosis
#' @param wise A string, the wise of parameters for theta1: global, variant, 
#' element.
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
#' 
cell_assign_Gibbs <- function(A, D, C, Psi=NULL, A_germ=NULL, D_germ=NULL, 
                              prior0 = c(0.3, 29.7), prior1 = c(2.25, 2.65),
                              model="binomial", Bern_threshold=1, 
                              max_iter=10000, min_iter=1000, wise="variant", 
                              verbose=TRUE){
  if(is.null(Psi)){Psi <- rep(1/ncol(C), ncol(C))}
  if (dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2] || 
      dim(A)[1] != dim(C)[1] || dim(C)[2] != length(Psi)) {
    stop(paste0("A and D must have the same size;\n ",
                "A and C must have the same number of variants;\n",
                "C and Psi must have the same number of clones"))
  }
  if(sum(c("element", "variant", "global") == wise) == 0){
    stop(paste0("Input wise mode: ", wise, 
                ", while only supporting: element, variant, global"))
  }
  
  ## preprocessing
  N <- dim(A)[1]             # number of variants 
  M <- dim(A)[2]             # number of cells
  K <- dim(C)[2]             # number of clones
  if(is.null(A_germ)){A_germ <- matrix(NA, nrow=N, ncol=M)}
  if(is.null(D_germ)){D_germ <- matrix(NA, nrow=N, ncol=M)}
  
  A[which(D == 0)] <- NA
  D[which(D == 0)] <- NA
  A[(D>0) & is.na(A)] <- 0
  A_germ[which(D_germ == 0)] <- NA
  D_germ[which(D_germ == 0)] <- NA
  A_germ[(D_germ>0) & is.na(A_germ)] <- 0
  
  ### test: hack variants phasing
  idx_tmp <- (A > 0) || (A_germ == 0) || (A_germ == D_germ)
  A_germ[idx_tmp] <- D_germ[idx_tmp] <- NA
  idx_tmp <- which(A_germ > (D_germ-A_germ))
  A_germ[idx_tmp] <- D_germ[idx_tmp] - A_germ[idx_tmp]
  
  C1 <- C
  C0 <- 1-C
  if (model == "Bernoulli"){
    A1 <- A >= Bern_threshold        # bool values
    A2 <- A_germ >= Bern_threshold
    B1 <- 1 - A1
    B2 <- 1 - A2
    W_log <- matrix(0, N, M)
  }else{
    A1 <- A                  #number of alteration reads
    A2 <- A_germ             #number of alteration reads in germline var
    B1 <- D - A              #number of reference reads
    B2 <- D_germ - A_germ    #number of reference reads in germline var
    W_log <- lchoose(D, A)   #log binomial coefficients, need for likelihood
    # TODO: dose the likelihood include the germline var? (Yes)
  }
  A1[is.na(A1)] <- 0
  B1[is.na(B1)] <- 0
  A2[is.na(A2)] <- 0
  B2[is.na(B2)] <- 0
  
  #reads number list for each clone
  S1_list <- list()
  S2_list <- list()
  S3_list <- list()
  S4_list <- list()
  for(k in seq_len(K)){
    S1_list[[k]] <- A1 * C0[,k]
    S2_list[[k]] <- B1 * C0[,k]
    S3_list[[k]] <- A1 * C1[,k] + A2   # A2 here: part of likelihood
    S4_list[[k]] <- B1 * C1[,k] + B2
  }
  
  ## Prepare for sampling
  if(wise == "global"){
    idx_vec <- seq_len(1)         # For: theta1_all[t, ] <- theta1[idx_vec]
    idx_mat <- seq_len(N*M)       # For: theta1[idx_mat] <- theta1_all[t, ]
  }else if(wise == "variant"){
    idx_vec <- seq_len(N)
    idx_mat <- seq_len(N*M)
  }else if(wise == "element"){
    idx_vec <- which(A1+B1 > 0)
    idx_mat <- which(A1+B1 > 0)
  }
  n_element <- length(idx_vec)
  
  if(is.null(dim(prior1)) && length(prior1) == 2){
    prior1 <- t(matrix(rep(prior1, n_element), nrow=2)) #two variable to a matrix
  }
  if(!is.matrix(prior1)){
    stop("prior1 need to be a matrix of n_element x 2")
  }
  
  prob_all <- matrix(0, nrow=max_iter, ncol=M*K)
  logLik_mat <- matrix(0, nrow=M, ncol=K)
  logLik_all <- matrix(0, nrow=max_iter, ncol=1)
  assign_all <- matrix(0, nrow=max_iter, ncol=M)
  theta0_all <- matrix(0, nrow=max_iter, ncol=1)
  theta1_all <- matrix(0, nrow=max_iter, ncol=n_element)
  
  ## Random initialization
  theta0_all[1,1] <- stats::rbeta(1, prior0[1], prior0[2])
  theta1_all[1, ] <- stats::rbeta(rep(1,n_element), prior1[,1], prior1[,2])
  theta1 <- matrix(NA, nrow=N, ncol=M)
  
  ## Gibbs sampling
  for(it in 2:max_iter){
    # Update prob_mat
    theta0 <- theta0_all[it-1, 1]
    theta1[idx_mat] <- theta1_all[it-1,  ]
    for(k in seq_len(K)){
      logLik_mat[,k] <- (colSums(S1_list[[k]] * log(theta0),   na.rm=T) + 
                         colSums(S2_list[[k]] * log(1-theta0), na.rm=T) + 
                         colSums(S3_list[[k]] * log(theta1),   na.rm=T) + 
                         colSums(S4_list[[k]] * log(1-theta1), na.rm=T))
      logLik_mat[,k] <- logLik_mat[,k] + log(Psi[k])
    }
    logLik_mat_amplify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
    prob_mat <- exp(logLik_mat_amplify) / rowSums(exp(logLik_mat_amplify))
    prob_all[it, ] <- prob_mat
    
    logLik_vec <- rep(NA, nrow(logLik_mat))
    for (i in seq_len(nrow(logLik_mat))) {
      logLik_vec[i] <- matrixStats::logSumExp(logLik_mat[i,], na.rm = TRUE)
    }
    logLik_all[it] <- sum(logLik_vec, na.rm = TRUE)
    #Update logLikelihood (TODO: add W0 and W1)

    # Sample assignment
    for (j in seq_len(M)){
      assign_all[it,j] <- sample(seq_len(K), 1, replace=TRUE, prob=prob_mat[j,])
    }

    # Sample theta with assigned clones
    S1_wgt <- S2_wgt <- 0 # weighted S1
    S3_wgt <- S4_wgt <- matrix(0, nrow=N, ncol=M)
    for(k in seq_len(K)){
      idx <- which(assign_all[it,] == k)
      S1_wgt <- S1_wgt + sum(S1_list[[k]][,idx], na.rm=T)
      S2_wgt <- S2_wgt + sum(S2_list[[k]][,idx], na.rm=T)
      S3_wgt[,idx] <- S3_wgt[,idx] + S3_list[[k]][,idx]
      S4_wgt[,idx] <- S4_wgt[,idx] + S4_list[[k]][,idx]
    }
    
    if(wise == "global"){
      S3_wgt[,] <- sum(S3_wgt, na.rm=T)
      S4_wgt[,] <- sum(S4_wgt, na.rm=T)
    }else if(wise == "variant"){
      S3_wgt[,] <- rowSums(S3_wgt, na.rm=T)
      S4_wgt[,] <- rowSums(S4_wgt, na.rm=T)
    }
    theta0_all[it, 1] <- stats::rbeta(1, prior0[1] + S1_wgt, prior0[2] + S2_wgt)
    theta1_all[it,  ] <- stats::rbeta(rep(1, n_element),
                                      prior1[,1] + S3_wgt[idx_vec], 
                                      prior1[,2] + S4_wgt[idx_vec])
    
    #Check convergence
    if ((it >= min_iter) && (it %% 100 == 0)){
      FLAG <- Geweke_Z(theta0_all[1:it, 1]) <= 2
      for(n in seq_len(ncol(theta1_all))){
        if(FLAG == FALSE){break}
        FLAG <- Geweke_Z(theta1_all[1:it, n]) <= 2
      }
      if(FLAG){break}
    }
  }
  if(verbose){print(paste("Converged in", it, "iterations."))}
  
  ## Return values
  n_buin = ceiling(it * 0.25)
  
  a <- A1[idx_mat]
  d <- A1[idx_mat] + B1[idx_mat]
  binom_pdf1 <- binom_pdf0 <- rep(0, n_element)
  for(i in seq(n_buin, it)){
    binom_pdf1 <- binom_pdf1 + stats::dbinom(a, size=d, prob=theta1_all[i,])
    binom_pdf0 <- binom_pdf0 + stats::dbinom(a, size=d, prob=theta0_all[i])
  }
  prob_variant <- matrix(NA, nrow=N, ncol=M)
  prob_variant[idx_mat] <- binom_pdf1 / (binom_pdf1 + binom_pdf0)
  row.names(prob_variant) <- row.names(A)
  colnames(prob_variant) <- colnames(A)
  
  # prob_mat <- get_Gibbs_prob(assign_all[1:it,], buin_in=0.25)
  prob_mat <- matrix(colMeans(prob_all[n_buin:it, ]), nrow = M)
  row.names(prob_mat) <- colnames(A)
  colnames(prob_mat) <- colnames(C)
  
  theta0 <-mean(theta0_all[n_buin:it,])
  theta1[idx_mat] <- colMeans(as.matrix(theta1_all[n_buin:it,]))
  
  return_list <- list("theta0"=theta0, "theta1"=theta1,
                      "theta0_all"=as.matrix(theta0_all[1:it,]), 
                      "theta1_all"=as.matrix(theta1_all[1:it,]),
                      "element"=idx_mat, "logLik"=logLik_all[1:it], 
                      "prob_all"=prob_all[1:it,], 
                      "prob"=prob_mat, "prob_variant"=prob_variant)
  return_list
}


#' Geweke diagnostic for MCMC sampling.
#' @param X A vector of MCMC samples. Note, it is one-dimentional.
#' @param first A float between 0 and 1. The initial region of MCMC chain.
#' @param last A float between 0 and 1. The final region of MCMC chain.
#' 
#' @return an absolute value of the Z score.
Geweke_Z <- function(X, first=0.1, last=0.5){
  N = length(X)
  A = X[1:floor(first*N)]
  B = X[ceiling(last*N):N]
  Z = abs(mean(A) - mean(B)) / sqrt(var(A) + var(B))
  
  return(Z)
}


#' Get the probability matrix from Gibbs sampling on assignment
#' @param assign_Gibbs A matrix of T iterations on M cells
#' @param buin_in A float of the fraction of buin-in iterations
#' 
#' @return A probability matrix of M cells and K clones
get_Gibbs_prob <- function(assign_Gibbs, buin_in=0.25){
  T <- nrow(assign_Gibbs)
  T_buin_in <- ceiling(T * buin_in)
  assign_Gibbs <- assign_Gibbs[T_buin_in:T, ]
  K <- max(assign_Gibbs)
  prob_mat <- matrix(0, nrow=ncol(assign_Gibbs), K)
  for(m in seq_len(nrow(prob_mat))){
    for(k in seq_len(K)){
      prob_mat[m, k] <- mean(assign_Gibbs[,m] == k)
    }
  }
  prob_mat
}


get_logLik <- function(A, D, Config, theta1, theta0){
  A[which(D == 0)] <- NA
  D[which(D == 0)] <- NA
  
  P0_mat <- dbinom(A, D, theta0, log=T)
  P1_mat <- dbinom(A, D, theta1, log=T)
  
  P0_mat[which(is.na(P0_mat))] <- 0
  P1_mat[which(is.na(P1_mat))] <- 0
  
  logLik_mat <- t(P0_mat) %*% (1-Config) + t(P1_mat) %*% Config
  logLik_mat_amplify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
  prob_mat <- exp(logLik_mat_amplify) / rowSums(exp(logLik_mat_amplify))
  
  logLik_vec <- rep(NA, nrow(logLik_mat))
  for (i in seq_len(nrow(logLik_mat))) {
    logLik_vec[i] <- matrixStats::logSumExp(logLik_mat[i,], na.rm = TRUE)
  }
  logLik <- sum(logLik_vec, na.rm = TRUE)
  
  list(prob_mat, logLik)
}
