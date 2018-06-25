## Donor id inference with Gibbs sampling.

#' Gibbs sampler of the cell donor id and genotypes
#' 
#' @param A A matrix of integers. Number of alteration reads in variant i cell j
#' @param D A matrix of integers. Number of reads depth in variant i cell j
#' @param K A integer. Number for expected donors
#' @param prior0 numeric(2), alpha and beta parameters for the Beta prior 
#' distribution of the ALT read when variant does not exist.
#' @param prior1 numeric(2), alpha and beta parameters for the Beta prior 
#' distribution of the ALT read when variant does exists.
#' @param min_iter A integer. The minimum number of iterations in the Gibbs 
#' sampling. The real iteration may be longer utile the convergence.
#' @param max_iter A integer. The maximum number of iterations in the Gibbs 
#' sampling, even haven't passed the convergence diagnosis
#' @param buin_in A float between 0 and 1. The fraction for buil-in in MCMC 
#' samplers.
#' 
#' @return a list containing 
#' \code{prob}, a matrix (M x K) of probability of cell j assign to donor k. 
#' \code{theta0}, a float, the averaged theta0 in samples.
#' \code{theta1}, a float, the averaged theta1 in samples.
#' \code{theta0_all}, a vector of all sampled theta0.
#' \code{theta1_all}, a vector of all sampled theta1.
#' \code{logLik}, a vector of all sampled log likelihood.
#' \code{assign}, a matrix of donor assignment of each cell in all samples.
#' \code{Config}, a matrix of genotype configuration, averaged all samples.
#' \code{Config_all}, a matrix of all sample genotype configuration. Each 
#' column is a melted N x K matrix.
#' 
#' @import matrixStats
#' 
#' @export
#' 
donor_id <- function(A, D, K, prior0 = c(0.02, 40), prior1 = c(2.4, 2.4),
                     min_iter=500, max_iter=50000, buin_in=0.25){
  ##TODO: 1) detect doublet with BF, 2) support inputed Config
  
  ## check input
  if (dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2]) {
    stop(paste0("A and D must have the same size;\n "))
  }
  
  ## preprocessing
  N <- dim(A)[1]             # number of variants 
  M <- dim(A)[2]             # number of cells
  A[which(D == 0)] <- NA
  D[which(D == 0)] <- NA
  A[(D > 0) & is.na(A)] <- 0
  
  A1 <- A                    # number of alteration reads
  B1 <- D - A                # number of reference reads
  W_log <- lchoose(D, A)     # log binomial coefficients, for likelihood
  # TODO: dose the likelihood include the germline var? (Yes)
  ### Double check
  # A1[is.na(A1)] <- 0
  # B1[is.na(B1)] <- 0
  
  logLik_all <- matrix(0, nrow = max_iter, ncol = 1)
  theta0_all <- matrix(0, nrow = max_iter, ncol = 1)
  theta1_all <- matrix(0, nrow = max_iter, ncol = 1)
  assign_all <- matrix(0, nrow = max_iter, ncol = M)
  Config_all <- matrix(0, nrow = max_iter, ncol = N*K)
  
  ## Random initialization
  #Potential extension: allele specific expression, K-Means warm initialization
  theta0_all[1,1] <- stats::rbeta(1, prior0[1], prior0[2])
  theta1_all[1,1] <- stats::rbeta(1, prior1[1], prior1[2])
  Iden_mat <- t(stats::rmultinom(M, size = 1, prob = rep(1/K, K)))
  
  ## Gibbs sampling
  for (it in 2:max_iter) {
    # Update element probability. 
    ## P0: prob without variant; P1: prob with variant
    P0_mat <- dbinom(A1, A1 + B1, theta0_all[it - 1, 1], log = TRUE)
    P1_mat <- dbinom(A1, A1 + B1, theta1_all[it - 1, 1], log = TRUE)
    P0_mat[which(is.na(P0_mat))] <- 0
    P1_mat[which(is.na(P1_mat))] <- 0
    
    # Sample configuration
    oddR_log <- P1_mat %*% Iden_mat - P0_mat %*% Iden_mat
    oddR_log[which(oddR_log > 50)] <- 50
    oddR_log[which(oddR_log < -50)] <- -50
    prob_mat <- exp(oddR_log) / (exp(oddR_log) + 1)    
    # oddR_log_amptify <- oddR_log - matrixStats::rowMaxs(oddR_log)
    # prob_mat <- exp(oddR_log_amptify) / (exp(oddR_log_amptify) + 1)
    
    Conf_mat <- matrix(stats::rbinom(N*K, size = 1, prob_mat), nrow = N)
    Config_all[it, ] <- Conf_mat
    
    # Sample assignment
    logLik_mat <- t(P0_mat) %*% (1 - Conf_mat) + t(P1_mat) %*% Conf_mat
    logLik_mat_amptify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
    prob_mat <- exp(logLik_mat_amptify) / rowSums(exp(logLik_mat_amptify))
    for (j in seq_len(M)) {
      Iden_mat[j, ] <- stats::rmultinom(1, size = 1, prob = prob_mat[j,])
    }
    assign_all[it, ] <- Iden_mat %*% seq_len(K)
    
    # Update logLikelihood (TODO: add W0 and W1)
    logLik_all[it] <- sum(log(rowSums(exp(logLik_mat), na.rm = TRUE)), 
                          na.rm = TRUE)
    
    # Sample parameters
    Geno_mat <- Conf_mat %*% t(Iden_mat) # genotype matrix: N * M
    S1 <- sum(A1 * (1 - Geno_mat), na.rm = TRUE)
    S2 <- sum(B1 * (1 - Geno_mat), na.rm = TRUE)
    S3 <- sum(A1 * Geno_mat, na.rm = TRUE)
    S4 <- sum(B1 * Geno_mat, na.rm = TRUE)
    theta0_all[it, 1] <- stats::rbeta(1, prior0[1] + S1, prior0[2] + S2)
    theta1_all[it, 1] <- stats::rbeta(1, prior1[1] + S3, prior1[2] + S4)
    
    # Check convergence
    if (it >= min_iter && it %% 100 == 0) {
      FLAG <- Geweke_Z(theta0_all[1:it,1]) <= 2
      for (n in seq_len(ncol(theta1_all))) {
        if (FLAG == FALSE) {break}
        FLAG <- Geweke_Z(theta1_all[1:it,n]) <= 2
      }
      if (FLAG) {break}
    }
  }
  print(paste("Converged in", it, "iterations."))
  
  ## return values
  n_buin = ceiling(it * buin_in)
  prob_mat <- get_Gibbs_prob(assign_all[1:it, ])
  row.names(prob_mat) <- colnames(A)
  colnames(prob_mat) <- paste0("clone", seq_len(ncol(prob_mat)))
  
  theta0 <- mean(theta0_all[n_buin:it,])
  theta1 <- mean(theta1_all[n_buin:it,])
  
  Config <- matrix(colMeans(Config_all[n_buin:it, ]), nrow = N)
  
  return_list <- list("prob" = prob_mat, "theta0" = theta0, "theta1" = theta1, 
                      "theta0_all" = as.matrix(theta0_all[1:it,]), 
                      "theta1_all" = as.matrix(theta1_all[1:it,]),
                      "logLik" = logLik_all[1:it], "assign" = assign_all[1:it,],
                      "Config" = Config, "Config_all" = Config_all[1:it,])
  return_list
}
