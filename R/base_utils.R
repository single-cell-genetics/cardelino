# A few basic utility functions

#' Get the probability gap between the highest and the second highest assignmnet
#' probabilities.
#' 
#' @export
#' 
#' @param prob_assign A matrix of floats. Clonal assignment probability of N 
#' cells to K clones.
get_prob_gap <- function(prob_assign){
  prob_gap <- rep(0, dim(prob_assign)[1])
  for (j in seq_len(length(prob_gap))){
    prob_sorted <- sort(prob_assign[j,], decreasing = TRUE)
    prob_gap[j] <- prob_sorted[1] - prob_sorted[2]
  }
  prob_gap
}

#' Get the clone label from the assignmnet probabilities. 
#' 
#' Note, when multiple clones have the same assignment probability, only the 
#' earliest clone will be return. Usually, these cells should be filtered for 
#' analysis.
#' 
#' @export
#' 
#' @param prob_assign A matrix of floats. Clonal assignment probability of N 
#' cells to K clones.
get_prob_label <- function(prob_assign){
  assign_clone <- rep(0, dim(prob_assign)[1])
  for (j in seq_len(length(assign_clone))){
    assign_clone[j] <- which.max(prob_assign[j,])
  }
  assign_clone
}


#' Calculate the SNV specific allelic frequency (True postive rate) and error
#' rate (False possitive rate)
#' 
#' @export
get_snv_theta <- function(C, D, A, prob_mat, threshold=1){
  ## preprocessing
  N <- dim(A)[1] # number of variants 
  M <- dim(A)[2] # number of cells
  K <- dim(C)[2] # number of clones
  
  A[is.na(A)] <- 0  # number of alteration reads
  D[is.na(D)] <- 0  # number of total reads
  
  snv_theta <- matrix(0, nrow = N, ncol = 2)
  for (n in seq_len(N)){
    a0 <- sum(A[n,] * prob_mat * (1-C[n,]))  #Alteration reads without variants
    d0 <- sum(D[n,] * prob_mat * (1-C[n,]))  #Total reads without variants (C=0)
    a1 <- sum(A[n,] * prob_mat * C[n,])  #Alteration reads with variants
    d1 <- sum(D[n,] * prob_mat * C[n,])  #Total reads with variants (C=1)
    
    snv_theta[n,1] <- a0 / d0
    snv_theta[n,2] <- a1 / d1
  }
  snv_theta
}