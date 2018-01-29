# A few basic utility functions

#' Get the probability gap between the highest and the second highest assignmnet
#' probabilities.
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
#' @param prob_assign A matrix of floats. Clonal assignment probability of N 
#' cells to K clones.
get_prob_label <- function(prob_assign){
  assign_clone <- rep(0, dim(prob_assign)[1])
  for (j in seq_len(length(assign_clone))){
    assign_clone[j] <- which.max(prob_assign[j,])
  }
  assign_clone
}
