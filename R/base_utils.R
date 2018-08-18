# A few basic utility functions

#' Get the collapsed probability value
#' 
#' the collapsed probability value includes 1) the best probability, 2) the
#' second best probability, and 3) the delta between the best and the second.
#' 
#' @export
#' 
#' @param prob_assign A matrix of floats. Clonal assignment probability of N 
#' cells to K clones.
#' @param mode A string, the mothod for defining scores for filtering cells: 
#' best, second and delta. Best: highest probability of a cell to K clones, 
#' similarly for second. delta is the difference between the best and second.
#' 
get_prob_value <- function(prob_assign, mode="best") {
    prob_val <- rep(0, nrow(prob_assign))
    for (i in seq_len(length(prob_val))) {
        prob_sorted <- sort(prob_assign[i,], decreasing = TRUE)
        if (mode == "delta") {
            prob_val[i] <- prob_sorted[1] - prob_sorted[2]
        } else if (mode == "second") {
            prob_val[i] <- prob_sorted[2]
        } else { #default mode: best
            prob_val[i] <- prob_sorted[1]
        }
    }
    prob_val
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
#' @param C A matrix of binary values. The clone-variant configuration, which 
#' encodes the phylogenetic tree structure as presence (1) or absence (0) of 
#' each variant (rows) in each clone (columns). This is the output Z of Canopy.
#' @param D A matrix of integers. Number of reads depth in variant i cell j.
#' @param A A matrix of integers. Number of reads supporting the alternative 
#' allele in variant i cell j.
#' @param prob_mat matrix of cell-clone assignment probabilities.
#' 
#' @export
get_snv_theta <- function(C, D, A, prob_mat){
    ## preprocessing
    N <- dim(A)[1] # number of variants 
    M <- dim(A)[2] # number of cells
    K <- dim(C)[2] # number of clones
    
    A[is.na(A)] <- 0  # number of alteration reads
    D[is.na(D)] <- 0  # number of total reads
    
    snv_theta <- matrix(0, nrow = N, ncol = 2)
    for (n in seq_len(N)) {
        a0 <- sum(A[n,] * prob_mat * (1 - C[n,]))  #Alteration reads without variants
        d0 <- sum(D[n,] * prob_mat * (1 - C[n,]))  #Total reads without variants (C=0)
        a1 <- sum(A[n,] * prob_mat * C[n,])  #Alteration reads with variants
        d1 <- sum(D[n,] * prob_mat * C[n,])  #Total reads with variants (C=1)
        
        snv_theta[n,1] <- a0 / d0
        snv_theta[n,2] <- a1 / d1
    }
    snv_theta
}
