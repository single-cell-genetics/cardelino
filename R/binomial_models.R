## This file is for binomial distributions

#' EM algorithm for estimating binomial mixture model
#' 
#' 
#' @param k A vector of integers. number of success
#' @param n A vector of integers. number of trials
#' @param n_components A number. number of components
#' @param tol numeric(1), tolerance value for convergence between iterations
#' 
#' @return a list containing \code{p}, a vector of floats between 0 and 1 giving
#' the estimated success probability for each component, \code{psi}, estimated
#' fraction of each component in the mixture, and \code{prob}, the matrix of 
#' fitted probabilities of each observation belonging to each component.
#' 
#' @import stats
#' 
#' @export
#' 
#' @examples
#' n1 <- array(sample(1:30, 50, replace = TRUE))
#' n2 <- array(sample(1:30, 200, replace = TRUE))
#' k1 <- apply(n1, 1, rbinom, n = 1, p = 0.5)
#' k2 <- apply(n2, 1, rbinom, n = 1, p = 0.01)
#' RV <- mixBinom(c(k1, k2), c(n1, n2))
mixBinom <- function(k, n, n_components = 2, tol = 1e-08) {
    if (length(k) != length(n)) {
        stop("n and k must be of the same length (number of observations)\n")
    }
    S <- length(n)
    ## Random initialzation on parameters
    p <- stats::runif(n_components, 0.0, 1.0)
    psi <- stats::runif(n_components, 0.0, 1.0)
    psi <- psi / sum(psi)
    
    p_new <- stats::runif(n_components, 0.0, 1.0)
    psi_new <- stats::runif(n_components, 0.0, 1.0)
    psi_new <- psi_new / sum(psi_new)
    
    prob_mat <- matrix(stats::runif(n_components * S, 0, 1), nrow = S, 
                       byrow = TRUE)
    prob_mat <- prob_mat / rowSums(prob_mat)
    
    ## Iterations
    while (!(all(abs(p - p_new) < tol))) {
        ## E-step:
        for (j in seq_len(n_components)) {
            prob_mat[,j] <- psi[j] * stats::dbinom(k, size = n, prob = p_new[j],
                                            log = FALSE)
        }
        prob_mat <- prob_mat / rowSums(prob_mat)
        
        ## M-step
        for (j in seq_len(n_components)) {
            p[j] <- p_new[j]
            p_new[j] <- sum(prob_mat[,j] * k) / sum(prob_mat[,j] * n)
            psi[j] <- psi_new[j]
            psi_new[j] <- sum(prob_mat[,j]) / S
        }
    }
    
    ## return values
    return_list <- list("p" = p_new, "psi" = psi_new, "prob" = prob_mat)
    return_list
}


#' Predicted probability from learned binomial mixture model
#' 
#' @param k A vector of integers. number of success
#' @param n A vector of integers. number of trials
#' @param p a vector of binomial success probabilities
#' @param psi A float between 0 and 1. fraction of each component
#' 
#' @export
#' 
#' @examples
#' n1 <- array(sample(1:30, 50, replace = TRUE))
#' n2 <- array(sample(1:30, 200, replace = TRUE))
#' k1 <- apply(n1, 1, rbinom, n = 1, p = 0.5)
#' k2 <- apply(n2, 1, rbinom, n = 1, p = 0.01)
#' RV <- mixBinom(c(k1, k2), c(n1, n2))
#' prob <- predMixBinom(3, 10, RV$p, RV$psi)
predMixBinom <- function(k, n, p, psi) {
    if (length(p) != length(psi)) {
        stop("p and psi must be of the same length (number of mixture components)\n")
    }
    if (length(k) != length(n)) {
        stop("n and k must be of the same length (number of observations)\n")
    }
    prob_test <- matrix(nrow = length(k), ncol = length(p))
    for (j in seq_along(p))
        prob_test[, j] <- psi[j] * stats::dbinom(k, size = n, prob = p[j], 
                                                 log = FALSE) 
    prob_test <- prob_test / rowSums(prob_test)
    prob_test
}

