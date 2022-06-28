## This file is for binomial distributions

#' EM algorithm for estimating binomial mixture model
#'
#'
#' @param k A vector of integers. number of success
#' @param n A vector of integers. number of trials
#' @param n_components A number. number of components
#' @param p_init A vector of floats with length n_components, the initial value
#'   of p
#' @param learn_p bool(1) or a vector of bool, wheter learn each p
#' @param max_iter integer(1). number of maximum iterations
#' @param min_iter integer(1). number of minimum iterations
#' @param logLik_threshold A float. The threshold of logLikelihood increase for
#' detecting convergence
#'
#' @return a list containing \code{p}, a vector of floats between 0 and 1
#' giving the estimated success probability for each component, \code{psi},
#' estimated fraction of each component in the mixture, and \code{prob}, the
#' matrix of fitted probabilities of each observation belonging to each
#' component.
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
mixBinom <- function(k, n, n_components = 2, p_init = NULL, learn_p = TRUE,
    min_iter = 10, max_iter = 1000, logLik_threshold = 1e-5) {
    if (length(k) != length(n)) {
        stop("n and k must be of the same length (number of observations)\n")
    }
    S <- length(n)

    ## Random initialzation on parameters
    if (is.null(p_init)) {
        p <- stats::runif(n_components, 0.0, 1.0)
    } else {
        p <- p_init
    }
    psi <- stats::runif(n_components, 0.0, 1.0)
    psi <- psi / sum(psi)

    if (length(learn_p) == 1) {
        learn_p <- rep(learn_p, n_components)
    }

    logLik_mat <- matrix(0, nrow = S, ncol = n_components)

    ## Iterations
    for (it in seq_len(max_iter)) {
        ## E-step:
        RV <- predMixBinom(k, n, p = p, psi = psi)
        prob_mat <- RV$prob
        logLik_new <- RV$logLik

        ## M-step
        for (j in seq_len(n_components)) {
            if (learn_p[j]) {
                p[j] <- sum(prob_mat[, j] * k) / sum(prob_mat[, j] * n)
            }
            psi[j] <- sum(prob_mat[, j]) / S
        }

        ## Check convergence
        if (it > min_iter && logLik_new - logLik_old < logLik_threshold) {
            break
        } else {
            logLik_old <- logLik_new
        }
    }

    ## return values
    return_list <- list(
        "p" = p, "psi" = psi, "prob" = prob_mat,
        "logLik" = logLik_new
    )
    return_list
}


#' Predicted probability from learned binomial mixture model
#'
#' @param k A vector of integers. number of success
#' @param n A vector of integers. number of trials
#' @param p a vector of binomial success probabilities
#' @param psi A float between 0 and 1. fraction of each component
#'
#' @return A list with two components: prob, a matrix representing the 
#'   probability of each of the passed values coming from each component of the
#'   mixture and logLik, the total log-likelihood of the new samples.
#'
#' @export
#'
#' @examples
#' n1 <- array(sample(1:30, 50, replace = TRUE))
#' n2 <- array(sample(1:30, 200, replace = TRUE))
#' k1 <- apply(n1, 1, rbinom, n = 1, p = 0.5)
#' k2 <- apply(n2, 1, rbinom, n = 1, p = 0.01)
#' RV <- mixBinom(c(k1, k2), c(n1, n2))
#' RV_pred <- predMixBinom(3, 10, RV$p, RV$psi)
predMixBinom <- function(k, n, p, psi) {
    if (length(p) != length(psi)) {
        stop(
            "p and psi must be of the same length ",
            "(number of mixture components)\n"
        )
    }
    if (length(k) != length(n)) {
        stop("n and k must be of the same length (number of observations)\n")
    }

    ## assignment probability
    logLik_mat <- matrix(nrow = length(n), ncol = length(p))
    for (j in seq_len(length(p))) {
        logLik_mat[, j] <- log(psi[j]) + dbinom(k,
            size = n, prob = p[j],
            log = TRUE
        )
    }
    logLik_mat_amplify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
    prob_mat <- exp(logLik_mat_amplify) / rowSums(exp(logLik_mat_amplify))

    ## logLikelihood
    logLik_vec <- rep(NA, nrow(logLik_mat))
    for (i in seq_len(nrow(logLik_mat))) {
        logLik_vec[i] <- matrixStats::logSumExp(logLik_mat[i, ], na.rm = TRUE)
    }
    logLik_val <- sum(logLik_vec, na.rm = TRUE)

    list("prob" = prob_mat, "logLik" = logLik_val)
}
