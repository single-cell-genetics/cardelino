# Functions for clone ID

#' Infer clonal identity of single cells
#'
#' @param A variant x cell matrix of integers; number of alternative allele
#' reads in variant i cell j
#' @param D variant x cell matrix of integers; number of total reads covering
#' variant i cell j
#' @param Config variant x clone matrix of binary values. The clone-variant
#' configuration, which encodes the phylogenetic tree structure. This is the
#' output Z of Canopy
#' @param n_clone integer(1), the number of clone to reconstruct. This is in use
#' only if Config is NULL
#' @param relax_Config logical(1), If TRUE, relaxing the Clone Configuration by
#' changing it from fixed value to act as a prior Config with a relax rate.
#' @param relax_rate_fixed numeric(1), If the value is between 0 to 1,
#' the relax rate will be set as a fix value during updating clone Config. If
#' NULL, the relax rate will be learned automatically with relax_rate_prior.
#' @param n_chain integer(1), the number of chains to run, which will be
#' averaged as an output result
#' @param n_proc integer(1), the number of processors to use. This parallel
#' computing can largely reduce time when using multiple chains
#' @param inference character(1), the method to use for inference, either
#' "sampling" to use Gibbs sampling (default) or "EM" to use
#' expectation-maximization (faster)
#' @param verbose logical(1), should the function output verbose information as
#' it runs?
#' @param ... arguments passed to \code{\link{clone_id_Gibbs}} or
#' \code{\link{clone_id_EM}} (as appropriate)
#' @param Psi A vector of float. The fractions of each clone, output P of Canopy
#' @param min_iter A integer. The minimum number of iterations in the Gibbs
#' sampling. The real iteration may be longer until the convergence.
#' @param max_iter A integer. The maximum number of iterations in the Gibbs
#' sampling, even haven't passed the convergence diagnosis
#'
#' @return
#' If inference method is "EM", a list containing \code{theta}, a vector of
#' two floats denoting the parameters of the two components of the base model,
#' i.e., mean of Bernoulli or binomial model given variant exists or not,
#' \code{prob}, the matrix of posterior probabilities of each cell belonging to
#' each clone with fitted parameters, and \code{logLik}, the log likelihood of
#' the final parameters.
#'
#' If inference method is "sampling", a list containing: \code{theta0}, the mean
#' of sampled false positive parameter values; \code{theta1} the mean of sampled
#' (1 - false negative rate) parameter values; \code{theta0_all}, all sampled
#' false positive parameter values; \code{theta1_all}, all sampled (1 - false
#' negative rate) parameter values; \code{element}; \code{logLik_all},
#' log-likelihood for model for all sampled parameter sets; \code{prob_all};
#' \code{prob}, matrix with mean of sampled cell-clone assignment posterior
#' probabilities (the key output of the model); \code{prob_variant}.
#'
#' @details
#' The two Bernoulli components correspond to false positive and false negative
#' rates. The two binomial components correspond to the read distributions
#' with and without the mutation present.
#'
#' @author Yuanhua Huang and Davis McCarthy
#'
#' @import stats
#' @name Clone ID
#' @rdname clone_id
#' @export
#'
#' @examples
#' data(example_donor)
#' assignments <- clone_id(A_clone, D_clone,
#'     Config = tree$Z,
#'     min_iter = 800, max_iter = 1200
#' )
#' prob_heatmap(assignments$prob)
#'
#' assignments_EM <- clone_id(A_clone, D_clone,
#'     Config = tree$Z,
#'     inference = "EM"
#' )
#' prob_heatmap(assignments_EM$prob)
clone_id <- function(A, D, Config = NULL, n_clone = NULL, Psi = NULL,
    relax_Config = TRUE, relax_rate_fixed = NULL,
    inference = "sampling", n_chain = 1, n_proc = 1,
    verbose = TRUE, ...) {
    inference <- match.arg(inference, c("sampling", "EM"))
    ## check input data
    if (!(all(rownames(A) == rownames(D)))) {
        stop("Rownames for A and D are not identical.")
    }
    if (!(all(colnames(A) == colnames(D)))) {
        stop("Colnames for A and D are not identical.")
    }
    if (is.null(Config) && is.null(n_clone)) {
        stop("Config and n_clone can't be NULL together.")
    }
    if (any(A %% 1 != 0, na.rm = TRUE) || 
        any(D %% 1 != 0, na.rm = TRUE)) {
        stop("A and/or D contain non-integer values.")
    }
    ## Cardelino-free mode if Config is NULL
    if (is.null(Config)) {
        cat("Config is NULL: de-novo mode is in use.\n")
        Config <- matrix(0, nrow = nrow(D), ncol = n_clone)
        rownames(Config) <- rownames(D)
        colnames(Config) <- paste0("Clone", seq_len(n_clone))
        relax_Config <- TRUE
        relax_rate_fixed <- 0.5
    }

    ## Match exome-seq and scRNA-seq data
    if (!any(rownames(D) %in% rownames(Config))) {
        stop("No matches in variant names between Config and D arguments.")
    }
    ## match variants
    common_vars <- intersect(rownames(Config), rownames(D))
    A <- A[common_vars, , drop = FALSE]
    D <- D[common_vars, , drop = FALSE]
    Config <- Config[common_vars, , drop = FALSE]
    if (verbose) {
        message(length(common_vars), " variants used for cell assignment.")
    }

    ## change sparse matrix to dense matrix
    A <- as.matrix(A)
    D <- as.matrix(D)

    ## pass data to specific functions
    if (inference == "sampling") {
        # Try to run sampling in parallel. 
        if (n_proc > 1) {
            
            # doMC is not available on Windows so there we fall back to 
            # running the chains using only one processor. 
            if (!requireNamespace("doMC", quietly = TRUE)) {
                message("Parallel sampling is not supported on Windows.")
                message("Falling back to sampling using one processor.")
                ids_list <- list(
                    clone_id_Gibbs(A, D, Config,
                                   Psi = Psi,
                                   relax_Config = relax_Config,
                                   relax_rate_fixed = relax_rate_fixed,
                                   verbose = verbose, ...
                    )
                )     
            } else {
                doMC::registerDoMC(n_proc)
                `%dopar%` <- foreach::`%dopar%`
                
                ids_list <- foreach::foreach(ii = seq_len(n_chain)) %dopar% {
                    clone_id_Gibbs(A, D, Config,
                                   Psi = Psi,
                                   relax_Config = relax_Config,
                                   relax_rate_fixed = relax_rate_fixed,
                                   verbose = verbose, ...
                    )
                }
            }
        } else {
            ids_list <- list(
                clone_id_Gibbs(A, D, Config,
                    Psi = Psi,
                    relax_Config = relax_Config,
                    relax_rate_fixed = relax_rate_fixed,
                    verbose = verbose, ...
                )
            )
        }

        ids_out <- ids_list[[1]]
        ids_out$n_chain <- 1
        if (n_chain > 1) {
            for (ii in seq(2, n_chain)) {
                ids_out$n_chain <- ids_out$n_chain + 1
                idx <- colMatch(ids_out$prob, ids_list[[ii]]$prob, force = TRUE)
                ids_out$prob <- ids_out$prob + ids_list[[ii]]$prob[, idx]
                ids_out$relax_rate <- ids_out$relax_rate +
                    ids_list[[ii]]$relax_rate
                ids_out$Config_prob <- (ids_out$Config_prob +
                    ids_list[[ii]]$Config_prob[, idx])
            }
            ids_out$prob <- ids_out$prob / n_chain
            ids_out$relax_rate <- ids_out$relax_rate / n_chain
            ids_out$Config_prob <- ids_out$Config_prob / n_chain
        }
        return(ids_out)
    }
    else {
        return(clone_id_EM(A, D, Config, verbose = verbose, ...))
    }
}


#' Assign cells to clones from cardelino results
#'
#' @param prob_mat numeric matrix (cells x clones) of clone posterior
#' probabilities as output by \code{\link{clone_id}}
#' @param threshold numeric(1), posterior probability threshold for cell-clone
#' assignment: if posterior probability is above threshold, assign cell to
#' clone, otherwise leave cell "unassigned"
#'
#' @return a \code{data.frame} with cell ID, assigned clone label and maximum
#' posterior probability across clones.
#'
#' @author Davis McCarthy
#' @export
#' @examples
#' data(example_donor)
#' assignments <- clone_id(A_clone, D_clone, Config = tree$Z, inference = "EM")
#' df <- assign_cells_to_clones(assignments$prob)
#' head(df)
#' table(df$clone)
assign_cells_to_clones <- function(prob_mat, threshold = 0.5) {
    assigns <- data.frame(
        cell = rownames(prob_mat),
        clone = rep("unassigned", nrow(prob_mat)),
        prob_max = rep(NA, nrow(prob_mat)),
        stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(prob_mat))) {
        assigns[i, "prob_max"] <- max(prob_mat[i, ])
        if (max(prob_mat[i, ]) > threshold) {
            assigns[i, "clone"] <- names(which.max(prob_mat[i, ]))
        }
    }
    assigns
}


#' Assign cells to clones using an EM algorithm
#'
#' @param logLik_threshold A float. The threshold of logLikelihood increase for
#' detecting convergence.
#'
#' @return a list containing \code{theta}, a vector of two floats denoting the
#' binomial rates given variant exists or not, \code{prob}, the matrix of
#' posterior probabilities of each cell belonging to each clone with fitted
#' parameters, and \code{logLik}, the log likelihood of the final parameters.
#'
#' @author Yuanhua Huang
#' @import stats
#' @rdname clone_id
#' @export
#'
clone_id_EM <- function(A, D, Config, Psi = NULL, min_iter = 10,
    max_iter = 1000, logLik_threshold = 1e-5,
    verbose = TRUE) {
    if (is.null(Psi)) {
        Psi <- rep(1 / ncol(Config), ncol(Config))
    }
    if (dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2] ||
        dim(A)[1] != dim(Config)[1] || dim(Config)[2] != length(Psi)) {
        stop(
            "A and D must have the same size;\n",
            "A and Config must have the same number of variants;\n",
            "Config and Psi must have the same number of clones"
        )
    }

    ## preprocessing
    N <- dim(A)[1] # number of variants
    M <- dim(A)[2] # number of cells
    K <- dim(Config)[2] # number of clones
    A[which(D == 0)] <- NA
    D[which(D == 0)] <- NA
    A[(D > 0) & is.na(A)] <- 0

    C1 <- Config
    C0 <- 1 - Config
    A1 <- A # number of alteration reads
    B1 <- D - A # number of reference reads
    W_log <- sum(lchoose(D, A), na.rm = TRUE) # log binomial coefficients

    A1[is.na(A1)] <- 0
    B1[is.na(B1)] <- 0

    S1 <- t(A1) %*% C0 # Alteration reads without variants (C=0): theta[1]
    S2 <- t(B1) %*% C0 # Reference reads without variants  (C=0): 1-theta[1]
    S3 <- t(A1) %*% C1 # Alteration reads with variants (C=1): theta[2]
    S4 <- t(B1) %*% C1 # Reference reads with variants  (C=1): 1-theta[2]

    ## random initialization for EM
    theta <- c(stats::runif(1, 0.001, 0.25), stats::runif(1, 0.3, 0.6))
    logLik <- 0
    logLik_new <- logLik + logLik_threshold * 2

    ## EM iterations
    for (it in seq_len(max_iter)) {
        # Check convergence
        if ((it > min_iter) && ((logLik_new - logLik) < logLik_threshold)) {
            break
        }
        logLik <- logLik_new

        # E-step
        logLik_mat <- (S1 * log(theta[1]) + S2 * log(1 - theta[1]) +
            S3 * log(theta[2]) + S4 * log(1 - theta[2]))
        logLik_mat <- t(t(logLik_mat) + log(Psi))

        logLik_vec <- rep(NA, nrow(logLik_mat))
        for (i in seq_len(nrow(logLik_mat))) {
            logLik_vec[i] <- matrixStats::logSumExp(logLik_mat[i, ],
                na.rm = TRUE
            )
        }
        logLik_new <- sum(logLik_vec, na.rm = TRUE) + W_log
        logLik_mat_amplify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
        prob_mat <- exp(logLik_mat_amplify) / rowSums(exp(logLik_mat_amplify))

        # M-step
        if (is.na(sum(prob_mat * (S1 + S2))) ||
            sum(prob_mat * (S1 + S2)) == 0) {
            theta[1] <- 0.02
        } else {
            theta[1] <- sum(prob_mat * S1) / sum(prob_mat * (S1 + S2))
        }
        if (is.na(sum(prob_mat * (S3 + S4))) ||
            sum(prob_mat * (S3 + S4)) == 0) {
            theta[2] <- 0.75
        } else {
            theta[2] <- sum(prob_mat * S3) / sum(prob_mat * (S3 + S4))
        }
    }
    if (verbose) {
        print(paste("Total iterations:", it))
    }

    ## return values
    return_list <- list("theta" = theta, "prob" = prob_mat, "logLik" = logLik)
    return_list
}


#' Assign cells to clones using a Gibbs sampling algorithm
#'
#' @param prior0 numeric(2), alpha and beta parameters for the Beta prior
#' distribution on the inferred false positive rate.
#' @param prior1 numeric(2), alpha and beta parameters for the Beta prior
#' distribution on the inferred (1 - false negative) rate.
#' @param relax_rate_prior numeric(2), the two parameters of beta prior
#' distribution of the relax rate for relaxing the clone Configuration. This
#' mode is used when relax_relax is NULL.
#' @param keep_base_clone bool(1), if TRUE, keep the base clone of Config to its
#' input values when relax mode is used.
#' @param buin_frac numeric(1), the fraction of chain as burn-in period
#' @param wise A string, the wise of parameters for theta1: global, variant,
#' element.
#' @param relabel bool(1), if TRUE, relabel the samples of both Config and prob
#' during the Gibbs sampling.
#'
#' @author Yuanhua Huang
#' @import matrixStats
#' @rdname clone_id
#' @export
#'
clone_id_Gibbs <- function(A, D, Config, Psi = NULL,
    relax_Config = TRUE, relax_rate_fixed = NULL,
    relax_rate_prior = c(1, 9), keep_base_clone = TRUE,
    prior0 = c(0.2, 99.8), prior1 = c(0.45, 0.55),
    min_iter = 5000, max_iter = 20000, buin_frac = 0.5,
    wise = "variant", relabel = FALSE, verbose = TRUE) {
    if (is.null(Psi)) {
        Psi <- rep(1 / ncol(Config), ncol(Config))
    }
    if (dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2] ||
        dim(A)[1] != dim(Config)[1] || dim(Config)[2] != length(Psi)) {
        stop(
            "A and D must have the same size;\n ",
            "A and Config must have the same number of variants;\n",
            "Config and Psi must have the same number of clones"
        )
    }
    if (sum(c("element", "variant", "global") == wise) == 0) {
        stop(
            "Input wise mode: ", wise,
            ", while only supporting: element, variant, global"
        )
    }

    ## preprocessing
    N <- dim(A)[1] # number of variants
    M <- dim(A)[2] # number of cells
    K <- dim(Config)[2] # number of clones

    A[which(D == 0)] <- NA
    D[which(D == 0)] <- NA
    A[(D > 0) & is.na(A)] <- 0

    C1 <- Config
    C0 <- 1 - Config
    A1 <- A # number of alteration reads
    B1 <- D - A # number of reference reads
    W_log <- sum(lchoose(D, A), na.rm = TRUE) # log binomial coefficients

    A1[is.na(A1)] <- 0
    B1[is.na(B1)] <- 0

    # reads number list for each clone
    S1_list <- list()
    S2_list <- list()
    S3_list <- list()
    S4_list <- list()
    for (k in seq_len(K)) {
        S1_list[[k]] <- A1 * C0[, k]
        S2_list[[k]] <- B1 * C0[, k]
        S3_list[[k]] <- A1 * C1[, k]
        S4_list[[k]] <- B1 * C1[, k]
    }

    ## Prepare for sampling
    if (wise == "global") {
        idx_vec <- seq_len(1) # For: theta1_all[t, ] <- theta1[idx_vec]
        idx_mat <- seq_len(N * M) # For: theta1[idx_mat] <- theta1_all[t, ]
    } else if (wise == "variant") {
        idx_vec <- seq_len(N)
        idx_mat <- seq_len(N * M)
    } else if (wise == "element") {
        idx_vec <- which(A1 + B1 > 0)
        idx_mat <- which(A1 + B1 > 0)
    }
    n_element <- length(idx_vec)

    if (is.null(dim(prior1)) && length(prior1) == 2) {
        # two variable to a matrix
        prior1 <- t(matrix(rep(prior1, n_element), nrow = 2))
    }
    if (!is.matrix(prior1)) {
        stop("prior1 need to be a matrix of n_element x 2")
    }

    prob_all <- matrix(0, nrow = max_iter, ncol = M * K)
    logLik_mat <- matrix(0, nrow = M, ncol = K)
    logLik_all <- matrix(0, nrow = max_iter, ncol = 1)
    assign_all <- matrix(0, nrow = max_iter, ncol = M)
    theta0_all <- matrix(0, nrow = max_iter, ncol = 1)
    theta1_all <- matrix(0, nrow = max_iter, ncol = n_element)
    Config_all <- matrix(0, nrow = max_iter, ncol = N * K)
    relax_rate_all <- matrix(0, nrow = max_iter, ncol = 1)

    if (!is.null(relax_Config) && relax_Config != FALSE) {
        if (!is.null(relax_rate_fixed)) {
            if (relax_rate_fixed > 1 || relax_rate_fixed < 0) {
                stop("relax_rate_fixed needs to be NULL or in [0, 1].")
            }
            relax_rate <- relax_rate_fixed ## fixed relax_rate
            relax_rate_all[, ] <- relax_rate
        } else if (!is.null(relax_rate_prior)) {
            relax_rate <- relax_rate_prior[1] / (relax_rate_prior[1] +
                relax_rate_prior[2])
        } else {
            stop("Require value for either relax_Config or relax_prior.")
        }

        Config_new <- Config
        Config_prior <- Config
        Config_prior[Config == 1] <- 1 - relax_rate
        Config_prior[Config == 0] <- relax_rate
        if (keep_base_clone) {
            Config_prior[, 1] <- Config[, 1]
        }
        Config_prior_oddlog <- log(Config_prior) - log(1 - Config_prior)
        Iden_mat <- matrix(0, nrow = M, ncol = K)
    }

    ## Random initialization
    theta0_all[1, 1] <- stats::rbeta(1, prior0[1], prior0[2])
    theta1_all[1, ] <- stats::rbeta(rep(1, n_element), prior1[, 1], prior1[, 2])
    theta1 <- matrix(NA, nrow = N, ncol = M)

    ## Gibbs sampling
    for (it in 2:max_iter) {
        # Update prob_mat
        theta0 <- theta0_all[it - 1, 1]
        theta1[idx_mat] <- theta1_all[it - 1, ]
        for (k in seq_len(K)) {
            logLik_mat[, k] <- (colSums(S1_list[[k]] * log(theta0),
                na.rm = TRUE
            ) +
                colSums(S2_list[[k]] * log(1 - theta0), na.rm = TRUE) +
                colSums(S3_list[[k]] * log(theta1), na.rm = TRUE) +
                colSums(S4_list[[k]] * log(1 - theta1), na.rm = TRUE))
            logLik_mat[, k] <- logLik_mat[, k] + log(Psi[k])
        }
        logLik_mat_amplify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
        prob_mat <- exp(logLik_mat_amplify) / rowSums(exp(logLik_mat_amplify))
        prob_all[it, ] <- prob_mat

        # logLik_vec <- rep(NA, nrow(logLik_mat))
        # for (i in seq_len(nrow(logLik_mat))) {
        #     logLik_vec[i] <- matrixStats::logSumExp(logLik_mat[i,],
        #                                             na.rm = TRUE)
        # }
        # logLik_all[it] <- sum(logLik_vec, na.rm = TRUE) + W_log
        # #Update logLikelihood (TODO: add W0 and W1)

        # Sample assignment
        for (j in seq_len(M)) {
            assign_all[it, j] <- sample(seq_len(K), 1,
                replace = TRUE,
                prob = prob_mat[j, ]
            )
        }

        ## Update Config
        if (!is.null(relax_Config) && relax_Config != FALSE) {
            if (it > (0.1 * min_iter + 5) && is.null(relax_rate_fixed)) {
                diff0 <- sum((Config == Config_new)[, 2:ncol(Config)])
                diff1 <- sum((Config != Config_new)[, 2:ncol(Config)])
                relax_rate <- stats::rbeta(
                    1, relax_rate_prior[1] + diff1,
                    relax_rate_prior[2] + diff0
                )
                relax_rate_all[it] <- relax_rate

                Config_prior <- Config
                Config_prior[Config == 1] <- 1 - relax_rate
                Config_prior[Config == 0] <- relax_rate
                if (keep_base_clone) {
                    Config_prior[, 1] <- Config[, 1]
                }
                Config_prior_oddlog <- log(Config_prior) - log(1 - Config_prior)
            }

            Iden_mat[, ] <- 0
            for (j in seq_len(M)) {
                Iden_mat[j, assign_all[it, j]] <- 1
            }

            # calculate log_probability matrix with genotype 0 and 1
            P0_mat <- A1 * log(theta0) + B1 * log(1 - theta0) + W_log
            P1_mat <- A1 * log(theta1) + B1 * log(1 - theta1) + W_log

            oddR_log <- P1_mat %*% Iden_mat - P0_mat %*% Iden_mat
            oddR_log <- oddR_log + Config_prior_oddlog
            oddR_log[which(oddR_log > 50)] <- 50
            oddR_log[which(oddR_log < -50)] <- -50
            Config_prob_tmp <- exp(oddR_log) / (exp(oddR_log) + 1)

            Config_new[, ] <- stats::rbinom(N * K, size = 1, Config_prob_tmp)
            Config_all[it, ] <- Config_new

            for (k in seq_len(K)) {
                S1_list[[k]] <- A1 * (1 - Config_new[, k])
                S2_list[[k]] <- B1 * (1 - Config_new[, k])
                S3_list[[k]] <- A1 * Config_new[, k]
                S4_list[[k]] <- B1 * Config_new[, k]
            }
        }

        # Sample theta with assigned clones
        S1_wgt <- S2_wgt <- 0 # weighted S1
        S3_wgt <- S4_wgt <- matrix(0, nrow = N, ncol = M)
        for (k in seq_len(K)) {
            idx <- which(assign_all[it, ] == k)
            S1_wgt <- S1_wgt + sum(S1_list[[k]][, idx], na.rm = TRUE)
            S2_wgt <- S2_wgt + sum(S2_list[[k]][, idx], na.rm = TRUE)
            S3_wgt[, idx] <- S3_wgt[, idx] + S3_list[[k]][, idx]
            S4_wgt[, idx] <- S4_wgt[, idx] + S4_list[[k]][, idx]
        }

        if (wise == "global") {
            S3_wgt[, ] <- sum(S3_wgt, na.rm = TRUE)
            S4_wgt[, ] <- sum(S4_wgt, na.rm = TRUE)
        } else if (wise == "variant") {
            S3_wgt[, ] <- rowSums(S3_wgt, na.rm = TRUE)
            S4_wgt[, ] <- rowSums(S4_wgt, na.rm = TRUE)
        }
        theta0_all[it, 1] <- stats::rbeta(
            1, prior0[1] + S1_wgt,
            prior0[2] + S2_wgt
        )
        theta1_all[it, ] <- stats::rbeta(
            rep(1, n_element),
            prior1[, 1] + S3_wgt[idx_vec],
            prior1[, 2] + S4_wgt[idx_vec]
        )

        # Calculate logLikelihood
        logLik_all[it] <- get_logLik(
            A1, B1, Config_new, assign_all[it, ],
            theta0, theta1
        )

        # Check convergence. TODO: consider relabel
        if ((it >= min_iter) && (it %% 100 == 0)) {
            Converged_all <- abs(Geweke_Z(prob_all[seq_len(it), ])) <= 2
            # Converged_all <- abs(Geweke_Z(Config_all[1:it, ])) <= 2
            if (verbose) {
                cat(paste0(
                    round(mean(Converged_all, na.rm = TRUE), 3) * 100,
                    "% converged.\n"
                ))
            }
            if (mean(Converged_all, na.rm = TRUE) > 0.995) {
                break
            }
        }
    }
    Converged_all <- abs(Geweke_Z(prob_all[seq_len(it), ])) <= 2
    percent_converged <- round(mean(Converged_all, na.rm = TRUE), 3) * 100
    if (percent_converged <= 95) {
        warning("Only ", percent_converged, "% of parameters converged. ", 
                "Consider increasing `max_iter`.", call. = FALSE)
    } else {
        print(paste("Converged in", it, "iterations.")) 
    }

    ## Return values
    n_buin <- ceiling(it * buin_frac)

    a <- A1[idx_mat]
    d <- A1[idx_mat] + B1[idx_mat]
    binom_pdf1 <- binom_pdf0 <- rep(0, n_element)
    for (i in seq(n_buin, it)) {
        binom_pdf1 <- binom_pdf1 + stats::dbinom(a,
            size = d,
            prob = theta1_all[i, ]
        )
        binom_pdf0 <- binom_pdf0 + stats::dbinom(a,
            size = d,
            prob = theta0_all[i]
        )
    }
    prob_variant <- matrix(NA, nrow = N, ncol = M)
    prob_variant[idx_mat] <- binom_pdf1 / (binom_pdf1 + binom_pdf0)
    row.names(prob_variant) <- row.names(A)
    colnames(prob_variant) <- colnames(A)

    if (relabel) {
        col_idx_use <- seq(K)
        for (ii in seq(n_buin, it)) {
            mat1 <- matrix(prob_all[ii - 1, ], nrow = M)
            mat2 <- matrix(prob_all[ii, ], nrow = M)
            # mat1 <- matrix(Config_all[ii - 1, ], nrow = N)
            # mat2 <- matrix(Config_all[ii, ], nrow = N)

            if (ncol(mat1) <= 5) {
                idx <- colMatch(mat1, mat2, force = TRUE)
            }
            else {
                idx <- colMatch(mat1, mat2, force = FALSE)
            }
            col_idx_use <- col_idx_use[idx]
            prob_all[ii, ] <- matrix(prob_all[ii, ], nrow = M)[, col_idx_use]
            Config_all[ii, ] <- matrix(Config_all[ii, ],
                nrow = N
            )[, col_idx_use]
        }
    }
    prob_mat <- matrix(colMeans(prob_all[n_buin:it, ]), nrow = M)
    row.names(prob_mat) <- colnames(A)
    colnames(prob_mat) <- colnames(Config)

    Config_prob <- Config
    Config_prob[, ] <- colMeans(Config_all[n_buin:it, ])

    theta0 <- mean(theta0_all[n_buin:it, ])
    theta1[idx_mat] <- colMeans(as.matrix(theta1_all[n_buin:it, ]))


    logLik_post <- get_logLik(A1, B1, Config_prob, prob_mat, theta0, theta1)
    DIC <- devianceIC(logLik_all[n_buin:it], logLik_post)

    #     DIC <- devianceIC(logLik_all[n_buin:it],
    #                       A1, B1, Config_prob, theta0, theta1, Psi, W_log)

    return_list <- list(
        "theta0" = theta0, "theta1" = theta1,
        "theta0_all" = as.matrix(theta0_all[seq_len(it), ]),
        "theta1_all" = as.matrix(theta1_all[seq_len(it), ]),
        "element" = idx_mat, "logLik" = logLik_all[seq_len(it)],
        "prob_all" = prob_all[seq_len(it), ],
        "prob" = prob_mat, "prob_variant" = prob_variant,
        "relax_rate" = mean(relax_rate_all[n_buin:it]),
        "Config_prob" = Config_prob,
        "Config_all" = Config_all[seq_len(it), ],
        "relax_rate_all" = relax_rate_all[seq_len(it)], "DIC" = DIC
    )
    return_list
}


#' Geweke diagnostic for MCMC sampling.
#'
#' @param X A matrix of MCMC samples for N samples per K variables
#' @param first A float between 0 and 1. The initial region of MCMC chain.
#' @param last A float between 0 and 1. The final region of MCMC chain.
#'
#' @author Yuanhua Huang
#' @return \code{Z}, a vector of absolute value of Z scores for each variable.
#' When |Z| <= 2, the sampling could be taken as converged.
#'
Geweke_Z <- function(X, first = 0.1, last = 0.5) {
    # The original Geweke’s diagnostics uses pectral densities to estimate
    # the sample variances (Geweke,1992), for adjusting the samples
    # as the two ‘samples’ are not independent
    # so as the implementation in coda::geweke.diag()

    # Converged_all = rep(NA, ncol(prob_all))
    # for (ic in seq_len(length(Converged_all))) {
    # Converged_all[ic] = abs(coda::geweke.diag(prob_all[, ic])$z) <= 2
    # }

    if (is.null(dim(X))) {
        X <- as.matrix(X, ncol = 1)
    }
    N <- nrow(X)
    A <- X[seq_len(floor(first * N)), , drop = FALSE]
    B <- X[ceiling(last * N):N, , drop = FALSE]

    A_col_mean <- colMeans(A)
    B_col_mean <- colMeans(B)
    A_col_var <- rowSums((t(A) - A_col_mean)^2) / (nrow(A) - 1) # ^2
    B_col_var <- rowSums((t(B) - B_col_mean)^2) / (nrow(A) - 1) # ^2

    min_var <- 10^(-50)
    Z <- (A_col_mean - B_col_mean) / sqrt(A_col_var + B_col_var + min_var)

    Z
}

#' Log likelihood of clone_id model
#' It returns P(A, D | C, I, theta0, theta1)
#'
#' @param A1 variant x cell matrix of integers; number of alternative allele
#' reads in variant i cell j
#' @param B1 variant x cell matrix of integers; number of reference allele
#' reads in variant i cell j
#' @param Config variant x clone matrix of float values. The clone-variant
#' configuration probability, averaged by posterior samples
#' @param Assign cells x clone matrix of float values. The cell-clone
#' assignment probability, averaged by posterior samples
#' @param theta0 the binomial rate for alternative allele from config = 0
#' @param theta1 the binomial rate for alternative allele from config = 1
#'
#' @author Yuanhua Huang
#' @return \code{logLik}, a float of log likelihood
#'
get_logLik <- function(A1, B1, Config, Assign, theta0, theta1) {
    if (is.null(dim(Assign)) || length(dim(Assign)) == 1) {
        Assign_prob <- matrix(0, length(Assign), ncol(Config))
        for (i in seq_len(length(Assign))) {
            Assign_prob[i, Assign[i]] <- 1
        }
    } else {
        Assign_prob <- Assign
    }

    prob_mat <- Config %*% t(Assign_prob)

    Lik_mat <- (exp(log(theta1) * A1 + log(1 - theta1) * B1) * prob_mat +
        exp(log(theta0) * A1 + log(1 - theta0) * B1) * (1 - prob_mat))

    logLik <- (sum(log(Lik_mat), na.rm = TRUE) +
        sum(lchoose(A1 + B1, A1), na.rm = TRUE))
    logLik
}


#' Deviance Information Criterion for cardelino model
#'
#' @param logLik_all A vector of numeric; the log likelihood of posterior
#' sample, i.e., posterior samples of deviance
#' @param logLik_post numeric(1); the log likelihood of mean posterior
#' parameters, i.e., deviance of posterior means
#'
#' @author Yuanhua Huang
#' @return \code{DIC}, a float of deviance information criterion
#'
devianceIC <- function(logLik_all, logLik_post) {
    # https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/DIC-slides.pdf
    # Spiegelhalter et al. The deviance information criterion: 12 years on, 2014
    # Spiegelhalter et al. Bayesian measures of model complexity and fit, 2002
    # Gelman et al. Bayesian Data Analysis. 3rd Edition, 2013
    # (Charpter 7.2, p173)

    logLik_mean <- mean(logLik_all)
    logLik_var <- var(logLik_all)

    p_D_Spiegelhalter <- -2 * logLik_mean - (-2 * logLik_post)
    DIC_Spiegelhalter <- -2 * logLik_post + 2 * p_D_Spiegelhalter

    p_D_Gelman <- 2 * logLik_var
    DIC_Gelman <- -2 * logLik_post + 2 * p_D_Gelman

    DIC <- DIC_Gelman

    cat(paste(
        "DIC:", round(DIC, 2),
        "D_mean:", round(-2 * logLik_mean, 2),
        "D_post:", round(-2 * logLik_post, 2),
        "logLik_var:", round(logLik_var, 2), "\n"
    ))

    list(
        "DIC" = DIC,
        "logLik_var" = logLik_var,
        "D_mean" = -2 * logLik_mean,
        "D_post" = -2 * logLik_post,
        "DIC_Gelman" = DIC_Gelman,
        "DIC_Spiegelhalter" = DIC_Spiegelhalter
    )
}
