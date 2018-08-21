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
#' @param model character(1), the model to use for inference; either "binomial"
#' to use the beta-binomial mixture model or "Bernoulli" to use the Bernoulli
#' mixure model (default: "binomial")
#' @param inference character(1), the method to use for inference, either
#' "sampling" to use Gibbs sampling (default) or "EM" to use
#' expectation-maximisation (faster)
#' @param verbose logical(1), should the function output verbose information as
#' it runs?
#' @param ... arguments passed to \code{\link{cell_assign_Gibbs}} or
#' \code{\link{cell_assign_EM}} (as appropriate)
#' @param Psi A vector of float. The fractions of each clone, output P of Canopy
#' @param min_iter A integer. The minimum number of iterations in the Gibbs
#' sampling. The real iteration may be longer utile the convergence.
#' @param max_iter A integer. The maximum number of iterations in the Gibbs
#' sampling, even haven't passed the convergence diagnosis
#'
#' @return
#' If inference method is "EM", a list containing \code{theta}, a vector of
#' two floats denoting the parameters of the two componets of the base model,
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
#' assignments <- clone_id(A_clone, D_clone, Config = tree$Z)
#' prob_heatmap(assignments$prob)
#'
#' assignments_EM <- clone_id(A_clone, D_clone, Config = tree$Z, inference = "EM")
#' prob_heatmap(assignments_EM$prob)
#'
#' assignments_bern <- clone_id(A_clone, D_clone, Config = tree$Z, model = "Bernoulli")
#' prob_heatmap(assignments_bern$prob)
#'
#' assignments_bern_EM <- clone_id(A_clone, D_clone, Config = tree$Z, model = "Bernoulli",
#'                                  inference = "EM")
#' prob_heatmap(assignments_bern_EM$prob)
#'
clone_id <- function(A, D, Config, model = "binomial", inference = "sampling",
                     verbose = TRUE, ...) {
    model <- match.arg(model, c("binomial", "Bernoulli"))
    inference <- match.arg(inference, c("sampling", "EM"))
    ## check input data
    if (!(all(rownames(A) == rownames(D))))
        stop("Rownames for A and D are not identical.")
    if (!(all(colnames(A) == colnames(D))))
        stop("Colnames for A and D are not identical.")
    ## Match exome-seq and scRNA-seq data
    if (!any(rownames(D) %in% rownames(Config)))
        stop("No matches in variant names between Config and D arguments.")
    ## match variants
    common_vars <- intersect(rownames(Config), rownames(D))
    A <- A[common_vars,, drop = FALSE]
    D <- D[common_vars,, drop = FALSE]
    Config <- Config[common_vars,, drop = FALSE]
    if (verbose)
        message(length(common_vars), " variants used for cell assignment.")
    ## pass data to specific functions
    if (inference == "sampling")
        return(cell_assign_Gibbs(A, D, Config, model = model, verbose = verbose,
                                 ...))
    else
        return(cell_assign_EM(A, D, Config, model = model, verbose = verbose,
                              ...))
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
#'
#' @export
#'
#' @examples
#' data(example_donor)
#' assignments <- clone_id(A_clone, D_clone, Config = tree$Z, inference = "EM")
#' df <- assign_cells_to_clones(assignments$prob)
#' head(df)
#' table(df$clone)
#'
assign_cells_to_clones <- function(prob_mat, threshold = 0.5) {
    assigns <- data.frame(cell = rownames(prob_mat),
                          clone = rep("unassigned", nrow(prob_mat)),
                          prob_max = rep(NA, nrow(prob_mat)),
                          stringsAsFactors = FALSE)
    for (i in seq_len(nrow(prob_mat))) {
        assigns[i, "prob_max"] <- max(prob_mat[i,])
        if (max(prob_mat[i,]) > threshold) {
            assigns[i, "clone"] <- names(which.max(prob_mat[i,]))
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
#' parameters of the two componets of the base model, i.e., mean of Bernoulli or
#' binomial model given variant exists or not, \code{prob}, the matrix of
#' posterior probabilities of each cell belonging to each clone with fitted
#' parameters, and \code{logLik}, the log likelihood of the final parameters.
#'
#' @import stats
#' @rdname clone_id
#' @export
#'
#' @examples
#' assignments_binom <- cell_assign_EM(A_clone, D_clone, Config = tree$Z, model = "binomial")
#' assignments_bern <- cell_assign_EM(A_clone, D_clone, Config = tree$Z, model = "Bernoulli")
#'
cell_assign_EM <- function(A, D, Config, Psi=NULL, model="binomial",
                           Bernoulli_threshold=1, min_iter=10, max_iter=1000,
                           logLik_threshold=1e-5, verbose=TRUE) {
    if (is.null(Psi))
        Psi <- rep(1/ncol(Config), ncol(Config))
    if (dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2] ||
        dim(A)[1] != dim(Config)[1] || dim(Config)[2] != length(Psi)) {
        stop(paste("A and D must have the same size;\n ",
                   "A and Config must have the same number of variants;\n",
                   "Config and Psi must have the same number of clones",sep = ""))
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
    if (model == "Bernoulli") {
        A1 <- A >= Bernoulli_threshold
        B1 <- 1 - A1
        W_log <- 0
    } else{
        A1 <- A                #number of alteration reads
        B1 <- D - A            #number of reference reads
        W_log <- sum(lchoose(D, A), na.rm = TRUE)  #log binomial coefficients
    }
    A1[is.na(A1)] <- 0
    B1[is.na(B1)] <- 0

    S1 <- t(A1) %*% C0  #Alteration reads without variants (C=0): theta[1]
    S2 <- t(B1) %*% C0  #Reference reads without variants  (C=0): 1-theta[1]
    S3 <- t(A1) %*% C1  #Alteration reads with variants (C=1): theta[2]
    S4 <- t(B1) %*% C1  #Reference reads with variants  (C=1): 1-theta[2]

    ## random initialization for EM
    theta <- c(stats::runif(1, 0.001, 0.25), stats::runif(1, 0.3, 0.6))
    logLik <- 0
    logLik_new <- logLik + logLik_threshold * 2

    ## EM iterations
    for (it in seq_len(max_iter)) {
        # Check convergence
        if ((it > min_iter) && ((logLik_new - logLik) < logLik_threshold))
            break
        logLik <- logLik_new

        #E-step
        logLik_mat <- (S1 * log(theta[1]) + S2 * log(1 - theta[1]) +
                       S3 * log(theta[2]) + S4 * log(1 - theta[2]))
        logLik_mat <- t(t(logLik_mat) + log(Psi))

        logLik_vec <- rep(NA, nrow(logLik_mat))
        for (i in seq_len(nrow(logLik_mat))) {
            logLik_vec[i] <- matrixStats::logSumExp(logLik_mat[i,], na.rm = TRUE)
        }
        logLik_new <- sum(logLik_vec, na.rm = TRUE) + W_log
        logLik_mat_amplify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
        prob_mat <- exp(logLik_mat_amplify) / rowSums(exp(logLik_mat_amplify))

        #M-step
        if (is.na(sum(prob_mat * (S1 + S2))) ||
            sum(prob_mat * (S1 + S2)) == 0)
            theta[1] <- 0.02
        else
            theta[1] <- sum(prob_mat * S1) / sum(prob_mat * (S1 + S2))
        if (is.na(sum(prob_mat * (S3 + S4))) ||
            sum(prob_mat * (S3 + S4)) == 0)
            theta[2] <- 0.75
        else
            theta[2] <- sum(prob_mat * S3) / sum(prob_mat * (S3 + S4))
    }
    if (verbose)
        print(paste("Total iterations:", it))

    ## return values
    return_list <- list("theta" = theta, "prob" = prob_mat, "logLik" = logLik)
    return_list
}


#' Assign cells to clones using a Gibbs sampling algorithm
#'
#' @param A_germ A matrix of integers. Number of alteration reads in germline
#' heterozygous site for variant i cell j
#' @param D_germ A matrix of integers. Number of reads depth in germline
#' heterozygous site for variant i cell j
#' @param prior0 numeric(2), alpha and beta parameters for the Beta prior
#' distribution on the inferred false positive rate.
#' @param prior1 numeric(2), alpha and beta parameters for the Beta prior
#' distribution on the inferred (1 - false negative) rate.
#' @param Bernoulli_threshold An integer. The count threshold of alteration
#' reads when using Bernoulli model.
#' @param wise A string, the wise of parameters for theta1: global, variant,
#' element.
#'
#' @import matrixStats
#' @rdname clone_id
#'
#' @export
#'
#' @examples
#' assignments_binom <- cell_assign_Gibbs(A_clone, D_clone, Config = tree$Z, model = "binomial")
#' assignments_bern <- cell_assign_Gibbs(A_clone, D_clone, Config = tree$Z, model = "Bernoulli")
#'
cell_assign_Gibbs <- function(A, D, Config, Psi=NULL, A_germ=NULL, D_germ=NULL,
                              prior0 = c(0.3, 29.7), prior1 = c(2.25, 2.65),
                              model="binomial", Bernoulli_threshold=1,
                              min_iter=1000, max_iter=10000, wise="variant",
                              verbose=TRUE) {
    if (is.null(Psi)) {Psi <- rep(1/ncol(Config), ncol(Config))}
    if (dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2] ||
        dim(A)[1] != dim(Config)[1] || dim(Config)[2] != length(Psi)) {
        stop(paste0("A and D must have the same size;\n ",
                    "A and Config must have the same number of variants;\n",
                    "Config and Psi must have the same number of clones"))
    }
    if (sum(c("element", "variant", "global") == wise) == 0) {
        stop(paste0("Input wise mode: ", wise,
                    ", while only supporting: element, variant, global"))
    }

    ## preprocessing
    N <- dim(A)[1]             # number of variants
    M <- dim(A)[2]             # number of cells
    K <- dim(Config)[2]             # number of clones
    if (is.null(A_germ))
        A_germ <- matrix(NA, nrow = N, ncol = M)
    if (is.null(D_germ))
        D_germ <- matrix(NA, nrow = N, ncol = M)

    A[which(D == 0)] <- NA
    D[which(D == 0)] <- NA
    A[(D > 0) & is.na(A)] <- 0
    A_germ[which(D_germ == 0)] <- NA
    D_germ[which(D_germ == 0)] <- NA
    A_germ[(D_germ > 0) & is.na(A_germ)] <- 0

    ### test: hack variants phasing
    idx_tmp <- (A > 0) || (A_germ == 0) || (A_germ == D_germ)
    A_germ[idx_tmp] <- D_germ[idx_tmp] <- NA
    idx_tmp <- which(A_germ > (D_germ - A_germ))
    A_germ[idx_tmp] <- D_germ[idx_tmp] - A_germ[idx_tmp]

    C1 <- Config
    C0 <- 1 - Config
    if (model == "Bernoulli") {
        A1 <- A >= Bernoulli_threshold        # bool values
        A2 <- A_germ >= Bernoulli_threshold
        B1 <- 1 - A1
        B2 <- 1 - A2
        W_log <- 0
    }else{
        A1 <- A                  #number of alteration reads
        A2 <- A_germ             #number of alteration reads in germline var
        B1 <- D - A              #number of reference reads
        B2 <- D_germ - A_germ    #number of reference reads in germline var
        W_log <- sum(lchoose(D, A), na.rm = TRUE)  #log binomial coefficients
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
    for (k in seq_len(K)) {
        S1_list[[k]] <- A1 * C0[,k]
        S2_list[[k]] <- B1 * C0[,k]
        S3_list[[k]] <- A1 * C1[,k] + A2   # A2 here: part of likelihood
        S4_list[[k]] <- B1 * C1[,k] + B2
    }

    ## Prepare for sampling
    if (wise == "global") {
        idx_vec <- seq_len(1)         # For: theta1_all[t, ] <- theta1[idx_vec]
        idx_mat <- seq_len(N*M)       # For: theta1[idx_mat] <- theta1_all[t, ]
    }else if (wise == "variant") {
        idx_vec <- seq_len(N)
        idx_mat <- seq_len(N*M)
    }else if (wise == "element") {
        idx_vec <- which(A1 + B1 > 0)
        idx_mat <- which(A1 + B1 > 0)
    }
    n_element <- length(idx_vec)

    if (is.null(dim(prior1)) && length(prior1) == 2) {
        #two variable to a matrix
        prior1 <- t(matrix(rep(prior1, n_element), nrow = 2))
    }
    if (!is.matrix(prior1)) {
        stop("prior1 need to be a matrix of n_element x 2")
    }

    prob_all   <- matrix(0, nrow = max_iter, ncol = M*K)
    logLik_mat <- matrix(0, nrow = M, ncol = K)
    logLik_all <- matrix(0, nrow = max_iter, ncol = 1)
    assign_all <- matrix(0, nrow = max_iter, ncol = M)
    theta0_all <- matrix(0, nrow = max_iter, ncol = 1)
    theta1_all <- matrix(0, nrow = max_iter, ncol = n_element)

    ## Random initialization
    theta0_all[1,1] <- stats::rbeta(1, prior0[1], prior0[2])
    theta1_all[1, ] <- stats::rbeta(rep(1,n_element), prior1[,1], prior1[,2])
    theta1 <- matrix(NA, nrow = N, ncol = M)

    ## Gibbs sampling
    for (it in 2:max_iter) {
        # Update prob_mat
        theta0 <- theta0_all[it - 1, 1]
        theta1[idx_mat] <- theta1_all[it - 1,  ]
        for (k in seq_len(K)) {
            logLik_mat[,k] <- (colSums(S1_list[[k]] * log(theta0), na.rm = TRUE) +
                        colSums(S2_list[[k]] * log(1 - theta0), na.rm = TRUE) +
                        colSums(S3_list[[k]] * log(theta1),     na.rm = TRUE) +
                        colSums(S4_list[[k]] * log(1 - theta1), na.rm = TRUE))
            logLik_mat[,k] <- logLik_mat[,k] + log(Psi[k])
        }
        logLik_mat_amplify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
        prob_mat <- exp(logLik_mat_amplify) / rowSums(exp(logLik_mat_amplify))
        prob_all[it, ] <- prob_mat

        logLik_vec <- rep(NA, nrow(logLik_mat))
        for (i in seq_len(nrow(logLik_mat))) {
            logLik_vec[i] <- matrixStats::logSumExp(logLik_mat[i,], na.rm = TRUE)
        }
        logLik_all[it] <- sum(logLik_vec, na.rm = TRUE) + W_log
        #Update logLikelihood (TODO: add W0 and W1)

        # Sample assignment
        for (j in seq_len(M)) {
            assign_all[it,j] <- sample(seq_len(K), 1, replace = TRUE,
                                       prob = prob_mat[j,])
        }

        # Sample theta with assigned clones
        S1_wgt <- S2_wgt <- 0 # weighted S1
        S3_wgt <- S4_wgt <- matrix(0, nrow = N, ncol = M)
        for (k in seq_len(K)) {
            idx <- which(assign_all[it,] == k)
            S1_wgt <- S1_wgt + sum(S1_list[[k]][,idx], na.rm = TRUE)
            S2_wgt <- S2_wgt + sum(S2_list[[k]][,idx], na.rm = TRUE)
            S3_wgt[,idx] <- S3_wgt[,idx] + S3_list[[k]][,idx]
            S4_wgt[,idx] <- S4_wgt[,idx] + S4_list[[k]][,idx]
        }

        if (wise == "global") {
            S3_wgt[,] <- sum(S3_wgt, na.rm = TRUE)
            S4_wgt[,] <- sum(S4_wgt, na.rm = TRUE)
        }else if (wise == "variant") {
            S3_wgt[,] <- rowSums(S3_wgt, na.rm = TRUE)
            S4_wgt[,] <- rowSums(S4_wgt, na.rm = TRUE)
        }
        theta0_all[it, 1] <- stats::rbeta(1, prior0[1] + S1_wgt,
                                             prior0[2] + S2_wgt)
        theta1_all[it,  ] <- stats::rbeta(rep(1, n_element),
                                          prior1[,1] + S3_wgt[idx_vec],
                                          prior1[,2] + S4_wgt[idx_vec])

        #Check convergence
        if ((it >= min_iter) && (it %% 100 == 0)) {
            FLAG <- Geweke_Z(theta0_all[1:it, 1]) <= 2
            for (n in seq_len(ncol(theta1_all))) {
                if (FLAG == FALSE) {break}
                FLAG <- Geweke_Z(theta1_all[1:it, n]) <= 2
            }
            if (FLAG) {break}
        }
    }
    if (verbose) {print(paste("Converged in", it, "iterations."))}

    ## Return values
    n_buin = ceiling(it * 0.25)

    a <- A1[idx_mat]
    d <- A1[idx_mat] + B1[idx_mat]
    binom_pdf1 <- binom_pdf0 <- rep(0, n_element)
    for (i in seq(n_buin, it)) {
        binom_pdf1 <- binom_pdf1 + stats::dbinom(a, size = d,
                                                 prob = theta1_all[i,])
        binom_pdf0 <- binom_pdf0 + stats::dbinom(a, size = d,
                                                 prob = theta0_all[i])
    }
    prob_variant <- matrix(NA, nrow = N, ncol = M)
    prob_variant[idx_mat] <- binom_pdf1 / (binom_pdf1 + binom_pdf0)
    row.names(prob_variant) <- row.names(A)
    colnames(prob_variant) <- colnames(A)

    # prob_mat <- get_Gibbs_prob(assign_all[1:it,], buin_in=0.25)
    prob_mat <- matrix(colMeans(prob_all[n_buin:it, ]), nrow = M)
    row.names(prob_mat) <- colnames(A)
    colnames(prob_mat) <- colnames(Config)

    theta0 <- mean(theta0_all[n_buin:it,])
    theta1[idx_mat] <- colMeans(as.matrix(theta1_all[n_buin:it,]))

    return_list <- list("theta0" = theta0, "theta1" = theta1,
                        "theta0_all" = as.matrix(theta0_all[1:it,]),
                        "theta1_all" = as.matrix(theta1_all[1:it,]),
                        "element" = idx_mat, "logLik" = logLik_all[1:it],
                        "prob_all" = prob_all[1:it,],
                        "prob" = prob_mat, "prob_variant" = prob_variant)
    return_list
}


#' Geweke diagnostic for MCMC sampling.
#'
#' @param X A matrix of MCMC samples for N samples per K variables
#' @param first A float between 0 and 1. The initial region of MCMC chain.
#' @param last A float between 0 and 1. The final region of MCMC chain.
#'
#' @return \code{Z}, a vector of absolute value of Z scores for each variable.
#' When |Z| <= 2, the sampling could be taken as converged.
#'
Geweke_Z <- function(X, first=0.1, last=0.5) {
    if (is.null(dim(X)))
        X <- as.matrix(X, ncol = 1)
    N <- nrow(X)
    A <- X[1:floor(first*N), , drop = FALSE]
    B <- X[ceiling(last*N):N, , drop = FALSE]

    A_col_mean <- colMeans(A)
    B_col_mean <- colMeans(B)
    A_col_var <- rowSums((t(A) - A_col_mean)^2) / (ncol(A) - 1)
    B_col_var <- rowSums((t(B) - B_col_mean)^2) / (ncol(A) - 1)

    min_var <- 10^(-10)
    Z <- abs(A_col_mean - B_col_mean) / sqrt(A_col_var + B_col_var + min_var)

    Z
}

