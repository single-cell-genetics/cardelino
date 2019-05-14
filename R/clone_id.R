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
#' @param n_proc integer(1), the numebr of processors to use. This parallel 
#' computing can largely reduce time when using multiple chains
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
clone_id <- function(A, D, Config = NULL, n_clone = NULL, Psi = NULL, 
                     relax_Config = FALSE, relax_rate_fixed = NULL,
                     n_chain = 1, model = "binomial", inference = "sampling", 
                     n_proc = 1, verbose = TRUE, ...) {
    model <- match.arg(model, c("binomial", "Bernoulli"))
    inference <- match.arg(inference, c("sampling", "EM"))
    ## check input data
    if (!(all(rownames(A) == rownames(D))))
        stop("Rownames for A and D are not identical.")
    if (!(all(colnames(A) == colnames(D))))
        stop("Colnames for A and D are not identical.")
    if (is.null(Config) && is.null(n_clone))
      stop("Config and n_clone can't be NULL together.")
    
    ## Cardelino-free mode if Config is NULL
    if (is.null(Config)) {
      cat("Config is NULL: de-novo mode is in use.")
      Config <- matrix(0, nrow = nrow(D), ncol = n_clone)
      rownames(Config) <- rownames(D)
      colnames(Config) <- paste0("Clone", seq_len(n_clone))
      relax_Config <-  TRUE
      relax_rate_fixed <- 0.5
    }
    
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
    if (inference == "sampling") {
      #TODO: solve the warnings for using foreach
      library(foreach)
      doMC::registerDoMC(n_proc)
      
      ids_list <- foreach::foreach(ii = 1:n_chain) %dopar% {
        cell_assign_Gibbs(A, D, Config, Psi = Psi, 
                          relax_Config = relax_Config, 
                          relax_rate_fixed = relax_rate_fixed,
                          model = model, verbose = verbose, ...)
      }
      
      ids_out <- ids_list[[1]]
      ids_out$n_chain <- 1
      if (n_chain > 1) {
        for (ii in seq(2, n_chain)) {
          ids_out$n_chain <- ids_out$n_chain + 1
          idx <- colMatch(ids_out$prob, ids_list[[ii]]$prob, force = TRUE)
          ids_out$prob <- ids_out$prob + ids_list[[ii]]$prob[, idx]
          ids_out$relax_rate <- ids_out$relax_rate + ids_list[[ii]]$relax_rate
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
        return(cell_assign_EM(A, D, Config, model = model, verbose = verbose,
                              ...))
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
#' @param relax_Config logical(1), If TRUE, relaxing the Clone Configuration by
#' changing it from fixed value to act as a prior Config with a relax rate.
#' @param relax_rate_fixed numeric(1), If the value is between 0 to 1, 
#' the relax rate will be set as a fix value during updating clone Config. If
#' NULL, the relax rate will be learned automatically with relax_rate_prior.
#' @param relax_rate_prior numeric(2), the two parameters of beta prior 
#' distribution of the relax rate for relaxing the clone Configuration. This 
#' mode is used when relax_relax is NULL.
#' @param keep_base_clone bool(1), if TRUE, keep the base clone of Config to its
#' input values when relax mode is used.
#' @param Bernoulli_threshold An integer. The count threshold of alteration
#' reads when using Bernoulli model.
#' @param wise A string, the wise of parameters for theta1: global, variant,
#' element.
#' @param relabel bool(1), if TRUE, relabel the samples of both Config and prob 
#' during the Gibbs sampleing.
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
                              relax_Config=FALSE, relax_rate_fixed=NULL, 
                              relax_rate_prior=c(1, 9), keep_base_clone=TRUE,
                              prior0=c(0.2, 99.8), prior1=c(0.45, 0.55),
                              model="binomial", Bernoulli_threshold=1,
                              min_iter=3000, max_iter=10000, wise="variant",
                              relabel=FALSE, verbose=TRUE) {
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
    K <- dim(Config)[2]        # number of clones
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

    ### test: hack variants phasing (TODO)
    # idx_tmp <- (A > 0) || (A_germ == 0) || (A_germ == D_germ)
    # A_germ[idx_tmp] <- D_germ[idx_tmp] <- NA
    # idx_tmp <- which(A_germ > (D_germ - A_germ))
    # A_germ[idx_tmp] <- D_germ[idx_tmp] - A_germ[idx_tmp]

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
    Config_all <- matrix(0, nrow = max_iter, ncol = N*K)
    relax_rate_all <- matrix(0, nrow = max_iter, ncol = 1)

    if (!is.null(relax_Config) && relax_Config != FALSE) {
        if (!is.null(relax_rate_fixed)) {
          if (relax_rate_fixed > 1 || relax_rate_fixed < 0) {
            stop("Error: relax_rate_fixed needs to be NULL or in [0, 1].")
          }
          relax_rate <- relax_rate_fixed ## fixed relax_rate
          relax_rate_all[,] <- relax_rate
        } else if (!is.null(relax_rate_prior)) {
          relax_rate <- relax_rate_prior[1] / (relax_rate_prior[1] + 
                                               relax_rate_prior[2])
        } else {
          stop("Error: require value for either relax_Config or relax_prior.")
        }

        Config_new <- Config
        Config_prior <- Config
        Config_prior[Config == 1] <- 1 - relax_rate
        Config_prior[Config == 0] <- relax_rate
        if (keep_base_clone) {
            Config_prior[, 1] <- Config[, 1]}
        Config_prior_oddlog <- log(Config_prior) - log(1 - Config_prior)
        Iden_mat <- matrix(0, nrow = M, ncol = K)
    }

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

        ## Update Config
        if (it > (0.1 * min_iter) && !is.null(relax_Config) &&
                  relax_Config != FALSE) {
            if (it > (0.1 * min_iter + 5) && is.null(relax_rate_fixed)) {
                diff0 <- sum((Config == Config_new)[, 2:ncol(Config)])
                diff1 <- sum((Config != Config_new)[, 2:ncol(Config)])
                relax_rate <- stats::rbeta(1, relax_rate_prior[1] + diff1,
                                              relax_rate_prior[2] + diff0)
                relax_rate_all[it] <- relax_rate

                Config_prior <- Config
                Config_prior[Config == 1] <- 1 - relax_rate
                Config_prior[Config == 0] <- relax_rate
                if (keep_base_clone) {
                    Config_prior[, 1] <- Config[, 1]}
                Config_prior_oddlog <- log(Config_prior) - log(1 - Config_prior)
            }

            Iden_mat[,] <- 0
            for (j in seq_len(M)) {
                Iden_mat[j, assign_all[it,j]] <- 1 }

            # calculate log_probability matrix with genotype 0 and 1
            P0_mat <- A1 * log(theta0) + B1 * log(1 - theta0) + W_log
            P1_mat <- A1 * log(theta1) + B1 * log(1 - theta1) + W_log

            oddR_log <- P1_mat %*% Iden_mat - P0_mat %*% Iden_mat
            oddR_log <- oddR_log + Config_prior_oddlog
            oddR_log[which(oddR_log > 50)] <- 50
            oddR_log[which(oddR_log < -50)] <- -50
            Conf_prob <- exp(oddR_log) / (exp(oddR_log) + 1)

            Config_new[,] <- stats::rbinom(N*K, size = 1, Conf_prob)
            Config_all[it, ] <- Config_new

            for (k in seq_len(K)) {
                S1_list[[k]] <- A1 * (1 - Config_new[,k])
                S2_list[[k]] <- B1 * (1 - Config_new[,k])
                S3_list[[k]] <- A1 * Config_new[,k] + A2
                S4_list[[k]] <- B1 * Config_new[,k] + B2
            }
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
            Converged_all <- rep(FALSE, ncol(prob_all))
            for (n in seq_len(ncol(prob_all))) {
                Converged_all[n] <- Geweke_Z(prob_all[1:it, n]) <= 2
            }
            cat(paste0(round(mean(Converged_all, na.rm = TRUE), 3) * 100, 
                       "% converged."))
            if (mean(Converged_all, na.rm = TRUE) > 0.999) {break}
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
    if (relabel) {
      for (ii in seq(n_buin, it)) {
        prob1 <- matrix(prob_all[ii - 1, ], nrow = M)
        prob2 <- matrix(prob_all[ii - 1, ], nrow = M)
        idx <- colMatch(prob1, prob2, force = TRUE)
        prob_all[ii, ] <- prob2[idx]
        Config_all[ii, ] <- matrix(Config_all[ii, ], nrow = N)[, idx]
      }
    }
    prob_mat <- matrix(colMeans(prob_all[n_buin:it, ]), nrow = M)
    row.names(prob_mat) <- colnames(A)
    colnames(prob_mat) <- colnames(Config)

    Config_prob <- Config
    Config_prob[,] <- colMeans(Config_all[n_buin:it,])

    theta0 <- mean(theta0_all[n_buin:it,])
    theta1[idx_mat] <- colMeans(as.matrix(theta1_all[n_buin:it,]))

    return_list <- list("theta0" = theta0, "theta1" = theta1,
                        "theta0_all" = as.matrix(theta0_all[1:it,]),
                        "theta1_all" = as.matrix(theta1_all[1:it,]),
                        "element" = idx_mat, "logLik" = logLik_all[1:it],
                        "prob_all" = prob_all[1:it,],
                        "prob" = prob_mat, "prob_variant" = prob_variant,
                        "relax_rate" = mean(relax_rate_all[n_buin:it]),
                        "Config_prob" = Config_prob,
                        "Config_all" = Config_all[1:it,],
                        "relax_rate_all" = relax_rate_all[1:it])
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


#' Get a clonal tree from a configuration matrix
#'
#' @param Config variant x clone matrix of binary values. The clone-variant
#' configuration, which encodes the phylogenetic tree structure. This is the
#' output Z of Canopy
#' @param P a one-column numeric matrix encoding the (observed or estimated)
#' prevalence (or frequency) of each clone
#' @param strictness character(1), a character string defining the strictness of
#' the function if there are all-zero rows in the Config matrix. If \code{"lax"}
#' then the function silently drops all-zero rows and proceeds. If \code{"warn"}
#' then the function warns of dropping all-zero rows and proceeds. If
#' \code{"error"} then the function throws an error is all-zero rows are
#' detected.
#'
#' @return
#' An object of class "phylo" describing the tree structure. The output object
#' also contains an element "sna" defining the clustering of variants onto the
#' branches of the tree, and if \code{P} is non-null it also contains VAF
#' (variant allele frequency), CCF (cell clone fraction) and clone prevalence
#' values (computed from the supplied \code{P} argument).
#'
#' @details
#' Output tree may be nonsensical if the input \code{Config} matrix does not
#' define a coherent tree structure.
#'
#' @author Davis McCarthy
#'
#' @import utils
#' @export
#'
#' @examples
#' Configk3 <- matrix(c(rep(0, 15), rep(1, 8), rep(0, 7), rep(1, 5), rep(0, 3), rep(1, 7)), ncol = 3)
#' tree_k3 <- get_tree(Config = Configk3, P = matrix(rep(1/3, 3), ncol = 1))
#' plot_tree(tree_k3)
get_tree <- function(Config, P = NULL, strictness = "lax") {
  if (!is.null(P)) {
    if (ncol(P) != 1)
      stop("P must be a matrix with one column encoding clone prevalence values")
  }
  all_zero_rows <- rowSums(Config) == 0
  strictness <- match.arg(strictness, c("lax", "warn", "error"))
  if (any(all_zero_rows)) {
    if (strictness == "error")
      stop("Config matrix contains all-zero rows.")
    else {
      if (strictness == "warn")
        warning(paste("Dropped", sum(all_zero_rows), "all-zero rows from Config matrix."))
      else
        message(paste("Dropped", sum(all_zero_rows), "all-zero rows from Config matrix."))
      Config <- Config[!all_zero_rows,]
    }
  }
  k <- ncol(Config) # number of clones
  varnames <- rownames(Config)
  sna <- matrix(nrow = nrow(Config), ncol = 3)
  sna[, 1] <- seq_len(nrow(sna))
  rownames(sna) <- varnames
  colnames(sna) <- c("sna", "sna.st.node", "sna.ed.node")
  tip_label <- seq_len(k)
  tip_vals <- 2^seq_len(k)
  ## Need to determine number of internal nodes in the tree
  Config2 <- t(t(Config) * tip_vals)
  var_bin_vals <- rowSums(Config2)
  node_vals <- names(table(var_bin_vals))
  node_vals <- node_vals[!(node_vals %in% tip_vals)]
  node_num <- sum(!(node_vals %in% tip_vals))
  ## define a list with subsets of edge matrices
  ## start with root node (k+1), which always connects to tip 1 (base clone)
  if (node_num > 0.5) {
    node_def_list <- list()
    edge_list <- list()
    for (i in 1:(k - 1)) {
      clone_combos <- utils::combn(2:k, (k - i), simplify = FALSE)
      for (j in seq_len(length(clone_combos))) {
        test_sum <- sum(tip_vals[clone_combos[[j]]])
        if (test_sum %in% node_vals)
          node_def_list[[
            paste0("node", paste0(clone_combos[[j]]), collapse = "_")]] <-
            clone_combos[[j]]
      }
    }
    if (node_num != length(node_def_list))
      stop("Conflict in computed number of internal nodes.")
    ## Sort out edges for the root node
    tip_nodes <- seq_len(k)
    root_to_tip <- tip_nodes[!(tip_nodes %in% unique(unlist(node_def_list)))]
    edge_list[["root_node_tips"]] <- matrix(
      c(rep(k + 1, length(root_to_tip)), root_to_tip),
      nrow = length(root_to_tip))
    el_counter <- 1
    for (i in seq_len(length(node_def_list))) {
      ## add edge from root to internal node if not already done
      if (i < 1.5) {
        el_counter <- el_counter + 1
        edge_list[[el_counter]] <- matrix(c(k + 1, k + 1 + i), nrow = 1)
        sna[var_bin_vals == sum(2^node_def_list[[i]]), 2] <- k + 1
        sna[var_bin_vals == sum(2^node_def_list[[i]]), 3] <- k + 1 + i
      } else {
        clones_in_this_node <- node_def_list[[i]]
        clones_in_prev_nodes <- unique(unlist(node_def_list[1:(i - 1)]))
        if (!any(clones_in_this_node %in% clones_in_prev_nodes)) {
          el_counter <- el_counter + 1
          edge_list[[el_counter]] <- matrix(c(k + 1, k + 1 + i), nrow = 1)
          sna[var_bin_vals == sum(2^node_def_list[[i]]), 2] <- k + 1
          sna[var_bin_vals == sum(2^node_def_list[[i]]), 3] <- k + 1 + i
        }
      }
      ## add edge from internal node to internal node
      ## if all of the clones for the node are present in the previous node in
      ## the tree (immediately above in the hierarchy), then add the edge
      ## check the size of previous nodes, and select the node that has minimum
      ## number of clones that is more than the number in this node
      if (i > 1.5) {
        prev_nodes <- 1:(i - 1)
        prev_node_sizes <- sapply(node_def_list[prev_nodes], length)
        prev_nodes <- prev_nodes[prev_node_sizes > length(node_def_list[[i]])]
        min_prev_node_size <- min(prev_node_sizes[prev_nodes])
        prev_nodes <- prev_nodes[prev_node_sizes[prev_nodes] ==
                                   min_prev_node_size]
        for (j in prev_nodes) {
          if (all(node_def_list[[i]] %in% node_def_list[[j]])) {
            el_counter <- el_counter + 1
            edge_list[[el_counter]] <- matrix(
              c(k + 1 + j, k + 1 + i), nrow = 1)
            sna[var_bin_vals == sum(2^node_def_list[[i]]), 2] <- k + 1 + j
            sna[var_bin_vals == sum(2^node_def_list[[i]]), 3] <- k + 1 + i
          }
        }
      }
      ## add edge from internal node to tip
      ## (if clone not present in any subsequent nodes)
      if (node_num < 1.5) {
        ## if only one internal node, there are edges from this node to all tips
        node_to_tip <- tip_nodes[-1]
        el_counter <- el_counter + 1
        edge_list[[el_counter]] <- matrix(
          c(rep(k + 1 + i, length(node_to_tip)), node_to_tip),
          nrow = length(node_to_tip))
        for (m in node_to_tip) {
          sna[var_bin_vals == sum(2^m), 2] <- k + 1 + i
          sna[var_bin_vals == sum(2^m), 3] <- m
        }
      } else {
        ## if more than one internal node, need to check if tips mentioned in
        ## this node appear in any subsequent nodes
        node_to_tip <- node_def_list[[i]]
        if (i < node_num) {
          node_to_tip <- node_to_tip[
            !(node_to_tip %in% unique(unlist(node_def_list[(i + 1):node_num])))]
        } ## else this is the last node; just connect edges from node to tips
        if (length(node_to_tip) > 0.5) {
          el_counter <- el_counter + 1
          edge_list[[el_counter]] <- matrix(
            c(rep(k + 1 + i, length(node_to_tip)), node_to_tip),
            nrow = length(node_to_tip))
          for (m in node_to_tip) {
            sna[var_bin_vals == sum(2^m), 2] <- k + 1 + i
            sna[var_bin_vals == sum(2^m), 3] <- m
          }
        }
      }
    }
  } else {
    edge_list <- list("root_node" = matrix(c(rep(k + 1, k), seq_len(k)),
                                           ncol = 2))
    for (j in 2:k) {
      sna[var_bin_vals == 2^j, 2] <- k + 1
      sna[var_bin_vals == 2^j, 3] <- j
    }
  }
  # node_def_list
  edge_mat <- do.call(rbind, edge_list)
  tree_out <- list(edge = edge_mat, Nnode = node_num + 1, tip.label = tip_label)
  class(tree_out) <- "phylo"
  tree_out$Z <- Config
  if (!is.null(P)) {
    tree_out$P <- P
    tree_out$VAF <- tree_out$Z %*% tree_out$P / 2
    tree_out$CCF <- tree_out$Z %*% tree_out$P
  }
  tree_out$sna <- sna
  tree_out
}


### Create a tree from a config matrix
# Configk3 <- matrix(c(rep(0, 15), rep(1, 8), rep(0, 7), rep(1, 5), rep(0, 3),
#                      rep(1, 7)), ncol = 3)
# Configk4_1 <- matrix(c(rep(0, 15), rep(1, 8), rep(0, 7), rep(1, 5), rep(0, 3),
#                        rep(1, 4), rep(0, 3), rep(1, 5), rep(0, 7), rep(1, 3)),
#                      ncol = 4)
# Configk4_2 <- matrix(c(rep(0, 15), rep(1, 8), rep(0, 7), rep(1, 5), rep(0, 3),
#                        rep(1, 4), rep(0, 3), rep(1, 5), rep(0, 3), rep(1, 2),
#                        rep(0, 2), rep(1, 3)),
#                      ncol = 4)
# Configk4_2_bad <- matrix(c(rep(0, 15), rep(1, 8), rep(0, 2), rep(1, 2),
#                            rep(0, 3), rep(1, 5), rep(0, 3),
#                            rep(1, 4), rep(0, 3), rep(1, 5), rep(0, 3), rep(1, 2),
#                            rep(0, 2), rep(1, 3)),
#                          ncol = 4)
#
# Configk5 <- matrix(c(rep(0, 15), rep(1, 8), rep(0, 7), rep(1, 5), rep(0, 3),
#                      rep(1, 4), rep(0, 3), rep(1, 5), rep(0, 3), rep(1, 2),
#                      rep(0, 2), rep(0, 3), rep(1, 8), rep(0, 4), rep(1, 3)),
#                    ncol = 5)
#
# tree_k3 <- get_tree(Config = Configk3, P = matrix(rep(1/3, 3), ncol = 1))
# tree_k3$sna
# tree_k3$Z
# tree_k3$edge
# cardelino::plot_tree(tree_k3)
#
# tree_k4_1 <- get_tree(Config = Configk4_1, P = matrix(rep(1/4, 4), ncol = 1))
# tree_k4_1$sna
# tree_k4_1$Z
# tree_k4_1$edge
# plot_tree(tree_k4_1)
#
# tree_k4_2 <- get_tree(Config = Configk4_2, P = matrix(rep(1/4, 4), ncol = 1))
# tree_k4_2$sna
# tree_k4_2$Z
# tree_k4_2$edge
# plot_tree(tree_k4_2)
#
# tree_k4_2_bad <- get_tree(Config = Configk4_2_bad, P = matrix(rep(1/4, 4), ncol = 1))
# tree_k4_2_bad$sna
# tree_k4_2_bad$Z
# tree_k4_2_bad$edge
# plot_tree(tree_k4_2_bad)
#
# tree_k5 <- get_tree(Config = Configk5, P = matrix(rep(1/5, 5), ncol = 1))
# tree_k5$sna
# tree_k5$Z
# tree_k5$edge
# plot_tree(tree_k5)
#
# ## joxm - works
# card_joxm <- readRDS("../fibroblast-clonality/data/cell_assignment/cardelino_results_carderelax.joxm.filt_lenient.cell_coverage_sites.rds")
# names(card_joxm)
# card_joxm_Config_prob <- card_joxm$tree$Z
# card_joxm_Config_prob[,] <- colMeans(card_joxm$Config_all)
# card_joxm_Config_best <- round(card_joxm_Config_prob)
#
# tree_joxm <- get_tree(card_joxm_Config_best,
#                       P = matrix(colMeans(card_joxm$prob_mat > 0.5), ncol = 1))
# p1 <- cardelino::plot_tree(tree_joxm, orient = "v") +
#   ggtitle("joxm: Cardelino output tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# p2 <- cardelino::plot_tree(card_joxm$tree, orient = "v") +
#   ggtitle("joxm: Canopy tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# cowplot::plot_grid(p2, p1, ncol = 1)
#
# ## zoxy - fails
# card_zoxy <- readRDS("../fibroblast-clonality/data/cell_assignment/cardelino_results_carderelax.zoxy.filt_lenient.cell_coverage_sites.rds")
# names(card_zoxy)
# card_zoxy_Config_prob <- card_zoxy$tree$Z
# card_zoxy_Config_prob[,] <- colMeans(card_zoxy$Config_all)
# card_zoxy_Config_best <- round(card_zoxy_Config_prob)
#
# tree_zoxy <- get_tree(card_zoxy_Config_best,
#                       P = matrix(colMeans(card_zoxy$prob_mat > 0.5), ncol = 1))
# tree_zoxy$edge
# tree_zoxy$sna
# tree_zoxy <- get_tree(card_zoxy_Config_best,
#                       P = matrix(colMeans(card_zoxy$prob_mat > 0.5), ncol = 1),
#                       strictness = "warn")
# tree_zoxy <- get_tree(card_zoxy_Config_best,
#                       P = matrix(colMeans(card_zoxy$prob_mat > 0.5), ncol = 1),
#                       strictness = "e")
# ## handles all-zero rows
# ## looks like a problem caused by no shared variants between cl2 and cl3
# p1 <- cardelino::plot_tree(tree_zoxy, orient = "v") +
#   ggtitle("zoxy: Cardelino output tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# p2 <- cardelino::plot_tree(card_zoxy$tree, orient = "v") +
#   ggtitle("zoxy: Canopy tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# cowplot::plot_grid(p2, p1, ncol = 1)
#
# ## lexy - some all-zero rows
# card_lexy <- readRDS("../fibroblast-clonality/data/cell_assignment/cardelino_results_carderelax.lexy.filt_lenient.cell_coverage_sites.rds")
# names(card_lexy)
# card_lexy_Config_prob <- card_lexy$tree$Z
# card_lexy_Config_prob[,] <- colMeans(card_lexy$Config_all)
# card_lexy_Config_best <- round(card_lexy_Config_prob)
#
# tree_lexy <- get_tree(card_lexy_Config_best,
#                       P = matrix(colMeans(card_lexy$prob_mat > 0.5), ncol = 1))
# tree_lexy$edge
# ##
# p1 <- cardelino::plot_tree(tree_lexy, orient = "v") +
#   ggtitle("lexy: Cardelino output tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# p2 <- cardelino::plot_tree(card_lexy$tree, orient = "v") +
#   ggtitle("lexy: Canopy tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# cowplot::plot_grid(p2, p1, ncol = 1)
#
# ## vils - works
# card_vils <- readRDS("../fibroblast-clonality/data/cell_assignment/cardelino_results_carderelax.vils.filt_lenient.cell_coverage_sites.rds")
# names(card_vils)
# card_vils_Config_prob <- card_vils$tree$Z
# card_vils_Config_prob[,] <- colMeans(card_vils$Config_all)
# card_vils_Config_best <- round(card_vils_Config_prob)
#
# tree_vils <- get_tree(
#   card_vils_Config_best,
#   P = matrix(colSums(card_vils$prob_mat > 0.5) /
#                sum(rowMax(card_vils$prob_mat) > 0.5), ncol = 1))
# tree_vils$edge
# tree_vils$P
# ##
# p1 <- cardelino::plot_tree(tree_vils, orient = "v") +
#   ggtitle("vils: Cardelino output tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# p2 <- cardelino::plot_tree(card_vils$tree, orient = "v") +
#   ggtitle("vils: Canopy tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# cowplot::plot_grid(p2, p1, ncol = 1)
#
# ## sehl - fails
# card_sehl <- readRDS("../fibroblast-clonality/data/cell_assignment/cardelino_results_carderelax.sehl.filt_lenient.cell_coverage_sites.rds")
# names(card_sehl)
# card_sehl_Config_prob <- card_sehl$tree$Z
# card_sehl_Config_prob[,] <- colMeans(card_sehl$Config_all)
# card_sehl_Config_best <- (card_sehl_Config_prob > 0.49) * 1L
#
# tree_sehl <- get_tree(
#   card_sehl_Config_best,
#   P = matrix(colSums(card_sehl$prob_mat > 0.5) /
#                sum(rowMax(card_sehl$prob_mat) > 0.5), ncol = 1))
# tree_sehl$edge
# tree_sehl$P
# ## Too many edges - looks like it finds an internal node that does not make sense
# p1 <- cardelino::plot_tree(tree_sehl, orient = "v") +
#   ggtitle("sehl: Cardelino output tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# p2 <- cardelino::plot_tree(card_sehl$tree, orient = "v") +
#   ggtitle("sehl: Canopy tree") +
#   theme(plot.title = element_text(hjust = 0.5))
# cowplot::plot_grid(p2, p1, ncol = 1)


#
#
# ## example data - fails; finds too many internal nodes
# data(example_donor)
# plot_tree(get_tree(tree$Z))
# assignments_binom <- cell_assign_Gibbs(A_clone, D_clone, Config = tree$Z,
#                                        relax_Config = TRUE, model = "binomial")
# eg_Config_prob <- assignments_binom$Config_prob
# eg_Config_best <- round(eg_Config_prob)
#
# tree_eg <- get_tree(eg_Config_best,
#                       P = matrix(colMeans(assignments_binom$prob > 0.5), ncol = 1))
# tree_eg$edge
# ## too many internal nodes
# ggtree::ggtree(tree_eg)
# plot_tree(tree_eg)
# ## gives a new and different error *shrug*
