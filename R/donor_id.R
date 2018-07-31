## Donor deconvolution in multiplexed scRNA-seq.

#' Donor deconvolution of scRNA-seq data
#'
#' @param cell_vcf_file character(1), path to a VCF file containing variant data
#' for cells
#' @param donor_vcf_file character(1), path to a VCF file containing reference
#' donor genotypes
#' @param check_doublet logical(1), should the function check for doublet cells?
#' @param n_vars_threshold integer(1), if the number of variants with coverage
#' in a cell is below this threshold, then the cell will be given an
#' "unassigned" donor ID (default: 10)
#' @param p_threshold numeric(1), threshold for posterior probability of
#' donor assignment (must be in [0, 1]); if best posterior probability for a
#' donor is greater then the threshold, then the cell is assigned to that donor
#' (as long as the cell is not determined to be a doublet) and if below the
#' threshold, then the cell's donor ID is "unassigned"
#' @param verbose logical(1), should the function output verbose information
#' while running?
#' @param genome character(1), string indicating the genome build used in the
#' VCF file(s) (default: "GRCh37")
#' @param seq_levels_style character(1), string passed to
#' \code{\link[GenomeInfoDb]{seqlevelsStyle}} the style to use for
#' chromosome/contig names (default: "Ensembl")
#' @param donors optional character vector providing a set of donors to use, by
#' subsetting the donors present in the \code{donor_vcf_file}; if \code{NULL}
#' (default) then all donors present in VCF will be used.
#' @param ... arguments passed to \code{donor_id_GT}
#'
#' @details This function reads in all elements of the provided VCF file(s) into
#' memory, so we highly recommend filtering VCFs to the minimal appropriate set
#' of variants (e.g. with the bcftools software) before applying them to this
#' function.
#'
#' @author Davis McCarthy and Yuanhua Huang
#'
#' @export
#'
#' @examples
#' ids <- donor_id(system.file("extdata", "cells.donorid.vcf.gz", package = "cardelino"),
#'                 system.file("extdata", "donors.donorid.vcf.gz", package = "cardelino"))
#' table(ids$assigned$donor_id)
#'
donor_id <- function(cell_vcf_file, donor_vcf_file = NULL, check_doublet = TRUE,
                     n_vars_threshold = 10, p_threshold = 0.95,
                     verbose = TRUE, genome = "GRCh37",
                     seq_levels_style = "Ensembl", donors = NULL, ...) {

    if (is.null(donor_vcf_file))
        stop("Donor deconvolution without reference genotypes not yet implemented.")
    else {
        vcf_cell <- read_vcf(cell_vcf_file, genome = genome,
                             seq_levels_style = seq_levels_style,
                             verbose = verbose)
        vcf_donor <- read_vcf(donor_vcf_file, genome = genome,
                             seq_levels_style = seq_levels_style,
                             verbose = verbose)
        in_data <- get_snp_matrices(vcf_cell, vcf_donor, verbose = verbose,
                                    donors = donors)
        if (verbose)
            message("Donor ID using ", nrow(in_data$A), " variants")
        out <- donor_id_GT(in_data$A, in_data$D, in_data$GT_donors,
                           check_doublet = check_doublet,
                           n_vars_threshold = n_vars_threshold,
                           p_threshold = p_threshold,
                           verbose = verbose, ...)
        out$A <- in_data$A
        out$D <- in_data$D
        out$GT <- in_data$GT_cells
    }
    out
}

## test example
ids <- donor_id("inst/extdata/cells.donorid.vcf.gz",
         "inst/extdata/donors.donorid.vcf.gz")
table(ids$assigned$donor_id)
head(ids$assigned)
ggplot(ids$assigned, aes(n_vars, prob_doublet, colour = donor_id)) +
    geom_point(size = 3, alpha = 0.5) +
    #ggthemes::scale_color_canva(palette = "Surf and turf")
    theme_bw()
dplyr::filter(ids$assigned, donor_id == "doublet")

sce <- readRDS("../hipsci-fibro/Data/processed/sce_merged_donors_cardelino_donorid_all_with_qc_labels.rds")
coldat <- colData(sce)
rm(sce)
ids$assigned$sample_id <- ids$assigned$cell
iddf <- dplyr::inner_join(ids$assigned, as.data.frame(coldat))
iddf %>% dplyr::filter(donor_id == "doublet") %>%
    dplyr::select(sample_id, well_condition, well_type, prob_max, prob_doublet) %>% as.data.frame


ids2_sing <- donor_id("inst/extdata/cells.donorid.vcf.gz",
                      "inst/extdata/donors.allfibro.vcf.gz",
                      check_doublet = FALSE)
table(ids2_sing$assigned$donor_id)
dplyr::filter(ids2_sing$assigned, grepl("laia", donor_id))
ids2_doub <- donor_id("inst/extdata/cells.donorid.vcf.gz",
                      "inst/extdata/donors.allfibro.vcf.gz",
                      check_doublet = TRUE,
                      donors = unique(ids2_sing$assigned$donor_id))
table(ids2_doub$assigned$donor_id)
ggplot(ids2_doub$assigned, aes(n_vars, prob_doublet, colour = donor_id)) +
    geom_point(size = 3, alpha = 0.5) +
    #ggthemes::scale_color_canva(palette = "Surf and turf")
    theme_bw()
ids2_doub$assigned$sample_id <- ids2_doub$assigned$cell
id2df <- dplyr::inner_join(ids2_doub$assigned, as.data.frame(coldat))
id2df %>% dplyr::filter(donor_id %in% c("doublet", "unassigned")) %>%
    dplyr::select(cell, donor_id, well_condition, well_type, prob_max, prob_doublet) %>% as.data.frame

dplyr::filter(ids2_doub$assigned, donor_id %in% c("doublet", "unassigned")) %>% as.data.frame
dplyr::filter(ids2_doub$assigned, grepl("laia", donor_id))

## ppca on cell genotypes
pp <- pcaMethods::ppca(t(ids2_doub$GT))
id2df$PPCA1 <- pp@scores[, 1]
id2df$PPCA2 <- pp@scores[, 2]
ggplot(id2df, aes(PPCA1, PPCA2, colour = donor_id)) +
    geom_point(size = 3, alpha = 0.5) +
    #ggthemes::scale_color_canva(palette = "Surf and turf")
    theme_bw()


#' EM algorithm for donor deconvolution in multiplexed scRNA-seq with genotypes.
#'
#' @param A A matrix of integers. Number of alteration reads in SNP i cell j
#' @param D A matrix of integers. Number of reads depth in SNP i cell j
#' @param GT A matrix of integers for genotypes. The donor-SNP configuration,
#' the element should be 0, 1, 2, but others are probably also appliable.
#' @param check_doublet logical(1), if TRUE, check doublet, otherwise ignore.
#' @param n_vars_threshold integer(1), if the number of variants with coverage
#' in a cell is below this threshold, then the cell will be given an
#' "unassigned" donor ID (default: 10)
#' @param p_threshold numeric(1), threshold for posterior probability of
#' donor assignment (must be in [0, 1]); if best posterior probability for a
#' donor is greater then the threshold, then the cell is assigned to that donor
#' (as long as the cell is not determined to be a doublet) and if below the
#' threshold, then the cell's donor ID is "unassigned"
#' @param Psi A vector of float. The fractions of each donor; default uniform.
#' @param min_iter A integer. The minimum number of iterations in EM algorithm.
#' @param max_iter A integer. The maximum number of iterations in EM algorithm.
#' The real iteration may finish earlier.
#' @param logLik_threshold A float. The threshold of logLikelihood increase for
#' detecting convergence.
#' @param verbose logical(1), If TRUE, output verbose information when running.
#'
#' @details Users should typically use \code{\link{donor_id}} rather than this
#' lower-level function.
#'
#' @return a list containing
#' \code{logLik}, the log likelihood.
#' \code{theta}, a vector denoting the binomial parameters for each genotype.
#' \code{prob}, a matrix of posterior probability of cell assignment to donors.
#' The summary may less than 1, as there are some probabilities go to doublets.
#' \code{prob_doublet}, a matrix of posterior probability of cell assignment to
#' each inter-donor doublet.
#'
#' @import stats
#'
#' @export
#'
#' @examples
#' data(example_donor)
#' ids <- donor_id_GT(A, D, GT = tree$Z[1:nrow(A),])
#' table(ids$assigned$donor_id)
#'
donor_id_GT <- function(A, D, GT, check_doublet=TRUE, n_vars_threshold = 10,
                        p_threshold = 0.95,
                        Psi=NULL, min_iter=10, max_iter=1000,
                        logLik_threshold=1e-5, verbose=TRUE) {
    ## Create doublet genotypes
    K_doublet <- 0
    K_singlet <- ncol(GT)
    if (check_doublet) {
        donor_ids <- colnames(GT)
        if (is.null(donor_ids)) { donor_ids <- paste0("donor", seq_len(ncol(GT))) }

        combn_idx <- utils::combn(ncol(GT), 2)
        K_doublet <- ncol(combn_idx)
        GT_doublet <- (GT[,combn_idx[1,]] + GT[,combn_idx[2,]]) / 2
        colnames(GT_doublet) <- paste0(donor_ids[combn_idx[1,]], ",",
                                       donor_ids[combn_idx[2,]])

        GT <- cbind(GT, GT_doublet)
    }
    ## Check input data
    if (is.null(Psi)) {
        Psi <- rep(1/ncol(GT), ncol(GT))
    }
    if (dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2]) {
        stop("A and D must have the same size.")
    }
    if (dim(A)[1] != dim(GT)[1]) {
        stop("nrow(A) and nrow(GT) must be the same.")
    }
    if (dim(GT)[2] != length(Psi)) {
        stop("ncol(GT) and length(Psi) must be the same.")
    }

    ## preprocessing
    N <- dim(A)[1]        # number of variants
    M <- dim(A)[2]        # number of cells
    K <- dim(GT)[2]       # number of clones

    n_vars <- matrixStats::colSums2(!is.na(D))
    A[is.na(A)] <- 0
    D[is.na(D)] <- 0
    W_log <- sum(lchoose(D, A), na.rm = TRUE)  #log binomial coefficients

    GT_uniq <- sort(unique(c(GT)))
    S1 <- S2 <- SS <- list()
    for (ig in seq_len(length(GT_uniq))) {
        #print(c(ig, GT_uniq[ig]))
        S1[[ig]] <- t(A) %*% (GT == GT_uniq[ig])
        SS[[ig]] <- t(D) %*% (GT == GT_uniq[ig])
        S2[[ig]] <- SS[[ig]] - S1[[ig]]
    }

    ## random initialization for EM
    theta <- matrix(0.98 * GT_uniq / max(GT_uniq) + 0.01, ncol = 1)
    row.names(theta) <- paste0("GT=", GT_uniq)
    logLik <- 0
    logLik_new <- logLik + logLik_threshold * 2

    ## EM iterations
    for (it in seq_len(max_iter)) {
        # Check convergence
        if ((it > min_iter) && ((logLik_new - logLik) < logLik_threshold))
            break
        logLik <- logLik_new

        # E-step: update assignment posterior probability
        logLik_mat <- matrix(0, nrow = M, ncol = K)
        for (ig in seq_len(length(GT_uniq))) {
            logLik_mat <- (logLik_mat + S1[[ig]] * log(theta[ig]) +
                               S2[[ig]] * log(1 - theta[ig]))
        }

        logLik_vec <- rep(NA, nrow(logLik_mat))
        for (i in seq_len(nrow(logLik_mat))) {
            logLik_vec[i] <- matrixStats::logSumExp(logLik_mat[i,], na.rm = TRUE)
        }
        logLik_new <- sum(logLik_vec, na.rm = TRUE) + W_log
        logLik_mat_amplify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
        prob_mat <- exp(logLik_mat_amplify) / rowSums(exp(logLik_mat_amplify))

        # M-step: maximumize the likelihood with theta
        for (ig in seq_len(length(GT_uniq))) {
            theta[ig] <- sum(prob_mat * S1[[ig]]) / sum(prob_mat * SS[[ig]])
        }

    }
    if (verbose)
        print(paste("Total iterations:", it))

    ## return values
    prob_singlet <- prob_mat[,1:K_singlet]
    prob_doublet <- NULL
    if (check_doublet)
        prob_doublet <- prob_mat[, (K_singlet + 1):ncol(prob_mat)]

    assigned <- assign_cells_to_clones(prob_singlet,
                                       threshold = p_threshold)
    colnames(assigned) <- c("cell", "donor_id", "prob_max")
    if (check_doublet) {
        assigned$prob_doublet <- matrixStats::rowSums2(prob_doublet)
        assigned$donor_id[assigned$prob_doublet > 0.05] <- "doublet"
    } else
        assigned$prob_doublet <- NA
    assigned$n_vars <- n_vars
    assigned$donor_id[n_vars < n_vars_threshold] <- "unassigned"

    list("logLik" = logLik, "theta" = theta, "prob" = prob_singlet,
         "prob_doublet" = prob_doublet, "assigned" = assigned)
}



#' Gibbs sampler for donor deconvolution without genotypes.
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
donor_id_Gibbs <- function(A, D, K, prior0 = c(0.02, 40), prior1 = c(2.4, 2.4),
                           min_iter=500, max_iter=50000, buin_in=0.25){
    ##TODO: 1) detect doublet with BF, 2) support inputed Config
    ## 3) multiple chain for local optima

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
    singlet_logLik <- matrix(0, nrow = max_iter, ncol = M)
    doublet_logLik <- matrix(0, nrow = max_iter, ncol = M)

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
            logLik_all[it] <- logLik_all[it] +
                matrixStats::logSumExp(logLik_mat[j,], na.rm = T)
        }
        assign_all[it, ] <- Iden_mat %*% seq_len(K)

        ## singlet loglikelihood
        for (j in seq_len(M)) {
            singlet_logLik[it, j] <- matrixStats::logSumExp(logLik_mat[j,], na.rm = T) /
                sum(!is.na(logLik_mat[j,]))
        }
        ## doublet loglikelihood
        combn_idx <- utils::combn(K, 2)
        Conf_mat2 <- Conf_mat[,combn_idx[1,]] + Conf_mat[,combn_idx[2,]]
        Conf_mat2[which(Conf_mat2 > 1)] <- 1
        logLik_mat2 <- t(P0_mat) %*% (1 - Conf_mat2) + t(P1_mat) %*% Conf_mat2
        for (j in seq_len(nrow(logLik_mat2))) {
            doublet_logLik[it, j] <- matrixStats::logSumExp(logLik_mat2[j,], na.rm = T) /
                sum(!is.na(logLik_mat2[j,]))
        }

        # Update logLikelihood (TODO: add W0 and W1)
        # logLik_all[it] <- sum(log(rowSums(exp(logLik_mat), na.rm = TRUE)),
        #                       na.rm = TRUE)
        logLik_all[it] <- sum(singlet_logLik[it])

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

    ## Bayes factor with harmonic mean of the likelihoods (Newton & Raftery,1994)
    log10_BF <- rep(0, M)
    for (j in seq_len(M)) {
        log10_BF[j] <- (matrixStats::logSumExp(-singlet_logLik[n_buin:it, j]) -
                            matrixStats::logSumExp(-doublet_logLik[n_buin:it, j]))
    }

    return_list <- list("prob" = prob_mat, "theta0" = theta0, "theta1" = theta1,
                        "theta0_all" = as.matrix(theta0_all[1:it,]),
                        "theta1_all" = as.matrix(theta1_all[1:it,]),
                        "logLik" = logLik_all[1:it], "assign" = assign_all[1:it,],
                        "Config" = Config, "Config_all" = Config_all[1:it,],
                        "log10_BF" = log10_BF)
    return_list
}
