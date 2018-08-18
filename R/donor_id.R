## Donor deconvolution in multiplexed scRNA-seq.

#' Donor deconvolution of scRNA-seq data
#'
#' @param cell_vcf either character(1), path to a VCF file containing variant
#' data for cells, or a \code{\link[VariantAnnotation]{CollapsedVCF}} object
#' @param donor_vcf either character(1), path to a VCF file containing genotype
#' data for donors, or a \code{\link[VariantAnnotation]{CollapsedVCF}} object
#' @param n_donor integer(1), number of donors to infer if not given genotypes
#' @param model character(1), model for estimating cell donor ids: EM or Gibbs
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
#' @param ... arguments passed to \code{donor_id_EM} or \code{donor_id_Gibbs}
#'
#' @details This function reads in all elements of the provided VCF file(s) into
#' memory, so we highly recommend filtering VCFs to the minimal appropriate set
#' of variants (e.g. with the bcftools software) before applying them to this
#' function.
#'
#' @return a list with elements: \code{logLik}, log-likelihood of the fitted
#' model; \code{theta}, ; \code{GT}, a matrix of inferred genotypes for each
#' donor; \code{GT_doublet}, a matrix of inferred genotypes for each possible
#' doublet (pairwise combinations of donors); \code{prob}, a matrix of posterior
#' probabilities of donor identities for each cell; \code{prob_doublet}, a
#' matrix of posterior probabilities for each possible doublet for each cell;
#' \code{A}, a variant x cell matrix of observed read counts supporting the
#' alternative allele; \code{D}, a variant x cell matrix of observed read depth;
#' \code{assigned}, a data.frame reporting the cell-donor assignments with
#' columns "cell" (cell identifier), "donor_id" (inferred donor, or "doublet" or
#' "unassigned"), "prob_max" (the maximum posterior probability across donors),
#' "prob_doublet" (the probability that the cell is a doublet), "n_vars" (the
#' number of variants with non-zero read depth used for assignment).
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
donor_id <- function(cell_vcf, donor_vcf = NULL, n_donor=NULL,
                     check_doublet = TRUE, model = "EM",
                     n_vars_threshold = 10, p_threshold = 0.95,
                     verbose = TRUE, genome = "GRCh37",
                     seq_levels_style = "Ensembl", donors = NULL, ...) {
    model <- match.arg(model, c("EM", "Gibbs"))
    if (!(class(cell_vcf) %in% c("character", "CollapsedVCF")))
        stop("cell_vcf argument must be a character (filename) or a CollapsedVCF object from the VariantAnnotation package.")
    if (!methods::is(cell_vcf, "CollapsedVCF"))
        cell_vcf <- read_vcf(cell_vcf, genome = genome,
                             seq_levels_style = seq_levels_style,
                             verbose = verbose)
    if (is.null(donor_vcf)) {
        in_data <- get_snp_matrices(cell_vcf, NULL, verbose = verbose,
                                    donors = donors)
        if (is.null(n_donor))
            stop("Need number of donors if not given donor genotypes.")
    }
    else {
        if (!(class(donor_vcf) %in% c("character", "CollapsedVCF")))
            stop("cell_vcf argument must be a character (filename) or a CollapsedVCF object from the VariantAnnotation package.")
        if (!methods::is(donor_vcf, "CollapsedVCF"))
            donor_vcf <- read_vcf(donor_vcf, genome = genome,
                                 seq_levels_style = seq_levels_style,
                                 verbose = verbose)
        in_data <- get_snp_matrices(cell_vcf, donor_vcf, verbose = verbose,
                                    donors = donors)
        if (ncol(in_data$GT_donors) < 2 && check_doublet) {
            warning("Only one donor with genotypes provided, so cannot check for doublets.")
            check_doublet <- FALSE
        }
    }
    if (verbose)
        message("Donor ID using ", nrow(in_data$A), " variants")
    if (model == "Gibbs") {
        if (!is.null(donor_vcf))
            stop("If donor genotype is known, we recommend using EM model.")
        out <- donor_id_Gibbs(in_data$A, in_data$D, K = n_donor,
                        check_doublet = check_doublet,
                        verbose = verbose, ...)
    } else {
        out <- donor_id_EM(in_data$A, in_data$D, GT = in_data$GT_donors,
                        K = n_donor, check_doublet = check_doublet,
                        verbose = verbose, ...)
    }

    ## output data
    out$A <- in_data$A
    out$D <- in_data$D
    # out$GT <- in_data$GT_cells #out has GT output

    ## assign data frame
    n_vars <- matrixStats::colSums2(!is.na(out$D))
    assigned <- assign_cells_to_clones(out$prob, threshold = p_threshold)
    colnames(assigned) <- c("cell", "donor_id", "prob_max")
    if (check_doublet) {
        assigned$prob_doublet <- matrixStats::rowSums2(out$prob_doublet)
        assigned$donor_id[assigned$prob_doublet > 0.05] <- "doublet"
    } else {
        assigned$prob_doublet <- NA
    }
    assigned$n_vars <- n_vars
    assigned$donor_id[n_vars < n_vars_threshold] <- "unassigned"

    out$assigned <- assigned
    out
}

# ## test example
# ids <- donor_id("inst/extdata/cells.donorid.vcf.gz",
#          "inst/extdata/donors.donorid.vcf.gz")
# table(ids$assigned$donor_id)
# head(ids$assigned)
# ggplot(ids$assigned, aes(n_vars, prob_doublet, colour = donor_id)) +
#     geom_point(size = 3, alpha = 0.5) +
#     #ggthemes::scale_color_canva(palette = "Surf and turf")
#     theme_bw()
# ggplot(ids$assigned, aes(n_vars, prob_max, colour = donor_id)) +
#     geom_point(size = 3, alpha = 0.5) +
#     scale_x_log10() +
#     theme_bw()
# dplyr::filter(ids$assigned, donor_id == "doublet")
#
# sce <- readRDS("../hipsci-fibro/Data/processed/sce_merged_donors_cardelino_donorid_all_with_qc_labels.rds")
# coldat <- colData(sce)
# rm(sce)
#
# dplyr::filter(coldat, run_lane %in% c("22782_5", "22782_6", "22666_7")) %>%
#     ggplot(aes(x = log10_total_counts_endogenous, y = total_features,
#                colour = well_type)) +
#     geom_point() +
#     facet_wrap(~run_lane, ncol = 1)


# ids$assigned$sample_id <- ids$assigned$cell
# iddf <- dplyr::inner_join(ids$assigned, as.data.frame(coldat))
# iddf %>% dplyr::filter(donor_id == "doublet") %>%
#     dplyr::select(sample_id, well_condition, well_type, prob_max, prob_doublet) %>% as.data.frame
#
#
# cell_vcf <- read_vcf("inst/extdata/cells.donorid.vcf.gz")
# donor_vcf <- read_vcf("inst/extdata/donors.allfibro.vcf.gz")
# ids2_sing <- donor_id(cell_vcf, donor_vcf, check_doublet = FALSE)
# table(ids2_sing$assigned$donor_id)
# dplyr::filter(ids2_sing$assigned, grepl("laia", donor_id))
# ids2_doub <- donor_id(cell_vcf, donor_vcf, check_doublet = TRUE,
#                       donors = unique(ids2_sing$assigned$donor_id))
# table(ids2_doub$assigned$donor_id)
# ggplot(ids2_doub$assigned, aes(n_vars, prob_doublet, colour = donor_id)) +
#     geom_point(size = 3, alpha = 0.5) +
#     #ggthemes::scale_color_canva(palette = "Surf and turf")
#     theme_bw()
# ids2_doub$assigned$sample_id <- ids2_doub$assigned$cell
# id2df <- dplyr::inner_join(ids2_doub$assigned, as.data.frame(coldat))
# id2df %>% dplyr::filter(donor_id %in% c("doublet", "unassigned")) %>%
#     dplyr::select(cell, donor_id, well_condition, well_type, prob_max, prob_doublet) %>% as.data.frame
# id2df %>% dplyr::filter(well_type %in% c("minibulk_10", "control")) %>%
#     dplyr::select(cell, donor_id, well_condition, well_type, prob_max, prob_doublet) %>% as.data.frame
#
# ## ppca on cell genotypes
# pp <- pcaMethods::ppca(t(ids2_doub$A / ids2_doub$D))
# id2df$PPCA1 <- pp@scores[, 1]
# id2df$PPCA2 <- pp@scores[, 2]
# ggplot(id2df, aes(PPCA1, PPCA2, colour = donor_id)) +
#     geom_point(size = 3, alpha = 0.5) +
#     theme_bw()


#' EM algorithm for donor deconvolution with or without genotypes.
#'
#' @param A A matrix of integers. Number of alteration reads in SNP i cell j
#' @param D A matrix of integers. Number of reads depth in SNP i cell j
#' @param GT A matrix of integers for genotypes. The donor-SNP configuration.
#' @param K An integer. The number of donors to infer if not given GT.
#' @param gt_singlet A vector of integers. The candidate elements of GT. Default
#' is c(0, 1, 2); but can be others. Note, gt will be averaged for doublet.
#' When GT is given, gt_singlet will be adapted to GT.
#' @param check_doublet logical(1), if TRUE, check doublet, otherwise ignore.
#' @param min_iter A integer. The minimum number of iterations in EM algorithm.
#' @param max_iter A integer. The maximum number of iterations in EM algorithm.
#' The real iteration may finish earlier.
#' @param n_init A integer. The number of random initializations for EM
#' algorithm, which can be useful to avoid local optima if not given genotypes.
#' Default: 1 if given GT, 5 if not given GT.
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
#' \code{GT}, the input GT or a point estimate of genotype of donors. Note,
#' this may be not very accurate, especially for lowly expressed SNPs.
#' \code{GT_doublet}, the pair-wise doublet genotype based on GT.
#'
#' @import stats
#'
#' @export
#'
#' @examples
#' data(example_donor)
#' res <- donor_id_EM(A, D, GT = tree$Z[1:nrow(A),])
#' head(res$prob)
#'
donor_id_EM <- function(A, D, GT=NULL, K=NULL, gt_singlet=c(0, 1, 2),
                        check_doublet=TRUE, min_iter=10, max_iter=200,
                        n_init=NULL, logLik_threshold=1e-5, verbose=TRUE) {
    ## TODO: use multiple cores for multiple initializations

    ## Check input data
    if (nrow(A) != nrow(D) || ncol(A) != ncol(D)) {
        stop("A and D must have the same size.")
    }
    if (!is.null(GT)) {
        if (nrow(A) != nrow(GT)) {
            stop("nrow(A) and nrow(GT) must be the same.")
        }
    }

    ## preprocessing
    N <- nrow(A)        # number of SNPs
    M <- ncol(A)        # number of cells
    A[which(is.na(A))] <- 0
    D[which(is.na(D))] <- 0
    B <- D - A
    W_log <- sum(lchoose(D, A), na.rm = TRUE)  #log binomial coefficients
    GT_in <- GT # keep the input GT

    ## multiple initializations
    if (is.null(n_init)) {
        if (is.null(GT)) {n_init <- 3}
        else {n_init <- 2}
    }
    if (verbose) {cat(paste0(n_init, " random initializations is running.\n"))}
    GT_list <- theta_list <- prob_mat_list <- logLik_list <- list()
    for (i_init in seq_len(n_init)) {
        # foreach::foreach (i_init = seq_len(n_init)) %dopar% {

        # generate default values
        if (is.null(GT_in)) {
            if (is.null(K)) {
                stop("GT and K cannot both be NULL.")
            }
            GT <- matrix(sample(gt_singlet, size = N * K, prob = NULL,
                                replace = TRUE), nrow = N, ncol = K)
            update_GT <- TRUE
        } else {
            GT <- GT_in
            update_GT <- FALSE
            K <- ncol(GT) # number of donors
            gt_singlet <- sort(unique(c(GT)))
        }
        if (is.null(colnames(GT))) {colnames(GT) <- paste0("donor", seq_len(K))}

        # check doublet;
        if (check_doublet) {
            K1 <- K             # K1: number of singlet donors
            GT <- cbind(GT, get_doublet_GT(GT))
            K2 <- ncol(GT)      # K2: number of singlet and doublet donors

            gt_pair <- colMeans(utils::combn(gt_singlet, 2))
            gt_uniq <- c(gt_singlet, subset(gt_pair, !(gt_pair %in% gt_singlet)))
        }else{
            K2 <- K1 <- K
            gt_uniq <- gt_singlet
        }

        S1 <- S2 <- SS <- list()
        for (ig in seq_len(length(gt_uniq))) {
            S1[[ig]] <- t(A) %*% (GT == gt_uniq[ig])
            SS[[ig]] <- t(D) %*% (GT == gt_uniq[ig])
            S2[[ig]] <- SS[[ig]] - S1[[ig]]
        }

        ## random initialization for EM
        theta <- matrix(0.98 * gt_uniq / max(gt_uniq) + 0.01, ncol = 1)
        row.names(theta) <- paste0("GT=", gt_uniq)
        logLik <- 0
        logLik_new <- logLik + logLik_threshold * 2

        ## EM iterations
        for (it in seq_len(max_iter)) {
            # Check convergence
            if (it > min_iter) {
                if (is.na(logLik_new) || (logLik_new == -Inf)) {break}
                if ((logLik_new - logLik) < logLik_threshold) {break}
            }
            logLik <- logLik_new

            # E-step: update assignment posterior probability
            logLik_mat <- matrix(0, nrow = M, ncol = K2)
            for (ig in seq_len(length(gt_uniq))) {
                logLik_mat <- (logLik_mat + S1[[ig]] * log(theta[ig]) +
                                   S2[[ig]] * log(1 - theta[ig]))
            }

            logLik_vec <- rep(NA, nrow(logLik_mat))
            for (i in seq_len(nrow(logLik_mat))) {
                logLik_vec[i] <- matrixStats::logSumExp(logLik_mat[i,],
                                                        na.rm = TRUE)
            }
            logLik_new <- sum(logLik_vec, na.rm = TRUE) + W_log
            logLik_mat_amplify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
            prob_mat <- exp(logLik_mat_amplify) / rowSums(exp(logLik_mat_amplify))

            # M-step: maximumize the likelihood with theta and GT
            # update theta, but use the default for longer initialization
            if (it > min_iter / 2) {
                for (ig in seq_len(length(gt_uniq))) {
                    if (!is.na(sum(prob_mat * SS[[ig]])) &&
                        sum(prob_mat * SS[[ig]]) > 0) {
                        theta[ig] <- (sum(prob_mat * S1[[ig]]) /
                                          sum(prob_mat * SS[[ig]]))
                        if (theta[ig] < 0.01)
                            theta[ig] <- 0.01
                        if (theta[ig] > 0.98)
                            theta[ig] <- 0.98
                    }
                }
            }

            # update GT
            if (update_GT) {
                S1_gt <- A %*% prob_mat[, 1:K1, drop = FALSE]
                S2_gt <- B %*% prob_mat[, 1:K1, drop = FALSE]
                GT_prob <- matrix(0, nrow = N * K1, ncol = length(gt_singlet))
                for (ig in seq_len(ncol(GT_prob))) {
                    GT_prob[, ig] <- S1_gt * log(theta[ig]) +
                        S2_gt * log(1 - theta[ig])
                }
                for (ik in seq_len(length(GT[, 1:K1]))) {
                    GT[ik] <- gt_uniq[which.max(GT_prob[ik,, drop = FALSE])]
                }
                if (check_doublet) {GT[, (K1 + 1):K2] <-
                    get_doublet_GT(GT[, 1:K1, drop = FALSE])}

                for (ig in seq_len(length(gt_uniq))) {
                    S1[[ig]] <- t(A) %*% (GT == gt_uniq[ig])
                    SS[[ig]] <- t(D) %*% (GT == gt_uniq[ig])
                    S2[[ig]] <- SS[[ig]] - S1[[ig]]
                }
            }
        }

        if (verbose)
            cat(sprintf("Total iterations: %d; logLik: %.2f\n", it, logLik))

        GT_list[[i_init]] <- GT
        theta_list[[i_init]] <- theta
        logLik_list[[i_init]] <- logLik
        prob_mat_list[[i_init]] <- prob_mat
    }

    ## choose the initialization with highest logLik
    max_idx <- which.max(logLik_list)
    GT <- GT_list[[max_idx]]
    theta <- theta_list[[max_idx]]
    logLik <- logLik_list[[max_idx]]
    prob_mat <- prob_mat_list[[max_idx]]

    ## return values
    prob_singlet <- prob_mat[, 1:K1, drop = FALSE]
    prob_doublet <- GT_doublet <- NULL
    if (check_doublet) {
        prob_doublet <- prob_mat[, (K1 + 1):K2, drop = FALSE]
        GT_doublet <- GT[, (K1 + 1):K2, drop = FALSE]
    }

    return_list <- list("logLik" = logLik,
                        "theta" = theta,
                        "GT" = GT[, 1:K1, drop = FALSE],
                        "GT_doublet" = GT_doublet,
                        "prob" = prob_singlet,
                        "prob_doublet" = prob_doublet)
    return_list
}



#' Gibbs sampler for donor deconvolution without genotypes.
#'
#' @param A A matrix of integers. Number of alteration reads in SNP i cell j
#' @param D A matrix of integers. Number of reads depth in SNP i cell j
#' @param K An integer. The number of donors to infer if not given GT.
#' @param gt_singlet A vector of integers. The candidate elements of GT. Default
#' is c(0, 1, 2); but can be others. Note, gt will be averaged for doublet.
#' @param check_doublet logical(1), if TRUE, check doublet, otherwise ignore.
#' @param min_iter A integer. The minimum number of iterations in the Gibbs
#' sampling. The real iteration may be longer utile the convergence.
#' @param max_iter A integer. The maximum number of iterations in the Gibbs
#' sampling, even haven't passed the convergence diagnosis
#' @param buin_in A float between 0 and 1. The fraction for buil-in in MCMC
#' samplers.
#' @param EM_initial logical(1), if TRUE, use the EM estimate for a warm
#' initialization on GT and theta; otherwise random GT.
#' @param verbose logical(1), If TRUE, output verbose information when running.
#' @param ... arguments passed to \code{donor_id_EM}
#'
#' @details Users should typically use \code{\link{donor_id}} rather than this
#' lower-level function.
#'
#' @return a list containing
#' \code{logLik_samples}, the log likelihood for all samples.
#' \code{theta}, a vector denoting the binomial parameters for each genotype;
#' mean of all samples in \code{theta_samples}.
#' \code{theta_samples}, all samples of theta vector.
#' \code{prob}, a matrix of posterior probability of cell assignment to donors.
#' The summary may less than 1, as there are some probabilities go to doublets.
#' \code{prob_doublet}, a matrix of posterior probability of cell assignment to
#' each inter-donor doublet.
#' \code{prob_samples}, all samples of singlet probability and doublet
#' probability. Each column should be reshaped into a M by K2 matrix.
#' \code{GT_prob}, the posterior distribution of genotypes for each donors.
#' It has size of (N*K, length(gt_singlet)); each row should be reshaped into
#' a N by K matrix via matrix(GT_prob[, t], nrow = N, ncol = K).
#' \code{GT_samples}, all samples of donor genotype configuration matrices.
#'
#' @import matrixStats
#'
#' @export
#' @examples 
#' data(example_donor)
#' res <- donor_id_Gibbs(A, D, K = 4)
#' head(res$prob)
#'
donor_id_Gibbs <- function(A, D, K, gt_singlet=c(0, 1, 2), check_doublet=TRUE,
                           min_iter=400, max_iter=1000, buin_in=0.25,
                           EM_initial=TRUE, verbose=TRUE, ...) {
    ##TODO: 1) BF seems not practically suitable; sparse matrix may be useful;
    ## future: support informative prior, e.g., bulk RNA-seq for donor genotype,
    ## and imputed genotype.

    ## Check input data
    if (nrow(A) != nrow(D) || ncol(A) != ncol(D)) {
        stop("A and D must have the same size.")
    }
    ## preprocessing
    N <- nrow(A)        # number of SNPs
    M <- ncol(A)        # number of cells
    A[which(is.na(A))] <- 0
    D[which(is.na(D))] <- 0
    B <- D - A
    W_log <- sum(lchoose(D, A), na.rm = TRUE)  #log binomial coefficients

    ## generate default values and check doublet
    if (EM_initial) {
        if (verbose) {cat("Running donor_id_EM for warm initialization.\n")}
        res_EM <- donor_id_EM(A, D, GT = NULL, K = K, gt_singlet = gt_singlet,
                        check_doublet = check_doublet, verbose = verbose, ...)
        GT <- res_EM$GT[, 1:K]
        theta_default <- res_EM$theta
    } else {
        GT <- matrix(sample(gt_singlet, size = N * K, prob = NULL,
                            replace = TRUE), nrow = N, ncol = K)
        theta_default <- NULL
    }
    ## check doubet or not
    if (check_doublet) {
        K1 <- K
        GT <- cbind(GT, get_doublet_GT(GT))
        K2 <- ncol(GT)

        gt_pair <- colMeans(utils::combn(gt_singlet, 2))
        gt_uniq <- c(gt_singlet, subset(gt_pair, !(gt_pair %in% gt_singlet)))
    } else {
        K2 <- K1 <- K
        gt_uniq <- gt_singlet
    }
    prior1 = rep(1, length(gt_uniq))
    prior2 = rep(1, length(gt_uniq))
    if (is.null(theta_default)) {
        theta_default <- 0.98 * gt_uniq / max(gt_uniq) + 0.01
    }

    S1 <- S2 <- SS <- list()
    for (ig in seq_len(length(gt_uniq))) {
        S1[[ig]] <- t(A) %*% (GT == gt_uniq[ig])
        SS[[ig]] <- t(D) %*% (GT == gt_uniq[ig])
        S2[[ig]] <- SS[[ig]] - S1[[ig]]
    }

    ## variables to sample
    GT_all <- matrix(0, nrow = max_iter, ncol = N * K1)
    prob_all <- matrix(0, nrow = max_iter, ncol = M * K2)
    logLik_all <- matrix(0, nrow = max_iter, ncol = 1)
    theta_all <- matrix(theta_default, nrow = max_iter, ncol = length(gt_uniq),
                        byrow = TRUE)

    ## Gibbs sampling
    if (verbose) {cat("Gibbs sampling starts.\n")}
    for (it in 2:max_iter) {
        # Sample assignment
        logLik_mat <- matrix(0, nrow = M, ncol = K2)
        for (ig in seq_len(length(gt_uniq))) {
            logLik_mat <- (logLik_mat + S1[[ig]] * log(theta_all[it, ig]) +
                            S2[[ig]] * log(1 - theta_all[it, ig]))
        }

        logLik_vec <- rep(NA, nrow(logLik_mat))
        for (i in seq_len(nrow(logLik_mat))) {
            logLik_vec[i] <- matrixStats::logSumExp(logLik_mat[i,],
                                                    na.rm = TRUE)
        }
        logLik_all[it] <- sum(logLik_vec, na.rm = TRUE) + W_log
        logLik_mat_amplify <- logLik_mat - matrixStats::rowMaxs(logLik_mat)
        prob_mat <- exp(logLik_mat_amplify) / rowSums(exp(logLik_mat_amplify))
        prob_all[it, ] <- prob_mat

        # Sample thetas
        if (it > 10) {
            for (ig in seq_len(length(gt_uniq))) {
                S1_theta = sum(prob_mat * S1[[ig]])
                S2_theta = sum(prob_mat * S2[[ig]])
                theta_all[it, ig] <- stats::rbeta(1, prior1[ig] + S1_theta,
                                                      prior2[ig] + S2_theta)
            }
        }

        # Sample genotype configuration
        S1_gt <- A %*% prob_mat[, 1:K1]
        S2_gt <- B %*% prob_mat[, 1:K1]
        GT_prob <- matrix(0, nrow = N * K1, ncol = length(gt_singlet))
        for (ig in seq_len(ncol(GT_prob))) {
            GT_prob[, ig] <- S1_gt * log(theta_all[it,ig]) +
                             S2_gt * log(1 - theta_all[it,ig])
        }
        GT_prob <- GT_prob - matrixStats::rowMaxs(GT_prob)
        GT_prob <- exp(GT_prob) / rowSums(exp(GT_prob))
        for (ik in seq_len(nrow(GT_prob))) {
            GT_all[it, ik] <- sample(gt_singlet, size = 1, prob = GT_prob[ik, ])
        }
        GT[, 1:K1] <- matrix(GT_all[it, ], nrow = N, ncol = K1)
        colnames(GT)[1:K1] <- paste0("donor", seq_len(K1))
        if (check_doublet) {GT <- cbind(GT[, 1:K1], get_doublet_GT(GT[, 1:K1]))}
        # if (check_doublet) {GT[, (K1+1) : K2] <- get_doublet_GT(GT[, 1:K1])}

        for (ig in seq_len(length(gt_uniq))) {
            S1[[ig]] <- t(A) %*% (GT == gt_uniq[ig])
            SS[[ig]] <- t(D) %*% (GT == gt_uniq[ig])
            S2[[ig]] <- SS[[ig]] - S1[[ig]]
        }

        # Check convergence
        if (it >= min_iter && it %% 100 == 0) {
            is_converged <- Geweke_Z(prob_all[1:it, 1:K1]) <= 2
            if (verbose) {cat(sprintf("%d iterations: %.1f%% cells converged.\n",
                                     it, mean(is_converged) * 100))}
            if (sum(is_converged) >= (length(is_converged) - 1)) {break}
        }
    }

    ## return values
    n_buin = ceiling(it * buin_in)
    prob_mat <- matrix(colMeans(prob_all[n_buin:it, ]), nrow = M, ncol = K2)
    row.names(prob_mat) <- colnames(A)
    colnames(prob_mat) <- colnames(GT)
    colnames(theta_all) <- paste0("GT=", gt_uniq)
    theta <- colMeans(theta_all[n_buin:it, ])

    GT_prob <- matrix(0, nrow = N * K1, ncol = length(gt_singlet))
    for (ig in seq_len(ncol(GT_prob))) {
        GT_prob[, ig] <- colMeans(GT_all[n_buin:it, ] == gt_singlet[ig])
    }
    prob_singlet <- prob_mat[, 1:K1]
    prob_doublet <- NULL
    if (check_doublet) {prob_doublet <- prob_mat[, (K1 + 1):K2]}

    return_list <- list("prob" = prob_singlet,
                        "prob_doublet" = prob_doublet,
                        "prob_samples" = prob_all,
                        "theta" = theta,
                        "theta_samples" = theta_all[1:it, ],
                        "GT_prob" = GT_prob[1:(N * K1), ],
                        "GT_samples" = GT_all[1:it, ],
                        "logLik_samples" = logLik_all[1:it])
    return_list
}

#' Generate genotype for doublets
#' @param GT A matrix of genotype for singlets
#' @return \code{GT_doublet}, a matrix genotype for doublets
#' @example 
#' GT <- matrix(sample(c(0,1,2), 150, replace = TRUE), nrow = 50)
#' GT_doublet <- get_doublet_GT(GT)
#'
get_doublet_GT <- function(GT) {
    donor_ids <- colnames(GT)
    if (is.null(donor_ids)) { donor_ids <- paste0("donor", seq_len(ncol(GT))) }

    combn_idx <- utils::combn(ncol(GT), 2)
    GT_doublet <- (GT[,combn_idx[1,]] + GT[,combn_idx[2,]]) / 2
    colnames(GT_doublet) <- paste0(donor_ids[combn_idx[1,]], ",",
                                   donor_ids[combn_idx[2,]])
    GT_doublet
}
