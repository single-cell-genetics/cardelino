# Functions to read data into cardelino

# parse_cell_vcf <- function(vcf_file, filter_variants = TRUE,
#                            filter_cells = FALSE, cell_nvars_threshold = 1.5,
#                            verbose = TRUE, ...) {
#     vcf <- vcfR::read.vcfR(vcf_file, verbose = verbose, ...)
#     ## get read count data
#     dp <- vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE)
#     ad <- vcfR::extract.gt(vcf, element = "AD")
#     rf <- vcfR::masplit(ad, record = 1L, sort = FALSE)
#     ad <- vcfR::masplit(ad, record = 2L, sort = FALSE)
#     rownames(ad) <- rownames(dp) <- rownames(rf) <-
#         paste0(vcf@fix[, 1], "_", vcf@fix[, 2])
#     ## fix up dodgy entries - if depth is zero, set A to NA
#     ad[which(dp == 0)] <- NA
#     dp[which(dp == 0)] <- NA
#     ## if depth is non-zero, but ad is NA, set ad to zero
#     ad[(dp > 0) & is.na(ad)] <- 0
#
#     ## filter variants with no cells
#     if (filter_variants) {
#         idx_var_use <- rowSums(dp, na.rm = TRUE) > 0.5
#         if (!any(idx_var_use))
#             stop("No variants are genotyped in at least one cell.")
#         ad <- ad[idx_var_use,, drop = FALSE]
#         dp <- dp[idx_var_use,, drop = FALSE]
#         rf <- rf[idx_var_use,, drop = FALSE]
#     }
#
#     ## filter cells
#     if (filter_cells) {
#         nvars_genotyped <- colSums(dp > 0.5, na.rm = TRUE)
#         if (all(nvars_genotyped < cell_nvars_threshold))
#             stop("No cells have more than one variant with read coverage.")
#         ad <- ad[, nvars_genotyped > cell_nvars_threshold, drop = FALSE]
#         dp <- dp[, nvars_genotyped > cell_nvars_threshold, drop = FALSE]
#         rf <- rf[, nvars_genotyped > cell_nvars_threshold, drop = FALSE]
#     }
#     list(A = ad, D = dp, R = rf)
# }


## Turn off the examples and imports for `read_vcf`
# #' @examples
# #' vcf <- read_vcf(system.file("extdata", "cells.donorid.vcf.gz",
# #'                 package = "cardelino"))
# #' @importFrom VariantAnnotation readVcf isSNV geno ref alt
# #'   genotypeToSnpMatrix
# #' @importFrom GenomeInfoDb seqlevelsStyle seqlengths
# #' @importFrom methods as is

#' Read a VCF file into R session
#'
#' @param vcf_file character(1), path to VCF file to read into R session as a
#' \code{\link[VariantAnnotation]{CollapsedVCF}} object
#' @param genome character(1), string indicating the genome build used in the
#' VCF file(s) (default: "GRCh37")
#' @param seq_levels_style character(1), string passed to
#' \code{\link[GenomeInfoDb]{seqlevelsStyle}} the style to use for
#' chromosome/contig names (default: "Ensembl")
#' @param verbose logical(1), should messages be printed as function runs?
#'
#' @return a vcf object
#'
#' @export
#'
read_vcf <- function(vcf_file, genome = "GRCh37",
    seq_levels_style = "Ensembl", verbose = TRUE) {
    ## Read in VCF from this sample
    if (verbose) {
          message("Reading VCF.")
      }
    vcf_sample <- VariantAnnotation::readVcf(vcf_file, genome)
    vcf_sample <- vcf_sample[VariantAnnotation::isSNV(vcf_sample)]
    if (length(vcf_sample) > 0) {
        GenomeInfoDb::seqlevelsStyle(vcf_sample) <- seq_levels_style
        new_snp_names <- gsub(
            "chr", "",
            gsub(
                ":", "_",
                gsub(
                    "_[ATCG]/[ATCG]", "",
                    names(vcf_sample)
                )
            )
        )
        names(vcf_sample) <- new_snp_names
    } else {
          stop("No SNVs present in VCF file.")
      }
    if (verbose) {
          message(
              "Read ", dim(vcf_sample)[1], " variants for ",
              dim(vcf_sample)[2], " samples."
          )
      }
    vcf_sample
}


## Turn off the examples and imports for `get_snp_matrices`
# #' @examples
# #' vcf_cell <- read_vcf(system.file("extdata", "cells.donorid.vcf.gz",
# #'                      package = "cardelino"))
# #' vcf_donor <-  read_vcf(system.file("extdata", "donors.donorid.vcf.gz",
# #'                        package = "cardelino"))
# #' snp_data <- get_snp_matrices(vcf_cell, vcf_donor)
# #' @importFrom GenomicRanges findOverlaps
# #' @importFrom S4Vectors queryHits subjectHits

#' Get SNP data matrices from VCF object(s)
#'
#' @param vcf_cell a \code{\link[VariantAnnotation]{CollapsedVCF}} object
#' containing variant data for cells
#' @param vcf_donor an optional \code{\link[VariantAnnotation]{CollapsedVCF}}
#' object containing genotype data for donors
#' @param verbose logical(1), should the function output verbose information as
#' it runs?
#' @param donors optional character vector providing a set of donors to use, by
#' subsetting the donors present in the \code{donor_vcf_file}; if \code{NULL}
#' (default) then all donors present in VCF will be used.
#'
#' @return a list containing
#' \code{A}, a matrix of integers. Number of alteration reads in SNP i cell j.
#' \code{D}, a matrix of integers. Number of reads depth in SNP i cell j.
#' \code{R}, a matrix of integers. Number of reference reads in SNP i cell j.
#' \code{GT_cells}, a matrix of integers for genotypes. The cell-SNP
#' configuration.
#' \code{GT_donors}, a matrix of integers for genotypes. The donor-SNP
#' configuration.
#'
#' @export
#'
get_snp_matrices <- function(vcf_cell, vcf_donor = NULL, verbose = TRUE,
    donors = NULL) {
    if (!methods::is(vcf_cell, "CollapsedVCF")) {
          stop("vcf_cell must be a CollapsedVCF object from the ",
               "VariantAnnotation package.")
      }
    vcf_cell <- GenomeInfoDb::sortSeqlevels(vcf_cell, X.is.sexchrom = TRUE)
    slengths_sample <- GenomeInfoDb::seqlengths(vcf_cell)
    if (!is.null(vcf_donor)) {
        if (!methods::is(vcf_donor, "CollapsedVCF")) {
              stop("vcf_donor must be a CollapsedVCF object from the ",
                   "VariantAnnotation package.")
          }
        ## filter sample VCF to those variants found in donor VCF
        if (!is.null(donors)) {
            if (sum(colnames(vcf_donor) %in% donors) < 1) {
                  stop("No donors in vcf_donor are found in the donors",
                       "argument supplied.")
              }
            vcf_donor <- vcf_donor[, colnames(vcf_donor) %in% donors]
        }
        vcf_donor <- GenomeInfoDb::sortSeqlevels(vcf_donor,
                                                 X.is.sexchrom = TRUE)
        GenomeInfoDb::seqlengths(vcf_donor) <-
            slengths_sample[GenomeInfoDb::seqlevels(vcf_donor)]
        ovlap <- GenomicRanges::findOverlaps(vcf_cell, vcf_donor)
        if (length(ovlap) < 1L) {
            stop("No variants overlapping in cell VCF and Donor VCF.")
        } else {
            if (verbose) {
                  message("Filtering cell VCF.")
              }
            vcf_cell <- vcf_cell[S4Vectors::queryHits(ovlap)]
            if (verbose) {
                  message("Filtering donor VCF.")
              }
            vcf_donor <- vcf_donor[S4Vectors::subjectHits(ovlap)]
            match_alleles <- unlist(
              VariantAnnotation::ref(vcf_cell) == VariantAnnotation::ref(vcf_donor) &
              VariantAnnotation::alt(vcf_cell) == VariantAnnotation::alt(vcf_donor)
            )
            if (sum(match_alleles) < 1L) {
                  stop("No variants with matching alleles in sample and ",
                       "donor VCFs")
              }
            vcf_cell <- vcf_cell[match_alleles]
            vcf_donor <- vcf_donor[match_alleles]
            if (verbose) {
                  message("Extracting donor SNP genotype matrix.")
              }
            sm_donor <- VariantAnnotation::genotypeToSnpMatrix(
                VariantAnnotation::geno(vcf_donor, "GT"),
                ref = VariantAnnotation::ref(vcf_donor),
                alt = VariantAnnotation::alt(vcf_donor)
            )
        }
    } else {
          sm_donor <- NULL
      }
    ## get snp matrices
    if (verbose) {
          message("Extracting cell SNP matrices.")
      }
    sm_sample <- VariantAnnotation::genotypeToSnpMatrix(
        VariantAnnotation::geno(vcf_cell, "GT"),
        ref = VariantAnnotation::ref(vcf_cell),
        alt = VariantAnnotation::alt(vcf_cell)
    )
    sm_sample_REF <- matrix(sapply(
        VariantAnnotation::geno(vcf_cell, "AD"),
        function(x) x[[1]]
    ), ncol = ncol(vcf_cell))
    sm_sample_ALT <- matrix(sapply(
        VariantAnnotation::geno(vcf_cell, "AD"),
        function(x) x[[2]]
    ), ncol = ncol(vcf_cell))
    sm_sample_ALT[is.na(sm_sample_ALT) & !is.na(sm_sample_REF)] <- 0
    sm_sample_DEP <- sm_sample_REF + sm_sample_ALT
    sm_sample_DEP[sm_sample_DEP == 0] <- NA
    sm_sample_REF[is.na(sm_sample_DEP)] <- NA
    sm_sample_ALT[is.na(sm_sample_DEP)] <- NA
    rownames(sm_sample_REF) <- rownames(sm_sample_ALT) <- rownames(sm_sample_DEP) <- rownames(vcf_cell)
    colnames(sm_sample_REF) <- colnames(sm_sample_ALT) <- colnames(sm_sample_DEP) <- colnames(vcf_cell)
    na_sample <- is.na(sm_sample_DEP)
    if (sum(!na_sample) < 1L) {
          stop("No variants with non-missing genotypes cells VCF and Donor VCF\n")
      }
    out <- list(
        A = sm_sample_ALT, D = sm_sample_DEP, R = sm_sample_REF,
        GT_cells = t(methods::as(sm_sample$genotypes, "numeric"))
    )
    if (!is.null(vcf_donor)) {
          out[["GT_donors"]] <- t(methods::as(sm_donor$genotypes, "numeric"))
      } else {
          out[["GT_donors"]] <- NULL
      }
    out
}


#' Load sparse matrices A and D from cellSNP VCF file with filtering SNPs
#'
#' @param vcf_file character(1), path to VCF file generated from cellSNP
#' @param max_other_allele maximum ratio of other alleles comparing to REF and
#' ALT alleles; for cellSNP vcf, we recommend 0.05
#' @param min_count minimum count across all cells, e.g., 20
#' @param min_MAF minimum minor allele fraction, e.g., 0.1
#' @param rowname_format the format of rowname: NULL is the default from vcfR,
#' short is CHROM_POS, and full is CHROM_POS_REF_ALT
#' @param keep_GL logical(1), if TRUE, check if GL (genotype probability) exists
#' it will be returned
#' @export
#'
#' @examples
#' vcf_file <- system.file("extdata", "cellSNP.cells.vcf.gz",
#'                         package = "cardelino")
#' input_data <- load_cellSNP_vcf(vcf_file)
load_cellSNP_vcf <- function(vcf_file, min_count = 0, min_MAF = 0,
    max_other_allele = NULL, rowname_format = "full",
    keep_GL = FALSE) {
    vcf_temp <- vcfR::read.vcfR(vcf_file)
    dp_full <- vcfR::extract.gt(vcf_temp, element = "DP", as.numeric = TRUE)
    ad_full <- vcfR::extract.gt(vcf_temp, element = "AD", as.numeric = TRUE)

    idx <- (rowSums(dp_full, na.rm = TRUE) >= min_count &
        ((rowSums(ad_full, na.rm = TRUE) /
            rowSums(dp_full, na.rm = TRUE)) >= min_MAF))

    idx <- idx & (!is.na(idx))

    if (!is.null(max_other_allele)) {
        dp_sum <- vcfR::extract.info(vcf_temp,
            element = "DP",
            as.numeric = TRUE
        )
        oth_sum <- vcfR::extract.info(vcf_temp,
            element = "OTH",
            as.numeric = TRUE
        )
        idx <- idx & (oth_sum / dp_sum < max_other_allele)
    }

    print(paste(sum(idx), "out of", length(idx), "SNPs passed."))

    A <- ad_full[idx, ]
    D <- dp_full[idx, ]
    A[is.na(A)] <- 0
    D[is.na(D)] <- 0

    if (!is.null(rowname_format)) {
        fix_val <- vcfR::getFIX(vcf_temp)
        if (rowname_format == "short") {
            var_ids <- paste0(fix_val[, 1], "_", fix_val[, 2])
        } else {
            var_ids <- paste0(
                fix_val[, 1], "_", fix_val[, 2],
                "_", fix_val[, 4], "_", fix_val[, 5]
            )
        }
    }
    row.names(A) <- row.names(D) <- var_ids[idx]

    GL_list <- list()
    if (keep_GL && sum(strsplit(vcf_temp@gt[1, 1], ":")$FORMAT == "GL") > 0) {
        GL_full <- vcfR::extract.gt(vcf_temp, element = "GL",
                                    as.numeric = FALSE)
        for (ii in seq_len(length(strsplit(GL_full[1, 1], ",")[[1]]))) {
            GL_tmp <- vcfR::masplit(GL_full, delim = ",", record = ii, sort = 0)
            GL_tmp[is.na(GL_tmp)] <- 0
            row.names(GL_tmp) <- var_ids
            GL_list[[ii]] <- Matrix::Matrix(GL_tmp[idx, ], sparse = TRUE)
        }
    }

    list(
        "A" = Matrix::Matrix(A, sparse = TRUE),
        "D" = Matrix::Matrix(D, sparse = TRUE),
        "GL" = GL_list
    )
}


#' Load genotype VCF into numeric values: 0, 1, or 2
#'
#' Note, the genotype VCF can be very big for whole genome. It would be more
#' efficient to only keep the wanted variants and samples. bcftools does such
#' jobs nicely.
#'
#' @param vcf_file character(1), path to VCF file for donor genotypes
#' @param rowname_format the format of rowname: NULL is the default from vcfR,
#' short is CHROM_POS, and full is CHROM_POS_REF_ALT
#' @param na.rm logical(1), if TRUE, remove the variants with NA values
#' @param keep_GP logical(1), if TRUE, check if GP (genotype probability) exists
#' it will be returned
#' @export
#'
#' @examples
#' vcf_file <- system.file("extdata", "cellSNP.cells.vcf.gz",
#'                         package = "cardelino")
#' GT_dat <- load_GT_vcf(vcf_file, na.rm = FALSE)
load_GT_vcf <- function(vcf_file, rowname_format = "full", na.rm = TRUE,
    keep_GP = TRUE) {
    GT_vcf <- vcfR::read.vcfR(vcf_file)

    ## variant ids
    if (!is.null(rowname_format)) {
        fix_val <- vcfR::getFIX(GT_vcf)
        if (rowname_format == "short") {
            var_ids <- paste0(fix_val[, 1], "_", fix_val[, 2])
        } else {
            var_ids <- paste0(
                fix_val[, 1], "_", fix_val[, 2],
                "_", fix_val[, 4], "_", fix_val[, 5]
            )
        }
    }

    ## GT values
    GT <- vcfR::extract.gt(GT_vcf, element = "GT")

    GT_num <- matrix(NA, nrow = nrow(GT), ncol = ncol(GT))
    colnames(GT_num) <- colnames(GT)
    rownames(GT_num) <- var_ids

    idx0 <- which(GT == "0/0" | GT == "0|0")
    idx2 <- which(GT == "1/1" | GT == "1|1")
    idx1 <- which(GT == "0/1" | GT == "1/0" |
        GT == "0|1" | GT == "1|0")

    GT_num[idx0] <- 0
    GT_num[idx1] <- 1
    GT_num[idx2] <- 2

    if (na.rm) {
        idx <- rowSums(is.na(GT_num)) == 0
    } else {
        idx <- seq_len(nrow(GT_num))
    }
    GT_num <- GT_num[idx, ]

    ## Genotype probability
    if (keep_GP && sum(strsplit(GT_vcf@gt[1, 1], ":")$FORMAT == "GP") > 0) {
        GP_val <- matrix(NA, nrow = length(GT_num), ncol = 3)
        rownames(GP_val) <- paste0(
            rep(colnames(GT_num), each = ncol(GT_num)), ":",
            rep(rownames(GT_num), times = ncol(GT_num))
        )
        colnames(GP_val) <- c("GT=0", "GT=1", "GT=2")
        GP <- vcfR::extract.gt(GT_vcf, element = "GP", as.numeric = FALSE)
        for (ii in seq_len(3)) {
            GP_val[, ii] <- vcfR::masplit(GP, delim = ",", record = ii,
                                          sort = 0)[idx]
        }
    } else {
        GP_val <- NULL
    }

    list("GT" = GT_num, "GP" = GP_val)
}


# #' Load sparse matrices A and D from cellSNP VCF file in HDF5 format
# #'
# #' @param filename character(1), path to HDF5 file generated from cellSNP
# #' @export
# #'
# load_vcf_h5 <- function(filename) {
#   if (!requireNamespace('hdf5r', quietly = TRUE)) {
#     stop("Please install hdf5r to read HDF5 files")
#   }
#   if (!file.exists(filename)) {
#     stop("File not found")
#   }
#   infile <- hdf5r::H5File$new(filename = filename, mode = 'r')

#   AD <- as.numeric(infile[['GenoINFO/AD']][])
#   DP <- as.numeric(infile[['GenoINFO/DP']][])
#   indices <- as.numeric(infile[['GenoINFO/indices']][])
#   indptr <- as.numeric(infile[['GenoINFO/indptr']][])
#   shp <- as.numeric(infile[['GenoINFO/shape']][])

#   A.mat <- Matrix::sparseMatrix(
#     i = indices + 1,
#     p = indptr,
#     x = as.numeric(x = AD),
#     dims = shp
#   )
#   D.mat <- Matrix::sparseMatrix(
#     i = indices + 1,
#     p = indptr,
#     x = as.numeric(x = DP),
#     dims = shp
#   )
#   rownames(A.mat) <- rownames(D.mat) <- infile[['samples']][]
#   colnames(A.mat) <- colnames(D.mat) <- infile[['variants']][]

#   list("A" = Matrix::t(A.mat), "D" = Matrix::t(D.mat))
# }
