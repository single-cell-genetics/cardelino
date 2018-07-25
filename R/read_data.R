# Functions to read data into cardelino

#' Parse with single-cell variant data from VCF file for cardelino
#'
#' @param vcf_file path to VCF file from which to read data
#' @param filter_variants logical(1), should variants that do not have any read
#' coverage in any cell be filtered out?
#' @param filter_cells logical(1), should cells with fewer than a threshold
#' number of variants with read coverage be filtered out?
#' @param cell_nvars_threshold numeric(1), threshold for the number of variants
#' with non-zero read coverage; if \code{filter_cells} is \code{TRUE} then cells
#' with a number of variants with read coverage below this threshold will be
#' filtered out.
#' @param verbose logical(1), should messages be printed as function runs?
#' @param ... further arguments passed to \code{\link[vcfR]{read.vcfR}}
#'
#' @return a list with elements: \code{A}, a variant x cell matrix of read
#' counts supporting the alternative allele; \code{D}, a variant x cell matrix
#' of total read depth for each variant (set to \code{NA} if depth is zero);
#' \code{R}, a variant x cell matrix of the read counts supporting the reference
#' allele.
#'
#' @author Davis McCarthy
#'
#' @export
#'
#' @examples
#' vcf <- system.file("extdata", "cell_example.mpileup.vcf.gz",
#'                    package = "cardelino")
#' input_data <- parse_cell_vcf(vcf, filter_variants = TRUE)
#'
parse_cell_vcf <- function(vcf_file, filter_variants = TRUE, filter_cells = FALSE,
                      cell_nvars_threshold = 1.5, verbose = TRUE, ...) {
    vcf <- vcfR::read.vcfR(vcf_file, verbose = verbose, ...)
    ## get read count data
    dp <- vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE)
    ad <- vcfR::extract.gt(vcf, element = "AD")
    rf <- vcfR::masplit(ad, record = 1L, sort = FALSE)
    ad <- vcfR::masplit(ad, record = 2L, sort = FALSE)
    rownames(ad) <- rownames(dp) <- rownames(rf) <-
        paste0(vcf@fix[, 1], "_", vcf@fix[, 2])
    ## fix up dodgy entries - if depth is zero, set A to NA
    ad[which(dp == 0)] <- NA
    dp[which(dp == 0)] <- NA
    ## if depth is non-zero, but ad is NA, set ad to zero
    ad[(dp > 0) & is.na(ad)] <- 0

    ## filter variants with no cells
    if (filter_variants) {
        idx_var_use <- rowSums(dp, na.rm = TRUE) > 0.5
        if (!any(idx_var_use))
            stop("No variants are genotyped in at least one cell.")
        ad <- ad[idx_var_use,, drop = FALSE]
        dp <- dp[idx_var_use,, drop = FALSE]
        rf <- rf[idx_var_use,, drop = FALSE]
    }

    ## filter cells
    if (filter_cells) {
        nvars_genotyped <- colSums(dp > 0.5, na.rm = TRUE)
        if (all(nvars_genotyped < cell_nvars_threshold))
            stop("No cells have more than one variant with read coverage.")
        ad <- ad[, nvars_genotyped > cell_nvars_threshold, drop = FALSE]
        dp <- dp[, nvars_genotyped > cell_nvars_threshold, drop = FALSE]
        rf <- rf[, nvars_genotyped > cell_nvars_threshold, drop = FALSE]
    }
    list(A = ad, D = dp, R = rf)
}
