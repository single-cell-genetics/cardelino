# Tests for donor_id methods
# library(cardelino); library(testthat); source("test-donor_id.R")

context("test donor ID")
cell_vcf <- system.file("extdata", "cells.donorid.vcf.gz",
                        package = "cardelino")
donor_vcf <- system.file("extdata", "donors.donorid.vcf.gz", 
                         package = "cardelino")

cell_data <- load_cellSNP_vcf(cell_vcf, 
                              max_other_allele = NULL, 
                              min_count = 0, min_MAF = 0)
donor_GT <- load_GT_vcf(donor_vcf)
rownames(donor_GT$GT) <- paste0("chr", rownames(donor_GT$GT)) #not always necessary

test_that("load_cellSNP_vcf works as expected", {
    expect_is(cell_data$A, "dgCMatrix")
    expect_is(cell_data$D, "dgCMatrix")
    expect_identical(dim(cell_data$D), dim(cell_data$A))
})

test_that("load_GT_vcf works as expected", {
    expect_is(donor_GT$GT, "matrix")
})

test_that("vireo without known donor GT works as expected", {
    ids1 <- vireo(cell_data = cell_data, n_donor = 3, n_init = 2)
    print(table(ids1$assigned$donor_id))
    expect_is(ids1, "list")
})

test_that("vireo with known donor GT works as expected", {
    ids2 <- vireo(cell_data = cell_data, donor_data = donor_GT$GT)
    print(table(ids2$assigned$donor_id))
    expect_is(ids2, "list")
})

# test_that("vireo with known donor GT works as expected", {
#     mm <- match(rownames(donor_GT$GT), rownames(cell_data$A))
#     print(paste(sum(!is.na(mm)), "out of", length(mm), 
#                 "variants mathed between GT and cells."))
#     
#     idx <- which(!is.na(mm))
#     GT_use <- donor_GT$GT[idx, ]
#     idx_GP <- c()
#     for (ii in seq_len(ncol(donor_GT$GT))) {
#         idx_GP <- c(idx_GP, idx + nrow(donor_GT$GT) * (ii - 1))
#     }
#     GP_use <- donor_GT$GP[idx_GP, ]
#     
#     cell_data$A <- cell_data$A[mm[!is.na(mm)], ]
#     cell_data$D <- cell_data$D[mm[!is.na(mm)], ]
#     
#     GP_use <- GT_to_prob(GT_use)
#     GP_use <- (GP_use - 0.5) * 0.999 + 0.5
#     
#     
#     ## It seems harmful to use GP
#     ids3 <- vireo(cell_data = cell_data, n_donor = 3, GT_prior = GP_use, 
#                   K_amplify = 1)
#     print(table(ids3$assigned$donor_id))
#     expect_is(ids3, "list")
# })
