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
