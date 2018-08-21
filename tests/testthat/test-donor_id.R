# Tests for donor_id methods
# library(cardelino); library(testthat); source("test-donor_id.R")

context("test donor ID")

data(example_donor, package = "cardelino")

cell_vcf <- system.file("extdata", "cells.donorid.vcf.gz",
                        package = "cardelino")
GT_vcf <- system.file("extdata", "donors.donorid.vcf.gz", package = "cardelino")

test_that("default donor_id works as expected", {
    expect_warning(ids1 <- donor_id(cell_vcf, donor_vcf = GT_vcf),
                   regexp = "non-diploid")
    expect_is(ids1, "list")
})

test_that("donor_id without doublet detection works as expected", {
    expect_warning(ids2 <- donor_id(cell_vcf, donor_vcf = GT_vcf,
                                    check_doublet = FALSE),
                   regexp = "non-diploid")
    expect_is(ids2, "list")
})

test_that("donor_id without genotypes works as expected", {
    expect_warning(ids3 <- donor_id(cell_vcf, n_donor = 3),
                   regexp = "non-diploid")
    expect_is(ids3, "list")
})

test_that("donor_id without genotypes and without doublet detection works as expected",
{
    expect_warning(ids4 <- donor_id(cell_vcf, n_donor = 3,
                                    check_doublet = FALSE),
                   regexp = "non-diploid")
    expect_is(ids4, "list")
})

test_that("donor_id_Gibbs without genotypes works as expected", {
    expect_warning(ids5 <- donor_id(cell_vcf, n_donor = 3, model = "Gibbs"),
                   regexp = "non-diploid")
    expect_is(ids5, "list")
})

test_that("donor_id_Gibbs without genotypes or doublet detection works as expected",
{
    expect_warning(ids6 <- donor_id(cell_vcf, n_donor = 3,
                                    check_doublet = FALSE, model = "Gibbs"),
                   regexp = "non-diploid")
    expect_is(ids6, "list")
})
