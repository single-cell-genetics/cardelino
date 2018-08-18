# Tests for donor_id methods
# library(cardelino); library(testthat); source("test-donor_id.R")

context("test donor ID")

data(example_donor)

cell_vcf <- system.file("extdata", "cells.donorid.vcf.gz", 
                        package = "cardelino")
GT_vcf <- system.file("extdata", "donors.donorid.vcf.gz", package = "cardelino")

test_that("default donor_id works as expected", {
    ids1 <- donor_id(cell_vcf, donor_vcf = GT_vcf)
    expect_is(ids1, "list")
})

test_that("donor_id without doublet detection works as expected", {
    ids2 <- donor_id(cell_vcf, donor_vcf = GT_vcf, check_doublet = FALSE)
    expect_is(ids2, "list")
})

test_that("donor_id without genotypes works as expected", {
    ids3 <- donor_id(cell_vcf, n_donor = 3)
    expect_is(ids3, "list")
})

test_that("donor_id without genotypes and doublet detection works as expected", 
{
    ids4 <- donor_id(cell_vcf, n_donor = 3, check_doublet = FALSE)
    expect_is(ids4, "list")
})