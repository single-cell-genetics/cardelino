# Tests for simulation methods
# library(cardelino); library(testthat); source("test-simulator.R")

context("test simulator")

data(simulation_input)

test_that("depths generator works as expected", {
    D1 <- sample_seq_depth(D_input, n_cells=500, n_sites=50, missing_rate=0.85)
    expect_is(D1, "matrix")
})

test_that("read counts generator works as expected", {
    D2 <- sample_seq_depth(D_input, n_cells=500, n_sites=nrow(tree_4clone$Z))
    simu <- sim_read_count(tree_4clone$Z, D2, Psi=NULL, cell_num=500)
    expect_is(simu, "list")
})

test_that("down sample variants in tree works as expected", {
    tree_lite <- sample_tree_SNV(tree_4clone, n_SNV=10)
    expect_is(tree_lite, "phylo")
})



context("test read_data.R")
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
