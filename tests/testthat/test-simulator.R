# Tests for simulation methods
# library(cardelino); library(testthat); source("test-simulator.R")

context("test simulator")

data(simulation_input)

test_that("depths generator works as expected", {
    D1 <- sample_seq_depth(D_input, n_cells = 500, n_sites = 50, missing_rate = 0.85)
    expect_is(D1, "matrix")
})

test_that("read counts generator works as expected", {
    D2 <- sample_seq_depth(D_input, n_cells = 500, n_sites = nrow(tree_4clone$Z))
    simu <- sim_read_count(tree_4clone$Z, D2, Psi = NULL, cell_num = 500)
    expect_is(simu, "list")
})

test_that("down sample variants in tree works as expected", {
    tree_lite <- sample_tree_SNV(tree_4clone, n_SNV = 10)
    expect_is(tree_lite, "phylo")
})
