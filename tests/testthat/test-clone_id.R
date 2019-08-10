# Tests for clone_id methods
# library(cardelino); library(testthat); source("test-clone_id.R")

context("test clone ID")
data("example_donor")

test_that("default inference works as expected", {
    assignments <- clone_id(A_clone, D_clone, Config = tree$Z, min_iter=1000)
    expect_is(assignments, "list")
})

test_that("binomial EM inference works as expected", {
    assignments_EM <- clone_id(A_clone, D_clone, Config = tree$Z, inference = "EM")
    expect_is(assignments_EM, "list")
})

# test_that("Bernoulli Gibbs inference works as expected", {
#     assignments_bern <- clone_id(A_clone, D_clone, Config = tree$Z,
#                                  model = "Bernoulli")
#     expect_is(assignments_bern, "list")
# })
# 
# test_that("Bernoulli EM inference works as expected", {
#     assignments_bern_EM <- clone_id(A_clone, D_clone, Config = tree$Z,
#                                     model = "Bernoulli", inference = "EM")
#     expect_is(assignments_bern_EM, "list")
# })

