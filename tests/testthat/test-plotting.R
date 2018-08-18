# Tests for plotting functions
# library(cardelino); library(testthat); source("test-plotting.R")

context("test plotting")

data(example_donor)

assignments <- clone_id(A, D, C = tree$Z)

test_that("heatmap for assignment probability works as expected", {
    fig1 <- prob_heatmap(assignments$prob)
    expect_is(fig1, "ggplot")
})

test_that("pheatmap for variants probability across cells works as expected", {
    fig2 <- vc_heatmap(assignments$prob_variant, assignments$prob, tree$Z)
    expect_is(fig2, "list")
})

test_that("plot tree works as expected", {
    fig3 <- plot_tree(tree, orient = "v")
    expect_is(fig3, "ggtree")
})

