# Tests for plotting functions
# library(cardelino); library(testthat); source("test-plotting.R")

context("test plotting")

data(example_donor)

test_that("plot tree works as expected", {
    fig3 <- plot_tree(tree, orient = "v")
    expect_is(fig3, "ggtree")
})

test_that("plot_config_diffs works as expected", {
    Config1 <- matrix(c(rep(0, 15), rep(1, 8), rep(0, 7), 
                        rep(1, 5), rep(0, 3), rep(1, 7)), ncol = 3)
    Config2 <- matrix(c(rep(0, 15), rep(1, 8), rep(1, 7), 
                        rep(0, 5), rep(1, 3), rep(1, 7)), ncol = 3)
    rownames(Config1) <- rownames(Config2) <- paste0("var", 1:nrow(Config1))
    colnames(Config1) <- colnames(Config2) <- paste0("clone", 1:ncol(Config1))
    fig4 <- plot_config_diffs(Config1, Config2) + pub.theme()
    expect_is(fig4, "ggplot")
})

