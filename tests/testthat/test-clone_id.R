# Tests for clone_id methods
# library(cardelino); library(testthat); source("test-clone_id.R")

context("test clone ID")
data("example_donor")

# test_that("Model selection as expected", {
#     doMC::registerDoMC(4)
#     `%dopar%` <- foreach::`%dopar%`
#     res_all <- foreach::foreach(ii = 2:5) %dopar% {
#         clone_id(A_clone, D_clone, n_clone = ii, 
#                  min_iter = 10000, max_iter = 10000, 
#                  prior1=c(45, 55), relabel = TRUE)
#     }
# 
#     n_clones <- seq(2,5)
#     DIC <- rep(0, 4)
#     for (i in seq_len(4)) {
#         DIC[i] <- res_all[[i]]$DIC$DIC
#     }
#     plot(n_clones, DIC, type = "b")
# })


test_that("binomial EM inference works as expected", {
    assignments_EM <- clone_id(A_clone, D_clone, Config = tree$Z, 
                               inference = "EM")
    expect_is(assignments_EM, "list")
})


assignments <- clone_id(A_clone, D_clone, Config = tree$Z, 
                        min_iter = 500, max_iter = 1000, 
                        relax_Config = TRUE, relabel=TRUE)

test_that("default inference works as expected", {
    expect_is(assignments, "list")
})


context("test plotting")

test_that("heatmap for assignment probability works as expected", {
    fig1 <- prob_heatmap(assignments$prob)
    expect_is(fig1, "ggplot")
})

test_that("pheatmap for variants probability across cells works as expected", {
    fig2 <- vc_heatmap(assignments$prob_variant, assignments$prob, tree$Z)
    expect_is(fig2, "pheatmap")
})

context("test assessment.R")
test_that("assign_scores, multiPRC, binaryPRC work as expected", {
    I_sim <- (assignments$prob == rowMax(assignments$prob))
    res <- assign_scores(assignments$prob, I_sim)
    expect_is(res, "list")
})

test_that("binaryROC work as expected", {
    I_sim <- (assignments$prob == rowMax(assignments$prob))
    res <- binaryROC(assignments$prob[,1], I_sim[, 1])
    expect_is(res, "list")
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

