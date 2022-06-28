# Tests for a set of methods
# library(cardelino); library(testthat); source("test-misc.R")

context("test assessment.R")

test_that("mtx_to_df works as expected", {
    df <- mtx_to_df(matrix(seq(12), nrow = 3))
    expect_is(df, "data.frame")
})

test_that("rowMax works as expected", {
    matA <- matrix(sample(seq(12)), nrow = 3)
    rMax <- rowMax(matA)
    expect_is(rMax, "numeric")
})

test_that("rowArgmax works as expected", {
    matA <- matrix(sample(seq(12)), nrow = 3)
    rArg <- rowArgmax(matA)
    expect_is(rArg, "numeric")
})

test_that("colMatch works as expected", {
    matA <- matrix(sample(seq(12)), nrow = 3)
    col_idx <- sample(4)
    matB <- matA[, col_idx]
    cIdx <- colMatch(matB, matA)
    expect_is(cIdx, "numeric")
})


context("test read_data.R")
cell_vcf <- system.file("extdata", "cellSNP.cells.vcf.gz", package = "cardelino")


cell_data <- load_cellSNP_vcf(cell_vcf,
    max_other_allele = NULL,
    min_count = 0, min_MAF = 0
)

test_that("load_cellSNP_vcf works as expected", {
    expect_is(cell_data$A, "dgCMatrix")
    expect_is(cell_data$D, "dgCMatrix")
    expect_identical(dim(cell_data$D), dim(cell_data$A))
})


# donor_vcf <- system.file("extdata", "donors.donorid.vcf.gz", package = "cardelino")
# donor_GT <- load_GT_vcf(donor_vcf)
# rownames(donor_GT$GT) <- paste0("chr", rownames(donor_GT$GT)) #not always necessary

donor_GT <- load_GT_vcf(cell_vcf, na.rm = FALSE)

test_that("load_GT_vcf works as expected", {
    expect_is(donor_GT$GT, "matrix")
})


context("test tree_utils.R")
data("example_donor")
test_that("get_tree works as expected", {
    tree_obj <- get_tree(tree$Z)
    expect_is(tree_obj, "phylo")
})
