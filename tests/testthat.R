Sys.setenv("R_TESTS" = "")
library(testthat)
library(cardelino)

test_check("cardelino")
