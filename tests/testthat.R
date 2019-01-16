Sys.setenv("R_TESTS" = "")
library(testthat)
library(cardelino)

test_check("cardelino")

# check test coverage with covr in Rstudio
# covr::package_coverage()
