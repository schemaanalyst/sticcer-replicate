#options(warn=-1)
#source("R/functions.R")

library(testthat)
context("testing_functions_file")

test_that("sorting a DF by schema name", {
  x <- data.frame("schema" = c("c","a","b"))
  x_sorted <- sortSchema(x)
  expect_true(all.equal(x_sorted$schema, c("\\aForTable", "\\bForTable", "\\cForTable")))
})

test_that("Testing assigement of sigifiance in data frame", {
  x <- sigBetweenTwoValues("R", 0, 1, 0.01, "large")
  expect_true(x == "$^{\\ast\\ast\\ast}$\\textbf{1}")
  x <- sigBetweenTwoValues("L", 0, 1, 0.01, "large")
  expect_true(x == "0")
})

test_that("mutationanalysis.dat is over or equal to 10200 rows", {
  results_path <- "../experiment_results/"
  mutationanalysis <- readr::read_csv(
    file.path(results_path, 
              dir(path = results_path, pattern = "mutationanalysis.dat")
              )
    )
  expect_true(nrow(mutationanalysis) >= 10200)
})

