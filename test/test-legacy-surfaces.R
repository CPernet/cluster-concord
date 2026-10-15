testthat::test_that("superseded analysis, preparation, and report artifacts are retired", {
  retired <- c(
    "../src/R/Analysis_1_code.R",
    "../src/R/Analysis_2_code.R",
    "../src/R/Prep_data.qmd",
    "../src/R/Prep_HC.Rmd",
    "../src/R/CV_new_final.Rmd",
    "../src/R/Analysis_after_clustering.qmd",
    "../src/R/Simulation-based validation.Rmd"
  )
  testthat::expect_false(any(file.exists(retired)))
})
