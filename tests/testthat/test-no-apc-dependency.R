test_that("repository sources contain no active apc package dependency", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  candidates <- unlist(lapply(
    c("R", "tests", "vignettes", ".github"),
    function(path) list.files(file.path(root, path), recursive = TRUE,
                              full.names = TRUE,
                              pattern = "\\.(R|Rmd|ya?ml)$")
  ))
  candidates <- c(candidates, file.path(root, c("DESCRIPTION", "NAMESPACE",
                                                "README.md")))
  candidates <- candidates[file.exists(candidates)]
  candidates <- candidates[!basename(candidates) %in%
                             c("test-no-apc-dependency.R",
                               "r-checkrelease.yml")]
  text <- unlist(lapply(candidates, readLines, warn = FALSE), use.names = FALSE)
  active <- paste0(
    "library\\s*\\(\\s*apc\\s*\\)|",
    "require\\s*\\(\\s*apc\\s*\\)|",
    "requireNamespace\\s*\\(\\s*['\"]apc['\"]|",
    "apc::|apc\\.(data|fit|forecast)|github::cran/apc"
  )
  expect_false(any(grepl(active, text, perl = TRUE)))
})
