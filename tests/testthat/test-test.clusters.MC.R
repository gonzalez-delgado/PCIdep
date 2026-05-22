# ==============================================================================
# Test suite: test.clusters.MC
#
# Verifies the behaviour of test.clusters.MC(), which performs post-clustering
# selective inference for the difference between two cluster means after any
# user-supplied clustering algorithm, under a general matrix normal model.
#
# Because the truncation region cannot be characterised analytically for an
# arbitrary algorithm, it is always approximated by Monte Carlo: ndraws
# perturbed data sets are re-clustered and the proportion that preserves the
# selected clusters gives the p-value.  The standard error of that estimate is
# returned as stderr.
#
# Tests cover:
#   - Basic output structure: named list with pvalue, stat, stderr, clusters.
#   - Validity of output values: pvalue in [0,1], stat > 0, clusters length == n.
#   - Execution with a precomputed clustering vs. calling cl_fun internally.
#   - Input validation: cluster labels absent from cl, wrong number of labels,
#     cl_fun missing the NC argument when NC is provided.
#   - Optional arguments: return_Sigma, auxiliary Y, sample splitting.
# ==============================================================================

set.seed(456)

# make_mc_data() generates a 30 x 5 matrix with three clearly separated groups
# (means -4, 0, +4) together with identity row and column covariances.  The
# strong separation keeps the Monte Carlo acceptance rate high even with small
# ndraws, making tests fast without sacrificing coverage.
make_mc_data <- function(seed = 1) {
  set.seed(seed)
  n <- 30; p <- 5
  X <- rbind(
    matrix(rnorm(10 * p, mean = -4), nrow = 10),
    matrix(rnorm(10 * p, mean =  0), nrow = 10),
    matrix(rnorm(10 * p, mean =  4), nrow = 10)
  )
  Sigma <- diag(p)
  U     <- diag(n)
  list(X = X, Sigma = Sigma, U = U, n = n, p = p)
}

# km_fun is a minimal wrapper around kmeans() that matches the cl_fun interface
# expected by test.clusters.MC() (takes X and NC, returns integer cluster labels).
km_fun <- function(X, NC) {
  km <- kmeans(X, centers = NC, nstart = 5)
  km$cluster
}

# ---- Basic output structure with precomputed cl -------------------------

# The return value must contain pvalue, stat, stderr, and clusters.  Passing a
# precomputed cl avoids running cl_fun a second time during the test.
test_that("test.clusters.MC returns list with pvalue, stat, stderr, clusters", {
  d <- make_mc_data()
  cl <- km_fun(d$X, NC = 3)
  res <- suppressMessages(
    test.clusters.MC(
      X = d$X, U = d$U, Sigma = d$Sigma,
      clusters = c(1, 3),
      cl_fun = km_fun,
      NC = 3,
      cl = cl,
      ndraws = 100
    )
  )
  expect_true(all(c("pvalue", "stat", "stderr", "clusters") %in% names(res)))
})

# The p-value is estimated as a weighted Monte Carlo proportion and must lie in
# [0, 1].
test_that("test.clusters.MC pvalue is in [0, 1]", {
  d <- make_mc_data()
  cl <- km_fun(d$X, NC = 3)
  res <- suppressMessages(
    test.clusters.MC(
      X = d$X, U = d$U, Sigma = d$Sigma,
      clusters = sort(unique(cl))[c(1, 3)],
      cl_fun = km_fun,
      NC = 3,
      cl = cl,
      ndraws = 100
    )
  )
  expect_gte(res$pvalue, 0)
  expect_lte(res$pvalue, 1)
})

# The test statistic is the V-norm of the difference of cluster means and must
# be a single positive number.
test_that("test.clusters.MC stat is a positive numeric scalar", {
  d <- make_mc_data()
  cl <- km_fun(d$X, NC = 3)
  res <- suppressMessages(
    test.clusters.MC(
      X = d$X, U = d$U, Sigma = d$Sigma,
      clusters = sort(unique(cl))[c(1, 3)],
      cl_fun = km_fun,
      NC = 3,
      cl = cl,
      ndraws = 100
    )
  )
  expect_true(is.numeric(res$stat))
  expect_length(res$stat, 1)
  expect_gt(res$stat, 0)
})

# clusters must record the assignment of every observation, so its length must
# equal the number of rows in X.
test_that("test.clusters.MC clusters output has length n", {
  d <- make_mc_data()
  cl <- km_fun(d$X, NC = 3)
  res <- suppressMessages(
    test.clusters.MC(
      X = d$X, U = d$U, Sigma = d$Sigma,
      clusters = sort(unique(cl))[c(1, 3)],
      cl_fun = km_fun,
      NC = 3,
      cl = cl,
      ndraws = 100
    )
  )
  expect_length(res$clusters, d$n)
})

# ---- Runs cl_fun internally when cl is NULL -----------------------------

# When no precomputed clustering is supplied, test.clusters.MC() must call
# cl_fun(X, NC) internally and proceed to a valid p-value.
test_that("test.clusters.MC runs cl_fun when cl is NULL", {
  d <- make_mc_data()
  res <- suppressMessages(
    test.clusters.MC(
      X = d$X, U = d$U, Sigma = d$Sigma,
      clusters = c(1, 3),
      cl_fun = km_fun,
      NC = 3,
      cl = NULL,
      ndraws = 100
    )
  )
  expect_gte(res$pvalue, 0)
  expect_lte(res$pvalue, 1)
})

# ---- Input validation ---------------------------------------------------

# Cluster label 99 is not present in cl; the function must stop with an
# informative error referencing the clusters argument.
test_that("test.clusters.MC stops when clusters not found in cl", {
  d <- make_mc_data()
  cl <- km_fun(d$X, NC = 3)
  expect_error(
    suppressMessages(
      test.clusters.MC(
        X = d$X, U = d$U, Sigma = d$Sigma,
        clusters = c(1, 99),
        cl_fun = km_fun,
        NC = 3,
        cl = cl,
        ndraws = 50
      )
    ),
    regexp = "clusters"
  )
})

# The inference procedure compares exactly two clusters; passing three labels
# must be caught before any Monte Carlo sampling is performed.
test_that("test.clusters.MC stops when clusters has wrong length", {
  d <- make_mc_data()
  cl <- km_fun(d$X, NC = 3)
  expect_error(
    suppressMessages(
      test.clusters.MC(
        X = d$X, U = d$U, Sigma = d$Sigma,
        clusters = c(1, 2, 3),
        cl_fun = km_fun,
        NC = 3,
        cl = cl,
        ndraws = 50
      )
    ),
    regexp = "length 2"
  )
})

# If NC is provided but cl_fun does not accept an NC argument, the function
# cannot pass NC to it and must stop with a clear error.
test_that("test.clusters.MC stops when cl_fun lacks NC argument but NC is given", {
  d <- make_mc_data()
  no_nc_fun <- function(X) kmeans(X, centers = 3)$cluster
  expect_error(
    suppressMessages(
      test.clusters.MC(
        X = d$X, U = d$U, Sigma = d$Sigma,
        clusters = c(1, 3),
        cl_fun = no_nc_fun,
        NC = 3,
        ndraws = 50
      )
    ),
    regexp = "NC"
  )
})

# ---- Optional arguments -------------------------------------------------

# When return_Sigma = TRUE, the column covariance used in the test must appear
# in the returned list with the correct p x p dimensions.
test_that("test.clusters.MC includes Sigma when return_Sigma = TRUE", {
  d <- make_mc_data()
  cl <- km_fun(d$X, NC = 3)
  res <- suppressMessages(
    test.clusters.MC(
      X = d$X, U = d$U, Sigma = d$Sigma,
      clusters = sort(unique(cl))[c(1, 3)],
      cl_fun = km_fun,
      NC = 3,
      cl = cl,
      ndraws = 100,
      return_Sigma = TRUE
    )
  )
  expect_true("Sigma" %in% names(res))
  expect_equal(dim(res$Sigma), c(d$p, d$p))
})

# When Sigma is not known, it can be over-estimated from an independent sample Y
# with the same variables.  The test must still produce a valid p-value.
test_that("test.clusters.MC works with Sigma estimated from auxiliary Y", {
  d <- make_mc_data()
  set.seed(55)
  Y <- matrix(rnorm(d$n * d$p), nrow = d$n)
  cl <- km_fun(d$X, NC = 3)
  res <- suppressMessages(
    test.clusters.MC(
      X = d$X, U = d$U, Sigma = NULL, Y = Y,
      clusters = sort(unique(cl))[c(1, 3)],
      cl_fun = km_fun,
      NC = 3,
      cl = cl,
      ndraws = 100
    )
  )
  expect_gte(res$pvalue, 0)
  expect_lte(res$pvalue, 1)
})

# Sample splitting diverts half of X to estimate Sigma; a larger n is used here
# so both the estimation and clustering subsamples are adequately sized.
test_that("test.clusters.MC works with sample splitting", {
  set.seed(22)
  n <- 60; p <- 5
  X <- rbind(
    matrix(rnorm(20 * p, mean = -4), nrow = 20),
    matrix(rnorm(20 * p, mean =  0), nrow = 20),
    matrix(rnorm(20 * p, mean =  4), nrow = 20)
  )
  res <- suppressWarnings(suppressMessages(
    test.clusters.MC(
      X = X, U = NULL, Sigma = NULL,
      clusters = c(1, 3),
      cl_fun = km_fun,
      NC = 3,
      sample_split = TRUE,
      ndraws = 100
    )
  ))
  expect_gte(res$pvalue, 0)
  expect_lte(res$pvalue, 1)
})
