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

# ---- n_preserved correctness --------------------------------------------

# Clusters with means at -100, +100, and 0 give stat_V ≈ 1000, which is far
# larger than sd(phi) ≈ sqrt(5) ≈ 2.24.  All ndraws phi values are therefore
# positive.  A threshold-based cl_fun that assigns labels purely by whether the
# first feature is below -50, above +50, or in between is immune to the tiny
# perturbations (at most ~0.3 units) induced by phi ≈ stat_V.  Consequently
# preserve.cl() returns TRUE for all ndraws iterations, so n_preserved must
# equal ndraws exactly.
test_that("n_preserved equals ndraws when all phi > 0 and all perturbations preserve clusters", {
  n <- 30; p <- 5
  X <- rbind(
    matrix(-100, nrow = 10, ncol = p),
    matrix( 100, nrow = 10, ncol = p),
    matrix(   0, nrow = 10, ncol = p)
  )
  U     <- diag(n)
  Sigma <- diag(p)
  cl    <- c(rep(1L, 10), rep(2L, 10), rep(3L, 10))

  # Deterministic cl_fun: assigns labels by thresholding the first feature.
  # The perturbations are ~0.1 * |phi - stat_V| << 50, so the label never changes.
  cl_fun_thresh <- function(X, NC) {
    ifelse(X[, 1] < -50, 1L, ifelse(X[, 1] > 50, 2L, 3L))
  }

  res <- suppressMessages(
    test.clusters.MC(
      X = X, U = U, Sigma = Sigma,
      clusters = c(1L, 2L),
      cl_fun = cl_fun_thresh,
      NC = 3, cl = cl,
      ndraws = 20,
      sample_split = FALSE
    )
  )

  expect_equal(res$n_preserved, 20L)
})

# A cl_fun that ignores the data entirely and always returns the same fixed
# partition guarantees preserve.cl() is TRUE for every phi >= 0.  By setting
# the same seed before the function call and before reproducing phi here, we
# obtain the exact phi vector used internally and can compute the expected
# n_preserved as sum(phi >= 0) without re-running the Monte Carlo loop.
#
# Note: with sample_split = FALSE and Sigma / cl supplied, setup.model()
# consumes no random numbers, so stats::rnorm() inside the function draws from
# the same state as the independently reproduced phi below.
test_that("n_preserved equals sum(phi >= 0) when cl_fun always preserves clusters", {
  n <- 30; p <- 5
  set.seed(7)
  X     <- matrix(rnorm(n * p), nrow = n, ncol = p)
  U     <- diag(n)
  Sigma <- diag(p)
  cl    <- c(rep(1L, 10), rep(2L, 10), rep(3L, 10))

  # Data-invariant cl_fun: returns the same partition for any input matrix.
  # preserve.cl(cl, cl_fun_fixed(Xphi), c(1,2)) is therefore always TRUE.
  cl_fun_fixed <- function(X, NC) rep(seq_len(NC), each = nrow(X) %/% NC)

  ndraws <- 50

  set.seed(99)
  res <- suppressMessages(
    test.clusters.MC(
      X = X, U = U, Sigma = Sigma,
      clusters = c(1L, 2L),
      cl_fun = cl_fun_fixed,
      NC = 3, cl = cl,
      ndraws = ndraws,
      sample_split = FALSE
    )
  )

  # Reproduce stat_V exactly as the function computes it:
  #   norm2U_nu = 1/n1 + 1/n2 = 0.2  (U = I, n1 = n2 = 10)
  #   V_g1g2    = 0.2 * Sigma = 0.2 * I_p
  #   stat_V    = sqrt(diff_means %*% (5 * I_p) %*% t(diff_means))
  diff_means_vec <- colMeans(X[1:10, ]) - colMeans(X[11:20, ])
  stat_V_rep     <- sqrt(sum(diff_means_vec^2) / 0.2)
  phi_sd         <- sqrt(sum(diag(Sigma)))   # sqrt(p) = sqrt(5)

  set.seed(99)
  phi_rep <- rnorm(ndraws, stat_V_rep, phi_sd)

  expect_equal(res$n_preserved, sum(phi_rep >= 0))
})
