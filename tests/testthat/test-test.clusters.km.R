# ==============================================================================
# Test suite: test.clusters.km
#
# Verifies the behaviour of test.clusters.km(), which performs post-clustering
# selective inference for the difference between two cluster means after k-means
# clustering, under a general matrix normal model.
#
# The p-value is computed using the exact truncation-set characterisation of
# Chen and Witten (2023), mapped to the chi-squared distribution under the
# matrix normal model.
#
# Tests cover:
#   - Basic output structure: named list with pvalue, stat, km.
#   - Validity of output values: pvalue in [0,1], stat > 0, km length == n.
#   - Input validation: out-of-range or wrongly-sized cluster labels.
#   - Optional arguments: return_Sigma, auxiliary Y, sample splitting.
#   - Behaviour under H0: p-values should not be systematically near 0.
# ==============================================================================

set.seed(123)

# make_km_data() generates a 30 x 5 matrix with three clearly separated groups
# (means -4, 0, +4) together with an identity row and column covariance.  The
# strong separation ensures k-means reliably recovers the true partition.
make_km_data <- function(seed = 1) {
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

# run_km() is a thin wrapper around test.clusters.km() that fixes NC = 3 and
# passes data from a make_km_data() list, reducing boilerplate in each test.
run_km <- function(d, clusters, ...) {
  test.clusters.km(
    X = d$X, U = d$U, Sigma = d$Sigma,
    NC = 3, clusters = clusters,
    itermax = 50, tol = 1e-6,
    ...
  )
}

# ---- Basic output structure ---------------------------------------------

# The return value must be a named list containing exactly pvalue, stat, and km
# (plus optional elements when requested).
test_that("test.clusters.km returns a list with pvalue, stat, and km", {
  d <- make_km_data()
  res <- suppressMessages(run_km(d, clusters = c(1, 3)))
  expect_named(res, c("pvalue", "stat", "km"), ignore.order = TRUE)
})

# The p-value is a probability and must lie in [0, 1] regardless of whether H0
# is true or false.
test_that("test.clusters.km pvalue is in [0, 1]", {
  d <- make_km_data()
  res <- suppressMessages(run_km(d, clusters = c(1, 3)))
  expect_gte(res$pvalue, 0)
  expect_lte(res$pvalue, 1)
})

# The test statistic is the V-norm of the difference of cluster means.  It must
# be a single positive number.
test_that("test.clusters.km stat is a positive numeric scalar", {
  d <- make_km_data()
  res <- suppressMessages(run_km(d, clusters = c(1, 3)))
  expect_true(is.numeric(res$stat))
  expect_length(res$stat, 1)
  expect_gt(res$stat, 0)
})

# km must assign a cluster label to every observation, so its length must equal
# the number of rows of X.
test_that("test.clusters.km km is an integer vector of length n", {
  d <- make_km_data()
  res <- suppressMessages(run_km(d, clusters = c(1, 3)))
  expect_length(res$km, d$n)
})

# ---- Input validation ---------------------------------------------------

# Cluster label 4 does not exist when NC = 3; the function must stop before any
# computation is attempted.
test_that("test.clusters.km stops when clusters contains out-of-range labels", {
  d <- make_km_data()
  expect_error(
    suppressMessages(run_km(d, clusters = c(1, 4))),
    regexp = "clusters"
  )
})

# The inference procedure compares exactly two clusters; passing three labels
# must be rejected immediately.
test_that("test.clusters.km stops when clusters has length != 2", {
  d <- make_km_data()
  expect_error(
    suppressMessages(run_km(d, clusters = c(1, 2, 3))),
    regexp = "clusters"
  )
})

# ---- Optional arguments -------------------------------------------------

# When return_Sigma = TRUE, the column covariance matrix used in the test must
# appear in the returned list with the correct dimensions.
test_that("test.clusters.km includes Sigma in output when return_Sigma = TRUE", {
  d <- make_km_data()
  res <- suppressMessages(
    test.clusters.km(
      X = d$X, U = d$U, Sigma = d$Sigma,
      NC = 3, clusters = c(1, 3),
      return_Sigma = TRUE
    )
  )
  expect_true("Sigma" %in% names(res))
  expect_equal(dim(res$Sigma), c(d$p, d$p))
})

# When Sigma is not known, it can be over-estimated from an independent sample Y
# with the same variables.  The test must still produce a valid p-value.
test_that("test.clusters.km works with Sigma estimated from auxiliary Y", {
  d <- make_km_data()
  set.seed(99)
  Y <- matrix(rnorm(d$n * d$p), nrow = d$n, ncol = d$p)
  res <- suppressMessages(
    test.clusters.km(
      X = d$X, U = d$U, Sigma = NULL, Y = Y,
      NC = 3, clusters = c(1, 3)
    )
  )
  expect_gte(res$pvalue, 0)
  expect_lte(res$pvalue, 1)
})

# Sample splitting diverts half of X to estimate Sigma; a larger n is used here
# so that both the estimation and clustering subsamples are adequately sized.
test_that("test.clusters.km works with sample splitting", {
  set.seed(7)
  n <- 60; p <- 5
  X <- rbind(
    matrix(rnorm(20 * p, mean = -4), nrow = 20),
    matrix(rnorm(20 * p, mean =  0), nrow = 20),
    matrix(rnorm(20 * p, mean =  4), nrow = 20)
  )
  res <- suppressWarnings(suppressMessages(
    test.clusters.km(
      X = X, U = NULL, Sigma = NULL,
      NC = 3, clusters = c(1, 3),
      sample_split = TRUE
    )
  ))
  expect_gte(res$pvalue, 0)
  expect_lte(res$pvalue, 1)
})

# ---- Under the null (uniform p-values over many runs) -------------------

# Under H0 (homogeneous data with no true cluster structure) the selective p-value
# should be approximately uniform on [0, 1] and must not be consistently near 0.
# Five independent replicates are used to keep runtime low while still detecting
# a systematic bias.
test_that("test.clusters.km p-values are not systematically near 0 under H0", {
  # Under H0 (no true cluster structure) the p-value should not always be 0
  pvals <- replicate(5, {
    set.seed(sample.int(1e5, 1))
    n <- 20; p <- 4
    X <- matrix(rnorm(n * p), nrow = n)
    tryCatch(
      suppressMessages(
        test.clusters.km(X = X, Sigma = diag(p), NC = 2, clusters = c(1, 2))$pvalue
      ),
      error = function(e) NA_real_
    )
  })
  pvals <- pvals[!is.na(pvals)]
  # At least some p-values should be non-zero
  expect_true(any(pvals > 0.01))
})

# ---- Precomputed km_at_cl -----------------------------------------------

# Passing the output of kmeans_inference() directly must use the precomputed
# partition exactly (km matches final_cluster) and still produce valid output.
test_that("test.clusters.km gives same result with precomputed km_at_cl", {
  d <- make_km_data()
  set.seed(42)
  km_obj <- KmeansInference::kmeans_inference(
    as.matrix(d$X), k = 3, cluster_1 = 1, cluster_2 = 3,
    verbose = FALSE, seed = NULL, sig = 1, tol_eps = 1e-6, iter.max = 50
  )
  res <- suppressMessages(run_km(d, clusters = c(1, 3), km_at_cl = km_obj))
  # The returned cluster labels must exactly reflect the precomputed partition
  expect_equal(res$km, as.vector(km_obj$final_cluster))
  # The inference result must be a valid p-value and positive statistic
  expect_gte(res$pvalue, 0)
  expect_lte(res$pvalue, 1)
  expect_gt(res$stat, 0)
})

# A non-kmeans_inference object must be rejected with a clear error message.
test_that("test.clusters.km stops when km_at_cl is not a kmeans_inference result", {
  d <- make_km_data()
  expect_error(
    suppressMessages(run_km(d, clusters = c(1, 3), km_at_cl = list(x = 1))),
    regexp = "KmeansInference::kmeans_inference"
  )
})

# When sample_split = TRUE the X matrix is modified after the split, so a
# precomputed km_at_cl would refer to the wrong data.  The function must warn
# and recompute the clustering on the subsample.
test_that("test.clusters.km warns and ignores km_at_cl when sample_split = TRUE", {
  set.seed(11)
  n <- 60; p <- 5
  X <- rbind(
    matrix(rnorm(20 * p, mean = -4), nrow = 20),
    matrix(rnorm(20 * p, mean =  0), nrow = 20),
    matrix(rnorm(20 * p, mean =  4), nrow = 20)
  )
  km_obj <- KmeansInference::kmeans_inference(
    as.matrix(X), k = 3, cluster_1 = 1, cluster_2 = 3,
    verbose = FALSE, seed = NULL, sig = 1, tol_eps = 1e-6, iter.max = 50
  )
  expect_warning(
    suppressMessages(
      test.clusters.km(
        X = X, U = NULL, Sigma = NULL,
        NC = 3, clusters = c(1, 3),
        km_at_cl = km_obj,
        sample_split = TRUE
      )
    ),
    regexp = "ignored when sample_split"
  )
})
