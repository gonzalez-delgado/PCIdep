# ==============================================================================
# Test suite: test.clusters.hc
#
# Verifies the behaviour of test.clusters.hc(), which performs post-clustering
# selective inference for the difference between two cluster means after
# hierarchical agglomerative clustering, under a general matrix normal model.
#
# Two inference strategies are exercised:
#   - Exact truncation set (all linkages except "complete"): the truncation
#     region is characterised analytically following Gao et al. (2022).
#   - Monte Carlo approximation ("complete" linkage): the truncation region is
#     approximated by re-clustering perturbed data sets, yielding an additional
#     stderr output element.
#
# Tests cover:
#   - Basic output structure for average linkage (the default exact path).
#   - All supported exact-computation linkages: single, centroid, ward.D,
#     mcquitty.
#   - Complete linkage (Monte Carlo path): output includes stderr.
#   - Input validation: invalid linkage string, out-of-range or wrongly-sized
#     cluster labels, invalid hcl object.
#   - Optional arguments: precomputed hcl, return_Sigma, auxiliary Y,
#     sample splitting.
# ==============================================================================

set.seed(321)

# make_hc_data() generates a 30 x 5 matrix with three clearly separated groups
# (means -4, 0, +4) together with identity row and column covariances.  The
# strong separation ensures hierarchical clustering reliably recovers three
# distinct groups regardless of the linkage criterion used.
make_hc_data <- function(seed = 1) {
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

# run_hc() is a thin wrapper around test.clusters.hc() that fixes NC = 3 and
# passes data from a make_hc_data() list, reducing boilerplate in each test.
# Extra arguments (including hcl) are forwarded via ....
run_hc <- function(d, clusters = c(1, 3), linkage = "average", ...) {
  test.clusters.hc(
    X = d$X, U = d$U, Sigma = d$Sigma,
    NC = 3, clusters = clusters,
    linkage = linkage,
    ...
  )
}

# ---- Basic output structure for non-complete linkage --------------------

# The return value for exact-computation linkages must contain at least pvalue,
# stat, and hcl (no stderr, which is reserved for the Monte Carlo path).
test_that("test.clusters.hc (average) returns list with pvalue, stat, hcl", {
  d <- make_hc_data()
  res <- suppressMessages(run_hc(d))
  expect_true(all(c("pvalue", "stat", "hcl") %in% names(res)))
})

# The p-value must lie in [0, 1] for the exact-computation path.
test_that("test.clusters.hc (average) pvalue is in [0, 1]", {
  d <- make_hc_data()
  res <- suppressMessages(run_hc(d))
  expect_gte(res$pvalue, 0)
  expect_lte(res$pvalue, 1)
})

# The test statistic is the V-norm of the difference of cluster means and must
# be a single positive number.
test_that("test.clusters.hc (average) stat is a positive numeric scalar", {
  d <- make_hc_data()
  res <- suppressMessages(run_hc(d))
  expect_true(is.numeric(res$stat))
  expect_length(res$stat, 1)
  expect_gt(res$stat, 0)
})

# hcl assigns a cluster label to every observation after cutting the dendrogram;
# it must have exactly n entries.
test_that("test.clusters.hc (average) hcl has length n", {
  d <- make_hc_data()
  res <- suppressMessages(run_hc(d))
  expect_length(res$hcl, d$n)
})

# ---- Other exact-computation linkages -----------------------------------

# All four remaining exact linkages must return a valid p-value in [0, 1].
# Tests are generated programmatically to avoid repetition.
for (lnk in c("single", "centroid", "ward.D", "mcquitty")) {
  local({
    lnk_local <- lnk
    test_that(paste("test.clusters.hc works with linkage =", lnk_local), {
      d <- make_hc_data()
      res <- suppressMessages(run_hc(d, linkage = lnk_local))
      expect_gte(res$pvalue, 0)
      expect_lte(res$pvalue, 1)
    })
  })
}

# ---- Complete linkage (Monte Carlo) -------------------------------------

# Complete linkage uses a Monte Carlo approximation of the truncation region.
# The output must include stderr (the Monte Carlo standard error of the p-value
# estimate) in addition to pvalue, stat, and hcl.  ndraws is kept small to
# limit test runtime while still exercising the full MC path.
test_that("test.clusters.hc (complete) returns stderr in addition to pvalue and stat", {
  d <- make_hc_data()
  res <- suppressMessages(run_hc(d, linkage = "complete", ndraws = 200))
  expect_true(all(c("pvalue", "stat", "stderr", "hcl") %in% names(res)))
  expect_gte(res$pvalue, 0)
  expect_lte(res$pvalue, 1)
  expect_gte(res$stderr, 0)
})

# ---- Input validation ---------------------------------------------------

# An unrecognised linkage string must be rejected before any clustering is run.
test_that("test.clusters.hc stops on invalid linkage", {
  d <- make_hc_data()
  expect_error(
    suppressMessages(run_hc(d, linkage = "invalid")),
    regexp = "linkage"
  )
})

# Cluster label 5 does not exist when NC = 3; the function must stop with an
# error before performing any clustering or testing.
test_that("test.clusters.hc stops when clusters are out of range", {
  d <- make_hc_data()
  expect_error(
    suppressMessages(run_hc(d, clusters = c(1, 5))),
    regexp = "clusters"
  )
})

# The inference procedure compares exactly two clusters; passing three labels
# must be rejected immediately.
test_that("test.clusters.hc stops when clusters has wrong length", {
  d <- make_hc_data()
  expect_error(
    suppressMessages(run_hc(d, clusters = c(1, 2, 3))),
    regexp = "clusters"
  )
})

# ---- Optional arguments -------------------------------------------------

# A precomputed hcl avoids rerunning hclust when testing multiple cluster pairs
# on the same data.  The result must be identical to computing hcl internally.
test_that("test.clusters.hc gives the same result with a precomputed hcl", {
  d <- make_hc_data()
  hcl_obj <- fastcluster::hclust(stats::dist(d$X)^2, method = "average")
  res_internal <- suppressMessages(run_hc(d, clusters = c(1, 3)))
  res_external <- suppressMessages(
    test.clusters.hc(
      X = d$X, U = d$U, Sigma = d$Sigma,
      NC = 3, clusters = c(1, 3),
      hcl = hcl_obj
    )
  )
  expect_equal(res_internal$pvalue, res_external$pvalue)
  expect_equal(res_internal$stat,   res_external$stat)
  expect_equal(res_internal$hcl,    res_external$hcl)
})

# Providing a precomputed hcl should allow testing a second cluster pair
# without any additional clustering cost.
test_that("test.clusters.hc accepts a precomputed hcl for a different cluster pair", {
  d <- make_hc_data()
  hcl_obj <- fastcluster::hclust(stats::dist(d$X)^2, method = "average")
  res <- suppressMessages(
    test.clusters.hc(
      X = d$X, U = d$U, Sigma = d$Sigma,
      NC = 3, clusters = c(1, 2),
      hcl = hcl_obj
    )
  )
  expect_gte(res$pvalue, 0)
  expect_lte(res$pvalue, 1)
})

# The linkage stored in hcl$method must be in the supported set; an unsupported
# method stored inside the hcl object must trigger a clear error.
test_that("test.clusters.hc stops when hcl contains an unsupported linkage method", {
  d <- make_hc_data()
  hcl_obj <- fastcluster::hclust(stats::dist(d$X)^2, method = "average")
  hcl_obj$method <- "unsupported"
  expect_error(
    suppressMessages(run_hc(d, hcl = hcl_obj)),
    regexp = "not supported"
  )
})

# When the user explicitly passes a linkage that disagrees with hcl$method, a
# warning must be emitted and the method stored in hcl must take precedence.
test_that("test.clusters.hc warns and overrides linkage when it mismatches hcl$method", {
  d <- make_hc_data()
  hcl_obj <- fastcluster::hclust(stats::dist(d$X)^2, method = "average")
  # Compute the reference result using the correct linkage (average)
  ref <- suppressMessages(run_hc(d, clusters = c(1, 3), linkage = "average"))
  # Passing linkage = "single" must warn but produce the same result as average
  expect_warning(
    res <- suppressMessages(
      test.clusters.hc(
        X = d$X, U = d$U, Sigma = d$Sigma,
        NC = 3, clusters = c(1, 3),
        linkage = "single",
        hcl = hcl_obj
      )
    ),
    regexp = "does not match"
  )
  expect_equal(res$pvalue, ref$pvalue)
  expect_equal(res$stat,   ref$stat)
})

# Passing a non-hclust object must be caught with a clear error.
test_that("test.clusters.hc stops when hcl is not an hclust object", {
  d <- make_hc_data()
  expect_error(
    suppressMessages(run_hc(d, hcl = list(method = "average"))),
    regexp = "class 'hclust'"
  )
})

# An hcl object built on a different number of observations than nrow(X) is
# incompatible and must be caught with a clear error.
test_that("test.clusters.hc stops when hcl has wrong number of observations", {
  d <- make_hc_data()
  # Build hcl on a dataset with a different number of rows
  X_other <- matrix(rnorm(10 * d$p), nrow = 10)
  hcl_wrong <- fastcluster::hclust(stats::dist(X_other)^2, method = "average")
  expect_error(
    suppressMessages(run_hc(d, hcl = hcl_wrong)),
    regexp = "nrow\\(X\\)|rows"
  )
})

# When sample_split = TRUE, X changes after splitting so a precomputed hcl
# would be stale.  The function must warn and recompute the clustering.
test_that("test.clusters.hc warns and ignores hcl when sample_split = TRUE", {
  set.seed(11)
  n <- 60; p <- 5
  X <- rbind(
    matrix(rnorm(20 * p, mean = -4), nrow = 20),
    matrix(rnorm(20 * p, mean =  0), nrow = 20),
    matrix(rnorm(20 * p, mean =  4), nrow = 20)
  )
  hcl_obj <- fastcluster::hclust(stats::dist(X)^2, method = "average")
  expect_warning(
    suppressMessages(
      withCallingHandlers(
        test.clusters.hc(
          X = X, U = NULL, Sigma = NULL,
          NC = 3, clusters = c(1, 3),
          hcl = hcl_obj,
          sample_split = TRUE
        ),
        warning = function(w) {
          if (grepl("Consider setting", conditionMessage(w)))
            invokeRestart("muffleWarning")
        }
      )
    ),
    regexp = "ignored when sample_split"
  )
})

# When return_Sigma = TRUE, the column covariance used in the test must appear
# in the returned list with the correct p x p dimensions.
test_that("test.clusters.hc includes Sigma when return_Sigma = TRUE", {
  d <- make_hc_data()
  res <- suppressMessages(
    test.clusters.hc(
      X = d$X, U = d$U, Sigma = d$Sigma,
      NC = 3, clusters = c(1, 3),
      return_Sigma = TRUE
    )
  )
  expect_true("Sigma" %in% names(res))
  expect_equal(dim(res$Sigma), c(d$p, d$p))
})

# When Sigma is not known, it can be over-estimated from an independent sample Y.
# The test must still produce a valid p-value.
test_that("test.clusters.hc works with Sigma estimated from auxiliary Y", {
  d <- make_hc_data()
  set.seed(77)
  Y <- matrix(rnorm(d$n * d$p), nrow = d$n)
  res <- suppressMessages(
    test.clusters.hc(
      X = d$X, U = d$U, Sigma = NULL, Y = Y,
      NC = 3, clusters = c(1, 3)
    )
  )
  expect_gte(res$pvalue, 0)
  expect_lte(res$pvalue, 1)
})

# Sample splitting diverts half of X to estimate Sigma; a larger n is used here
# so both the estimation and clustering subsamples are adequately sized.
test_that("test.clusters.hc works with sample splitting", {
  set.seed(11)
  n <- 60; p <- 5
  X <- rbind(
    matrix(rnorm(20 * p, mean = -4), nrow = 20),
    matrix(rnorm(20 * p, mean =  0), nrow = 20),
    matrix(rnorm(20 * p, mean =  4), nrow = 20)
  )
  res <- suppressWarnings(suppressMessages(
    test.clusters.hc(
      X = X, U = NULL, Sigma = NULL,
      NC = 3, clusters = c(1, 3),
      sample_split = TRUE
    )
  ))
  expect_gte(res$pvalue, 0)
  expect_lte(res$pvalue, 1)
})
