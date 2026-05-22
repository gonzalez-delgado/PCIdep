# ==============================================================================
# Test suite: setup.model
#
# Verifies the behaviour of setup.model(), which prepares the data matrix X and
# the dependence structures U and Sigma before clustering and testing.  It
# handles four distinct responsibilities that are each tested here:
#
#   1. U handling     — defaults to identity when NULL; validates positive-
#                       definiteness; warns when U is not compound symmetry.
#   2. Sigma handling — preserves a user-supplied Sigma; estimates Sigma from
#                       an auxiliary sample Y; performs sample splitting to
#                       generate Y from X; errors when Sigma cannot be obtained.
#   3. Dimension checks — errors when X, U, and Sigma have incompatible sizes,
#                         or when Y has a different number of columns than X.
#   4. Return structure — the returned list always contains X, U, and Sigma.
# ==============================================================================

set.seed(42)

n <- 10
p <- 4

make_X <- function(n, p) matrix(rnorm(n * p), nrow = n, ncol = p)
make_cs_U <- function(n, rho = 0.3) {
  U <- matrix(rho, nrow = n, ncol = n)
  diag(U) <- 1
  U
}

# ---- U handling ---------------------------------------------------------

# When U is not supplied, row-wise independence is assumed and U is replaced by
# the n x n identity matrix before any further processing.
test_that("setup.model sets U to identity when U is NULL", {
  X <- make_X(n, p)
  Sigma <- diag(p)
  res <- suppressMessages(setup.model(X = X, Sigma = Sigma))
  expect_equal(dim(res$U), c(n, n))
  expect_equal(as.numeric(res$U), as.numeric(diag(n)))
})

# A valid compound-symmetry U must be accepted and stored with correct dimensions.
test_that("setup.model preserves a valid CS matrix U", {
  X <- make_X(n, p)
  U <- make_cs_U(n)
  Sigma <- diag(p)
  res <- suppressMessages(setup.model(X = X, U = U, Sigma = Sigma))
  expect_equal(dim(res$U), c(n, n))
})

# A Toeplitz U is positive-definite but not CS.  The function should proceed but
# warn the user that selective type I error control may be compromised.
test_that("setup.model emits a warning for a non-CS U", {
  X <- make_X(n, p)
  U <- toeplitz(seq(1, 0.1, length.out = n))
  Sigma <- diag(p)
  expect_warning(
    suppressMessages(setup.model(X = X, U = U, Sigma = Sigma)),
    regexp = "Compound Symmetry"
  )
})

# A matrix of all -1s is not positive-definite; the function must stop with an
# informative error rather than silently producing invalid results.
test_that("setup.model stops on a non-positive-definite U", {
  X <- make_X(n, p)
  U <- matrix(-1, nrow = n, ncol = n)   # not PD
  Sigma <- diag(p)
  expect_error(
    suppressMessages(setup.model(X = X, U = U, Sigma = Sigma)),
    regexp = "positive-definite"
  )
})

# ---- Sigma handling -----------------------------------------------------

# A user-supplied positive-definite Sigma (Toeplitz here) must be stored with
# correct dimensions and returned unchanged.
test_that("setup.model preserves a valid known Sigma", {
  X <- make_X(n, p)
  Sigma <- toeplitz(seq(1, 0.1, length.out = p))
  res <- suppressMessages(setup.model(X = X, Sigma = Sigma))
  expect_equal(dim(res$Sigma), c(p, p))
})

# Negating the identity matrix gives a negative-definite matrix.  The function
# must reject it with a clear error.
test_that("setup.model stops on a non-positive-definite Sigma", {
  X <- make_X(n, p)
  Sigma <- -diag(p)
  expect_error(
    suppressMessages(setup.model(X = X, Sigma = Sigma)),
    regexp = "positive-definite"
  )
})

# When Sigma is not supplied but an independent sample Y is provided, Sigma is
# estimated from Y.  The estimate must be a p x p symmetric matrix.
test_that("setup.model estimates Sigma from an auxiliary sample Y", {
  X <- make_X(n, p)
  Y <- make_X(n, p)
  res <- suppressMessages(setup.model(X = X, Y = Y))
  expect_equal(dim(res$Sigma), c(p, p))
  # Estimated Sigma should be symmetric
  expect_equal(as.matrix(res$Sigma), t(as.matrix(res$Sigma)), tolerance = 1e-10)
})

# If neither Sigma nor Y is given and sample_split is FALSE, there is no way to
# obtain Sigma; the function must stop with an informative error.
test_that("setup.model stops when Sigma is NULL and Y is absent without sample splitting", {
  X <- make_X(n, p)
  expect_error(
    suppressMessages(setup.model(X = X, sample_split = FALSE)),
    regexp = "i.i.d. sample Y"
  )
})

# With sample_split = TRUE, half the rows of X are diverted to Y for Sigma
# estimation, so the returned X must have fewer rows than the original.
test_that("setup.model performs sample splitting when sample_split is TRUE", {
  X <- make_X(20, p)
  res <- suppressMessages(setup.model(X = X, sample_split = TRUE))
  # After splitting, X should have fewer rows than 20
  expect_lt(nrow(res$X), 20)
  expect_equal(dim(res$Sigma), c(p, p))
})

# nY controls how many rows go to Y; with nY = 5 out of 20, the clustering
# matrix X must have the remaining 15 rows.
test_that("setup.model respects nY when sample splitting", {
  X <- make_X(20, p)
  res <- suppressMessages(setup.model(X = X, sample_split = TRUE, nY = 5))
  expect_equal(nrow(res$X), 15)
})

# ---- Dimension checks ---------------------------------------------------

# Sigma has p+1 columns but X has only p; the dimension mismatch must be caught.
test_that("setup.model stops on incompatible dimensions between X and Sigma", {
  X <- make_X(n, p)
  Sigma <- diag(p + 1)
  expect_error(
    suppressMessages(setup.model(X = X, Sigma = Sigma)),
    regexp = "compatible dimensions"
  )
})

# Y must share the same variables (columns) as X; a column mismatch must be
# caught before any estimation is attempted.
test_that("setup.model stops when Y and X have different numbers of columns", {
  X <- make_X(n, p)
  Y <- make_X(n, p + 1)
  expect_error(
    suppressMessages(setup.model(X = X, Y = Y)),
    regexp = "same number of variables"
  )
})

# ---- Return structure ---------------------------------------------------

# The return value must always be a named list with exactly the elements X, U,
# and Sigma, regardless of which inputs were supplied.
test_that("setup.model return value contains X, U, and Sigma", {
  X <- make_X(n, p)
  Sigma <- diag(p)
  res <- suppressMessages(setup.model(X = X, Sigma = Sigma))
  expect_named(res, c("X", "U", "Sigma"), ignore.order = TRUE)
})
