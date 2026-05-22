# ==============================================================================
# Test suite: is.CS
#
# Verifies the behaviour of is.CS(), which checks whether a square matrix has
# compound symmetry (CS) structure, i.e. all diagonal entries are equal and all
# off-diagonal entries are equal.  CS is required for guaranteed selective type
# I error control in the post-clustering tests.
#
# Tests cover:
#   - Positive detection: exact CS matrices (uniform off-diagonal, identity).
#   - Negative detection: Toeplitz, non-square, and partially asymmetric matrices.
# ==============================================================================

# A 4x4 matrix with diagonal 1 and off-diagonal 0.3 is the canonical CS example.
test_that("is.CS returns TRUE for a compound symmetry matrix", {
  U <- matrix(0.3, nrow = 4, ncol = 4)
  diag(U) <- 1
  expect_true(is.CS(U))
})

# The identity matrix is CS: all diagonal entries equal 1, all off-diagonal equal 0.
test_that("is.CS returns TRUE for an identity matrix", {
  expect_true(is.CS(diag(5)))
})

# A Toeplitz matrix has decreasing off-diagonal entries, so it is not CS.
test_that("is.CS returns FALSE for a Toeplitz matrix", {
  U <- toeplitz(c(1, 0.5, 0.25, 0.1))
  expect_false(is.CS(U))
})

# is.CS requires a square matrix; a rectangular matrix must return FALSE.
test_that("is.CS returns FALSE for a non-square matrix", {
  U <- matrix(1, nrow = 3, ncol = 4)
  expect_false(is.CS(U))
})

# Even a single differing off-diagonal entry breaks CS structure.
test_that("is.CS returns FALSE when off-diagonal entries are not all equal", {
  U <- diag(4)
  U[1, 2] <- 0.5
  U[2, 1] <- 0.5
  expect_false(is.CS(U))
})
