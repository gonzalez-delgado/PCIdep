# ==============================================================================
# Test suite: ARI
#
# Verifies the behaviour of ARI(), which computes the Adjusted Rand Index
# between two cluster assignment vectors.
#
# The ARI equals 1 for identical partitions (up to label permutation), has an
# expected value of 0 under independence, and equals 0 for the specific
# degenerate case where every observation is in its own singleton cluster vs.
# every observation in one cluster.
#
# Tests cover:
#   - ARI = 1 for identical partitions.
#   - ARI = 1 when labels are permuted but the partition is the same.
#   - ARI near 0 for two independent random partitions (large n).
#   - Correct value for a small hand-computable example.
#   - Error for inputs of different lengths.
# ==============================================================================

# Identical partitions must yield ARI = 1.
test_that("ARI returns 1 for identical partitions", {
  cl <- c(1, 1, 2, 2, 3, 3)
  expect_equal(ARI(cl, cl), 1)
})

# A consistent relabeling is the same partition, so ARI must still be 1.
test_that("ARI returns 1 when labels are permuted consistently", {
  cl1 <- c(1, 1, 2, 2, 3, 3)
  cl2 <- c(3, 3, 1, 1, 2, 2)
  expect_equal(ARI(cl1, cl2), 1)
})

# For two large, independent random partitions the ARI should be close to 0.
test_that("ARI is close to 0 for two independent random partitions", {
  set.seed(42)
  cl1 <- sample(1:5, 500, replace = TRUE)
  cl2 <- sample(1:5, 500, replace = TRUE)
  expect_lt(abs(ARI(cl1, cl2)), 0.05)
})

# Hand-computed reference case.
# cl1 = (1,1,1,2,2,2), cl2 = (1,1,2,2,3,3)
# Contingency table:
#      1  2  3
#   1  2  1  0   row sum: 3
#   2  0  1  2   row sum: 3
# col sums: 2  2  2
# n = 6, C(6,2) = 15
# sum C(n_ij,2): C(2,2)+C(1,2)+C(1,2)+C(2,2) = 1+0+0+1 = 2
# sum C(a_i,2):  C(3,2)+C(3,2) = 3+3 = 6
# sum C(b_j,2):  C(2,2)+C(2,2)+C(2,2) = 1+1+1 = 3
# expected = 6*3/15 = 1.2
# numerator   = 2 - 1.2 = 0.8
# denominator = 0.5*(6+3) - 1.2 = 4.5 - 1.2 = 3.3
# ARI = 0.8/3.3 ≈ 0.2424...
test_that("ARI matches hand-computed reference value", {
  cl1 <- c(1, 1, 1, 2, 2, 2)
  cl2 <- c(1, 1, 2, 2, 3, 3)
  expect_equal(ARI(cl1, cl2), 0.8 / 3.3, tolerance = 1e-10)
})

# Mismatched lengths must raise an error.
test_that("ARI raises an error when inputs have different lengths", {
  expect_error(ARI(c(1, 2, 3), c(1, 2)), "same length")
})
