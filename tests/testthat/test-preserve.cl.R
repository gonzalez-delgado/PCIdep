# ==============================================================================
# Test suite: preserve.cl
#
# Verifies the behaviour of preserve.cl(), which checks whether the two selected
# clusters in a partition cl remain invariant (up to label permutation) in a
# perturbed partition cl_phi.  This predicate is used in the Monte Carlo loop of
# test.clusters.MC() and the complete-linkage branch of test.clusters.hc() to
# decide whether a perturbed data set falls in the same selection event.
#
# Invariance requires:
#   (i)  All observations in each selected cluster stay together in cl_phi.
#   (ii) No outside observation joins either selected cluster in cl_phi.
#
# Tests cover:
#   - TRUE cases: consistent relabeling, identical partitions.
#   - FALSE cases: cluster split, outside intrusion, two clusters merged.
# ==============================================================================

# All observations keep their relative grouping; only the integer labels change.
test_that("preserve.cl returns TRUE when clusters are relabeled consistently", {
  cl     <- c(1, 1, 2, 2, 3, 3)
  cl_phi <- c(4, 4, 5, 5, 6, 6)
  expect_true(preserve.cl(cl, cl_phi, clusters = c(1, 2)))
})

# An identical partition trivially satisfies the invariance condition.
test_that("preserve.cl returns TRUE when clusters are identical", {
  cl <- c(1, 1, 2, 2, 3, 3)
  expect_true(preserve.cl(cl, cl, clusters = c(1, 2)))
})

# Observation 1 is moved to a different label than observation 2, splitting cluster 1.
test_that("preserve.cl returns FALSE when a selected cluster is split", {
  cl     <- c(1, 1, 2, 2, 3, 3)
  cl_phi <- c(4, 5, 5, 5, 6, 6)   # obs 1 moved away from its cluster
  expect_false(preserve.cl(cl, cl_phi, clusters = c(1, 2)))
})

# Observations 5-6 (originally cluster 3) are assigned to the label of cluster 1,
# violating condition (ii).
test_that("preserve.cl returns FALSE when outside observations join a selected cluster", {
  cl     <- c(1, 1, 2, 2, 3, 3)
  cl_phi <- c(4, 4, 5, 5, 4, 4)   # obs 5-6 (originally cluster 3) join cluster 4
  expect_false(preserve.cl(cl, cl_phi, clusters = c(1, 2)))
})

# Merging the two selected clusters into a single label does not preserve the
# selection event: the original clustering identified clusters 1 and 2 as
# distinct groups to compare, so they must map to different labels in cl_phi.
test_that("preserve.cl returns FALSE when both selected clusters merge into one label", {
  cl     <- c(1, 1, 2, 2, 3, 3)
  cl_phi <- c(4, 4, 4, 4, 5, 5)   # clusters 1 and 2 merged into label 4
  expect_false(preserve.cl(cl, cl_phi, clusters = c(1, 2)))
})
