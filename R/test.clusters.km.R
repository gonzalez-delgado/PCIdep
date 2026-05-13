#' @title Post-clustering inference after k-means clustering
#' @family post-clustering inference functions
#' @description
#' Performs post-clustering inference for the difference between the means of two clusters obtained by k-means clustering, under a general matrix normal
#' model.
#' 
#' The covariance structure between features \code{Sigma} is estimated if not provided, respecting the selective type I error control. The covariance
#' matrix between observations \code{U} is assumed known and should have a compound symmetry structure (see details).
#'
#' @param X A numeric \eqn{n \times p} matrix assumed to arise from a matrix normal distribution \eqn{\mathcal{MN}(M, U, \Sigma)}.
#' @param U An \eqn{n \times n} positive-definite matrix describing the dependence structure between the rows of \code{X}. If \code{NULL},
#'   observations are assumed to be independent and \code{U} is set to the \eqn{n \times n} identity matrix.
#' @param Sigma A \eqn{p \times p} positive-definite matrix describing the dependence structure between the columns of \code{X}. If \code{NULL},
#'   \code{Sigma} is over-estimated from an auxiliary independent sample \code{Y} (in the sense of the Loewner partial order).
#' @param Y If \code{Sigma} is \code{NULL}, an independent copy of \code{X} used to estimate \code{Sigma}. It must have the same number of columns as
#'   \code{X}.
#' @param UY If \code{Sigma} is \code{NULL}, an \eqn{n_Y \times n_Y} positive-definite matrix describing the dependence structure between the
#'   rows of \code{Y}. If \code{NULL} and \code{precUY} is not provided, the identity matrix is used by default.
#' @param precUY The inverse of \code{UY}. Providing \code{precUY} may improve computational efficiency. If \code{UY} is provided but \code{precUY} is
#'   \code{NULL}, it is computed internally.
#' @param NC Integer scalar giving the number of clusters.
#' @param clusters Integer vector of length 2 specifying the pair of clusters to compare. Entries must belong to \code{1:NC}.
#' @param itermax Integer. Maximum number of iterations for the k-means algorithm, passed to \code{KmeansInference::kmeans_inference}.
#' @param tol Numeric tolerance parameter controlling convergence of the k-means algorithm, passed as \code{tol_eps}.
#' @param sample_split Logical. Whether to use sample splitting to estimate \code{Sigma} when \code{Sigma = NULL}. Ignored when \code{Sigma} is provided by the user.
#' @param nY Integer. If \code{Y} is not provided and \code{sample_split = TRUE}, the number of rows of the auxiliary sample \code{Y} used to estimate \code{Sigma}. If \code{nY} is \code{NULL}, half of the rows of \code{X} are used for estimation. Ignored when \code{Sigma} is provided by the user.
#' 
#' @details
#' Selective type I error control is guaranteed when the row covariance matrix
#' \eqn{\mathbf{U}} has a compound symmetry (CS) structure, i.e.,
#' \deqn{
#'   \mathbf{U} = (a - b)\mathbf{I}_n + b\mathbf{1}_n,
#' }
#' for some \eqn{a/(n - 1) < b < a}. In practice, the method is robust to moderate deviations from the CS structure. In particular, reliable 
#' performance is typically observed when \eqn{\mathbf{U}} remains close to CS, for example in autoregressive (AR(1)), diagonal, banded, or Toeplitz 
#' covariance structures, depending on their parameters. We therefore recommend applying the method primarily in settings where \eqn{\mathbf{U}} is known and
#' does not deviate substantially from the compound symmetry structure.
#'  
#' When \code{Sigma = NULL}, the column covariance matrix is estimated from an auxiliary independent sample \code{Y}. When \code{U = NULL}, row-wise
#' independence is assumed.
#' 
#' The clustering step is performed using the function \code{KmeansInference::kmeans_inference}, which provides both the clustering
#' partition and the truncation region associated with the selection event. The p-value is computed using an exact characterization of the truncation
#' region for the Euclidean norm, following Chen and Witten (2022), and then mapped to the corresponding statistic under the general matrix normal model.
#'
#' @return
#' A named list with the following components:
#' \describe{
#'   \item{pvalue}{The p-value for testing equality of the two selected cluster
#'   means.}
#'   \item{stat}{The observed test statistic.}
#'   \item{km}{An integer vector of length \eqn{n} giving the cluster
#'   membership of each observation returned by the k-means algorithm.}
#'   \item{Sigma}{The column covariance matrix used in the test, either provided
#'   by the user or estimated from \code{Y}.}
#' }
#'
#' @examples
#' n <- 50
#' p <- 20
#'
#' # Simulatin under the null hypothesis
#' M <- Matrix::Matrix(0, nrow = n, ncol = p)
#' Sigma <- stats::toeplitz(seq(1, 0.1, length.out = p))
#' U <- matrixNormal::I(n)
#'
#' X <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
#' Y <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
#'
#' # k-means with known Sigma
#' test.km <- test.clusters.km(
#'   X = X, U = U, Sigma = Sigma,
#'   NC = 3, clusters = sample(1:3, 2)
#' )
#' test.km$pvalue
#'
#' # k-means with over-estimation of Sigma
#' test.km <- test.clusters.km(
#'   X = X, U = U, Sigma = NULL, Y = Y,
#'   NC = 3, clusters = sample(1:3, 2)
#' )
#' test.km$pvalue
#'
#' @references
#' Gao, L. L., Bien, J., and Witten, D. (2022).
#' Selective inference for hierarchical clustering.
#' \emph{Journal of the American Statistical Association}, 117(540), 2533--2547.
#'
#' Chen, Y. T., and Witten, D. M. (2023).
#' Selective inference for k-means clustering.
#' \emph{Journal of Machine Learning Research}, 24(152), 1-41.
#' 
#' González-Delgado, J., Deronzier, M. Cortés, J., and Neuvial, P. (2023)
#' Post-clustering Inference under Dependence. 
#' \emph{arXiv.2310.11822}.
#'
#' @export

test.clusters.km <- function(X, U = NULL, Sigma = NULL, Y = NULL, UY = NULL, precUY = NULL, NC, clusters, itermax = 10, tol = 1e-6, sample_split = FALSE, nY = NULL){
  
  # --------------- Initial checks and pre-processing ---------------

  # Check for correct clustering specification
  if(!all(clusters %in% c(1:NC)) | length(clusters) != 2){stop('clusters must be a vector of two integers between 1 and NC.')}
  
  # Set up data and dependency structures, estimate Sigma if needed
  setup_model <- setup.model(X = X, U = U, Sigma = Sigma, Y = Y, UY = UY, precUY = precUY, sample_split = sample_split, nY = nY)
  X <- setup_model$X
  U <- setup_model$U
  Sigma <- setup_model$Sigma
 
  # --------------- Cluster data ---------------
  
  # K-means clustering from KmeansInference package 
  test_kmeans <- KmeansInference::kmeans_inference(as.matrix(X), k = NC, cluster_1 = clusters[1], cluster_2 = clusters[2], verbose = FALSE, seed = NULL, sig = 1, tol_eps = tol, iter.max = itermax)
  
  # --------------- Test for the difference of cluster means ---------------
  
  # Select individuals from each cluster
  km_at_cl <- as.vector(test_kmeans$final_cluster) # Clustering partition
  n1 <- sum(km_at_cl == clusters[1])
  n2 <- sum(km_at_cl == clusters[2])
  nu <- Matrix::Matrix((km_at_cl == clusters[1])/n1 - (km_at_cl == clusters[2])/n2)
  norm2_nu <- 1/n1 + 1/n2
  
  # Difference of cluster means 
  diff_means <- Matrix::Matrix(Matrix::t(nu)%*%X) 
  
  # Computation of the norm \norm{x}_V = \sqrt{x^T V^{-1} x}, where V = nu^T U nu Sigma
  norm2U_nu <- Matrix::crossprod(nu, U)%*%nu
  V_g1g2 <- norm2U_nu[1]*Sigma
  V_g1g2_inv <- Matrix::solve(V_g1g2)
  stat_V <- as.numeric(sqrt(diff_means %*% Matrix::tcrossprod(V_g1g2_inv,diff_means))) # Test statistic (norm_V{diff_means})
  
  # Computation of the truncation set for the 2-norm in Chen et al. 2022
  S2 <- test_kmeans$final_interval # Truncation set for the 2-norm
  stat_2 <- test_kmeans$test_stat # Test statistic for the 2-norm (2-norm of diff_means)
  SV <- S2*as.numeric(exp(log(stat_V)-log(stat_2))) # Truncation set for norm_V.
  
  # Computation of the p-value using the chisq distribution
  I1 <- intervals::Intervals(c(stat_V^2, Inf))
  I2 <- suppressWarnings(intervals::interval_intersection(I1, SV^2))
  pv <- clusterpval::TChisqRatioApprox(dim(Sigma)[1], I2, SV^2)
  
  return(list(pvalue = pv, stat = stat_V, km = km_at_cl, Sigma = Sigma))}


