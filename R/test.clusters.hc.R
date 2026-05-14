#' @title Post-clustering inference after hierarchical clustering
#' @family post-clustering inference functions
#' @description
#' Performs post-clustering inference for the difference between the means of two clusters obtained by hierarchical agglomerative clustering, under a
#' general matrix normal model.
#'
#' The covariance structure between features \code{Sigma} is estimated if not provided, respecting the selective type I error control. The covariance
#' matrix between observations \code{U} is assumed known and should have a compound symmetry structure (see details).
#'
#' Supported linkage criteria (as in Gao et al. 2022) are \code{"single"}, \code{"average"}, \code{"centroid"},
#' \code{"ward.D"}, \code{"median"}, \code{"mcquitty"}, and \code{"complete"}.
#'
#' @param X An \eqn{n \times p} data matrix assumed to arise from a matrix normal distribution
#'   \eqn{\mathcal{MN}(M, U, \Sigma)}.
#' @param U An \eqn{n \times n} positive-definite matrix describing the dependence structure between the rows of \code{X}. If \code{NULL},
#'   observations are assumed to be independent and \code{U} is set to the \eqn{n \times n} identity matrix.
#' @param Sigma A \eqn{p \times p} positive-definite matrix describing the dependence structure between the columns of \code{X}. If \code{NULL},
#'   \code{Sigma} is over-estimated from an auxiliary independent sample \code{Y} (in the sense of the Loewner partial order).
#' @param Y If \code{Sigma} is \code{NULL}, an independent copy of \code{X} used to estimate \code{Sigma}. It must have the same number of columns as
#'   \code{X}.
#' @param UY If \code{Sigma} is \code{NULL}, an \eqn{n_Y \times n_Y} positive-definite matrix describing the dependence structure between the
#'   rows of \code{Y}. If \code{NULL} and \code{precUY} is not provided, the identity matrix is used by default.
#' @param precUY The inverse of \code{UY}. Supplying \code{precUY} may improve computational efficiency. If \code{UY} is provided but \code{precUY} is
#'   \code{NULL}, \code{precUY} is computed internally by matrix inversion.
#' @param NC Integer. Number of clusters in the partition obtained by cutting the hierarchical clustering tree.
#' @param clusters Integer vector of length 2 containing the labels of the two clusters to compare. Its entries must belong to \code{1:NC}.
#' @param linkage Character string specifying the linkage criterion used in hierarchical clustering. Must be one of \code{"single"},
#'   \code{"average"}, \code{"centroid"}, \code{"ward.D"}, \code{"median"}, \code{"mcquitty"}, or \code{"complete"}.
#' @param ndraws Integer. Number of Monte Carlo samples used to approximate the p-value when \code{linkage = "complete"}. Ignored otherwise.
#' @param sample_split Logical. Whether to use sample splitting to estimate \code{Sigma} when \code{Sigma = NULL}. Ignored when \code{Sigma} is provided by the user.
#' @param nY Integer. If \code{Y} is not provided and \code{sample_split = TRUE}, the number of rows of the auxiliary sample \code{Y} used to estimate \code{Sigma}. If \code{nY} is \code{NULL}, half of the rows of \code{X} are used for estimation. Ignored when \code{Sigma} is provided by the user.
#' 
#' @details
#' Selective type I error control is guaranteed when the row covariance matrix \eqn{\mathbf{U}} has a compound symmetry (CS) structure, i.e.,
#' \deqn{
#'   \mathbf{U} = (a - b)\mathbf{I}_n + b\mathbf{1}_n,
#' }
#' for some \eqn{a/(n - 1) < b < a}. In practice, the method is robust to moderate deviations from the CS structure. In particular, reliable 
#' performance is typically observed when \eqn{\mathbf{U}} remains close to CS, for example in autoregressive (AR(1)), diagonal, banded, or Toeplitz 
#' covariance structures, depending on their parameters. We therefore recommend applying the method primarily in settings where \eqn{\mathbf{U}} is known and
#' does not deviate substantially from the compound symmetry structure.
#'
#' When \code{Sigma = NULL}, the column covariance matrix is estimated from an auxiliary independent sample \code{Y}. When \code{U = NULL}, row-wise independence is assumed.
#'
#' For linkage criteria other than \code{"complete"}, the p-value is computed using an exact characterization of the truncation region, adapted from
#' \pkg{clusterpval}. For \code{"complete"} linkage, the truncation region is not computed explicitly and the p-value is approximated by Monte Carlo. In
#' this case, the Monte Carlo procedure can be automatically parallelized using the \pkg{future} package (see examples).
#'
#' @return
#' A named list with the following components:
#' \describe{
#'   \item{pvalue}{The p-value for testing equality of the two selected cluster
#'   means.}
#'   \item{stat}{The observed test statistic.}
#'   \item{stderr}{Monte Carlo standard error of the estimated p-value when
#'   \code{linkage = "complete"}. This component is omitted otherwise.}
#'   \item{hcl}{An integer vector of length \eqn{n} giving the cluster
#'   membership of each observation in the partition with \code{NC} clusters.}
#'   \item{Sigma}{The column covariance matrix used in the test, either provided
#'   by the user or estimated from \code{Y}.}
#' }
#'
#' @examples
#' n <- 50
#' p <- 20
#' # Simulating under the null hypothesis
#' M <- Matrix::Matrix(0, nrow = n, ncol = p) 
#' Sigma <- stats::toeplitz(seq(1, 0.1, length.out = p))
#' U <- matrixNormal::I(n)
#'
#' X <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
#' Y <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
#'
#' # Hierarchical clustering with average linkage and known Sigma
#' test.hc <- test.clusters.hc(
#'   X = X, U = U, Sigma = Sigma,
#'   NC = 3, clusters = sample(1:3, 2),
#'   linkage = "average"
#' )
#' test.hc$pvalue
#'
#' # Hierarchical clustering with complete linkage and estimated Sigma
#' test.hc <- test.clusters.hc(
#'   X = X, U = U, Sigma = NULL, Y = Y,
#'   NC = 3, clusters = sample(1:3, 2),
#'   linkage = "complete"
#' )
#' test.hc$pvalue
#'
#' # Hierarchical clustering with complete linkage, estimated Sigma, and parallelization
#' library(future)
#' plan(multisession, workers = 4) # Set parallelization plan (adjust workers according to your machine)
#' test.hc <- test.clusters.hc(
#'   X = X, U = U, Sigma = NULL, Y = Y,
#'   NC = 3, clusters = sample(1:3, 2),
#'   linkage = "complete", ndraws = 500
#' )
#' test.hc$pvalue
#'
#' @references
#' Gao, L. L., Bien, J., and Witten, D. (2022).
#' Selective inference for hierarchical clustering.
#' \emph{Journal of the American Statistical Association}, 117(540), 2533--2547.
#'
#' González-Delgado, J., Deronzier, M. Cortés, J., and Neuvial, P. (2023)
#' Post-clustering Inference under Dependence. 
#' \emph{arXiv.2310.11822}.
#'
#' @export

test.clusters.hc <- function(X, U = NULL, Sigma = NULL, Y = NULL, UY = NULL, precUY = NULL, NC, clusters, linkage = 'average', ndraws = 2000, sample_split = FALSE, nY = NULL){
  
  # --------------- Initial checks and pre-processing ---------------

  # Check for correct input linkage
  if(!linkage%in%c("single", "average", "centroid", "ward.D", "median", "mcquitty", "complete")){
    stop('linkage must be one in "single", "average", "centroid", "ward.D", "median", "mcquitty", "complete".')}

  # Check for correct clustering specification
  if(!all(clusters %in% c(1:NC)) | length(clusters) != 2){stop('clusters must be a vector of two integers between 1 and NC.')}
  
  # Set up data and dependency structures, estimate Sigma if needed
  setup_model <- setup.model(X = X, U = U, Sigma = Sigma, Y = Y, UY = UY, precUY = precUY, sample_split = sample_split, nY = nY)
  X <- setup_model$X
  U <- setup_model$U
  Sigma <- setup_model$Sigma
  
  # --------------- Cluster data ---------------

  # Hierarchical clustering
  dismat <- stats::dist(X, method = "euclidean")^2
  hcl <- fastcluster::hclust(dismat, method = linkage) 

  # --------------- Test for the difference of cluster means ---------------
  
  # Select individuals from each cluster
  hcl_at_K <- stats::cutree(hcl, NC)  
  n1 <- sum(hcl_at_K == clusters[1])
  n2 <- sum(hcl_at_K == clusters[2])
  nu <- Matrix::Matrix((hcl_at_K == clusters[1])/n1 - (hcl_at_K == clusters[2])/n2)
  norm2_nu <- 1/n1 + 1/n2

  # Difference of cluster means 
  diff_means <- Matrix::Matrix(Matrix::t(nu)%*%X) 
  
  # Computation of the norm \norm{x}_V = \sqrt{x^T V^{-1} x}, where V = nu^T U nu Sigma
  norm2U_nu <- Matrix::crossprod(nu, U)%*%nu
  V_g1g2 <- norm2U_nu[1]*Sigma
  V_g1g2_inv <- Matrix::solve(V_g1g2)
  stat_V <- as.numeric(sqrt(diff_means %*% Matrix::tcrossprod(V_g1g2_inv,diff_means))) # Test statistic (norm_V{diff_means})
  
  if(!linkage == "complete"){ # Code adapted from clusterpval package (Gao et al. 2022)
    
    # Computation of the truncation set for the 2-norm in Gao et al. 2022
    test_gao <- clusterpval::test_hier_clusters_exact(as.matrix(X), link = linkage, K = NC, k1 = clusters[1], k2 = clusters[2], hcl = hcl, dist = dismat)
    S2 <- test_gao$trunc # Truncation set for the 2-norm and positive phi's
    stat_2 <- test_gao$stat # Test statistic for the 2-norm (2-norm of diff_means)
    SV <- S2*as.numeric(exp(log(stat_V)-log(stat_2))) # Truncation set for the V-norm
    
    # Computation of the p-value using the chisq distribution
    I1 <- intervals::Intervals(c(stat_V^2, Inf))
    I2 <- suppressWarnings(intervals::interval_intersection(I1, SV^2))
    pv <- clusterpval::TChisqRatioApprox(dim(Sigma)[1], I2, SV^2)
    
    return(list(pvalue = pv, stat = stat_V, hcl = hcl_at_K, Sigma = Sigma))}else{ # Complete linkage
      
      cat('Clustering with complete linkage. Monte-Carlo approximation of the p-value.\n')
      
      # Monte-Carlo approximation of the p-value for complete linkage, without explicit computation of the truncation set.
      # Code adapted from clusterval package (Gao et al. 2022).
      
      prop_k2 <- n2/(n1+n2)
      log_survives <- rep(NA, ndraws)
      phi <- stats::rnorm(ndraws) + stat_V # N(stat, 1)
      
      diff_means <- as.numeric(diff_means)
      k1_constant <- prop_k2*exp(log(abs(diff_means)) - log(stat_V))*sign(diff_means)
      k2_constant <- (prop_k2 - 1)*exp(log(abs(diff_means)) - log(stat_V))*sign(diff_means)
      orig_k1 <- Matrix::t(X[hcl_at_K == clusters[1], , drop = FALSE])
      orig_k2 <- Matrix::t(X[hcl_at_K == clusters[2], , drop = FALSE])
      
      Xphi <- X

      log_survives <- unlist(future.apply::future_lapply(X = 1:ndraws, FUN = function(j) {

        if(phi[j] < 0){return(NA)}
        
        # Compute perturbed data set for positive phi's
        Xphi <- X
        phi_minus_stat <- phi[j] - stat_V 
        Xphi[hcl_at_K == clusters[1], ] <- Matrix::t(orig_k1 + sign(k1_constant)*sign(phi_minus_stat)*exp(log(abs(k1_constant)) + log(abs(phi_minus_stat))))
        Xphi[hcl_at_K == clusters[2], ] <- Matrix::t(orig_k2 + sign(k2_constant)*sign(phi_minus_stat)*exp(log(abs(k2_constant)) + log(abs(phi_minus_stat))))
        
        # Recluster the perturbed data sets
        hcl_Xphi <- fastcluster::hclust(stats::dist(Xphi)^2, method = "complete")
        clusters_Xphi <- stats::cutree(hcl_Xphi, NC)
        
        if((sum(table(hcl_at_K, clusters_Xphi) != 0) == NC)) { # Check for same cluster
          log_survives[j] <- -phi[j]^2/2 + (dim(Sigma)[1]-1)*log(phi[j]) - (dim(Sigma)[1]/2 - 1)*log(2) - lgamma(dim(Sigma)[1]/2) -
            stats::dnorm(phi[j], mean = stat_V, sd = 1, log = TRUE)
        }

        return(NA)

      }, future.seed=TRUE))
      
      # Trim down to only survives
      phi <- phi[!is.na(log_survives)]
      log_survives <- log_survives[!is.na(log_survives)]
      
      survives <- length(log_survives)
      
      # Return nothing if nothing survives
      if(survives == 0) {
        warning("No samples that preserved the clusters were generated. Try re-running with a larger value of ndraws.")
        return(list(pvalue = NA, stat = NA, stderr = NA))
      }
      
      #  Approximate p-values
      log_survives_shift <- log_survives - max(log_survives)
      props <- exp(log_survives_shift)/sum(exp(log_survives_shift))
      pv <- sum(props[phi >= stat_V])
      var_pv <- (1 - pv)^2*sum(props[phi >= stat_V]^2) + pv^2*sum(props[phi < stat_V]^2)
      
      return(list(pvalue = pv, stat = stat_V, stderr = sqrt(var_pv), hcl = hcl_at_K, Sigma = Sigma))
      
    }
}
