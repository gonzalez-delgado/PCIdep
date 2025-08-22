#' Test for the difference of two cluster means after hierarchical clustering, for matrix normal model with arbitrary scale matrices.
#' Supported linkages (as in Gao et al. 2022)  are "single", "average", "centroid", "ward.D", "median", "mcquitty" and "complete".
#' 
#' @param X A \eqn{n \times p} matrix drawn from a \eqn{n \times p} matrix normal distribution \eqn{\mathcal{MN}(}\code{M}, \code{U}, \code{Sigma}\eqn{)}. \code{X} must have \eqn{n} rows and \eqn{p} columns.
#' @param U A \eqn{n \times n} positive-definite matrix describing the dependence structure between the rows in \code{X}. If \code{NULL}, observations are considered independent and \code{U} is set to the \eqn{n \times n} identity matrix.
#' @param Sigma A \eqn{p \times p} positive-definite matrix describing the dependence structure between the columns in \code{X}. If \code{NULL}, \code{Sigma} is over-estimated (in the sense of the Loewner partial order).
#' @param Y If \code{Sigma} is \code{NULL}, an i.i.d. copy of \code{X} allowing its estimation. \code{Y} must have the same number of columns as \code{X}.
#' @param UY If \code{Sigma} is \code{NULL}, a positive-definite matrix describing the dependence structure between the rows in \code{Y}. If \code{NULL} and its inverse is not provided, set to the identity matrix by default.
#' @param precUY The inverse matrix of \code{UY}, that can be provided to increase computational efficiency. If \code{UY} is not \code{NULL} and \code{precUY} is \code{NULL}, \code{precUY} is obtained by inverting \code{UY}.
#' @param NC The number of clusters to choose.
#' @param clusters A vector of two integers from \eqn{1} to \code{NC} indicating the pair of clusters whose means have to be compared.
#' @param linkage The type of linkage for hierarchical clustering. Must be either \code{single}, \code{average}, \code{centroid}, \code{ward.D}, \code{median}, \code{mcquitty} or \code{complete}.
#' @param ndraws If linkage is \code{complete}, the number of Monte Carlo iterations.
#'
#' @return 
#' \itemize{
#'   \item pvalue - The p-value for the difference of cluster means.
#'   \item stat - The test statistic.
#'   \item stdrr - If linkage is \code{complete}, the Monte Carlo standard error.
#'   \item hcl - The partition of the \code{n} observations retrieved by the clustering algorithm.
#' }
#'
#' @examples
#' 
#' n <- 50
#' p <- 20
#' M <- Matrix::Matrix(0, nrow = n , ncol = p) # Mean matrix
#' Sigma <- stats::toeplitz(seq(1, 0.1, length = p)) # Sigma: dependence between features
#' U <- matrixNormal::I(n) # U: dependence between observations
#' X <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
#' Y <- matrixNormal::rmatnorm(s = 1, M, U, Sigma) # i.i.d. copy of X
#'
#' # HAC with average linkage under the global null hypothesis
#' test.clusters.hc(X, U, Sigma, NC = 3, clusters = sample(1:3, 2), linkage = "average")
#' # HAC with complete linkage under the global null hypothesis and over-estimation of Sigma
#' test.clusters.hc(X, U, Sigma = NULL, Y = Y, NC = 3, clusters = sample(1:3, 2), linkage = "complete")
#'
#' @references [1] L. L. Gao, J. Bien, and D. Witten. Selective inference for hierarchical clustering. Journal of the American Statistical Association, 0(0):1–11, 2022.
#' 
#' @export


test.clusters.hc <- function(X, U = NULL, Sigma = NULL, Y = NULL, UY = NULL, precUY = NULL, NC, clusters, linkage = 'average', ndraws = 2000){
  
  #### Initial checks and pre-processing #######################################
  
  # Check for correct input linkage
  if(!linkage%in%c("single", "average", "centroid", "ward.D", "median", "mcquitty", "complete")){
    stop('linkage must be one in "single", "average", "centroid", "ward.D", "median", "mcquitty", "complete".')
  }
  
  # Check U
  if(is.null(U)){
    
    U <- Matrix::Diagonal(dim(X)[1]) # Independent observations by default
    cat('U is not provided: observations are considered independent with unit variance.')  
  
  }else{ 
    
    if(!matrixNormal::is.positive.definite(U)){ # Check for positive-definiteness of U
      
      stop('U must be positive-definite.')}else{ 
      
      if(!is.CS(U)){warning('U is not Compound Symmetry: selective type I error control might be lost if the deviation from the CS structure is large.')}
      U <- Matrix::Matrix(U) # Memory efficiency
      
    }}
  
  # Check for correct clustering specification
  if(!all(clusters %in% c(1:NC)) | length(clusters) != 2){stop('clusters must be a vector of two integers between 1 and NC.')}
  
  # Check Sigma
  # If Sigma is not provided, estimate it
  if(is.null(Sigma)){
    
    if(is.null(Y)){
      
      stop('Sigma is not provided. An i.i.d. sample Y must be provided to allow its over-estimation.')} # Need to provide i.i.d. copy of Y if Sigma is NULL
    
    if(dim(Y)[2] != dim(X)[2]){
      
      stop('Y and X must have the number of variables')} # X and Y must have the same number of features
    
    cat('Sigma not provided: plugging an over-estimate.\n')
    
    if(is.null(UY) & is.null(precUY)){
      
      precUY <- Matrix::Diagonal(dim(Y)[1])} # If the matrix U for Y and its inverse are not provided, independent observations by default
    
    if(is.null(precUY) & !is.null(UY)){ # Provide U for Y but not its inverse
      
      UY <- Matrix::Matrix(UY) 
      precUY <- Matrix::solve(UY)}
    
    if(!is.null(precUY)){
      
      precUY <- Matrix::Matrix(precUY)} # Provide the inverse of U for Y
    
    # Estimate Sigma
    Y <- Matrix::Matrix(Y) # Memory efficiency
    Ybar <- Matrix::colMeans(Y)
    Sigma <- Matrix::crossprod(Matrix::t(1/(nrow(Y) - 1)*((Matrix::t(Y) - Ybar))),  Matrix::tcrossprod(precUY, Matrix::t(Y) - Ybar)) # Sigma estimate
    
  }else{# Sigma is known and provided by the user
    
    if(!matrixNormal::is.positive.definite(Sigma)){ # Check for positive-definiteness for Sigma
      
      stop('Sigma must be positive-definite.')}else{
        
        Sigma <- Matrix::Matrix(Sigma) # Memory efficiency
    }
  }
  
  #### Cluster data ############################################################
  
  # Hierarchical clustering
  dismat <- stats::dist(X, method = "euclidean")^2
  hcl <- fastcluster::hclust(dismat, method = linkage) 

  #### Test for the difference of cluster means ################################
  
  # Select individuals from each cluster
  hcl_at_K <- stats::cutree(hcl, NC)  
  n1 <- sum(hcl_at_K == clusters[1])
  n2 <- sum(hcl_at_K == clusters[2])
  nu <- Matrix::Matrix((hcl_at_K == clusters[1])/n1 - (hcl_at_K == clusters[2])/n2)
  norm2_nu <- 1/n1 + 1/n2
  
  # Difference of cluster means 
  diff_means <- Matrix::Matrix(Matrix::t(nu)%*%X) 
  
  # Computation of the norm \norm{x}_V = \sqrt{x^T V^{-1} x},
  # where V = nu^T U nu Sigma
  norm2U_nu <- Matrix::crossprod(nu, U)%*%nu
  V_g1g2 <- norm2U_nu[1]*Sigma
  V_g1g2_inv <- Matrix::solve(V_g1g2)
  stat_V <- as.numeric(sqrt(diff_means %*% Matrix::tcrossprod(V_g1g2_inv,diff_means))) # Test statistic (norm_V{diff_means})
  
  if(!linkage == "complete"){ # Code adapted from clusterpval package (Gao et al. 2022)
    
    # Computation of the truncation set for the 2-norm in Gao et al. 2022
    test_gao <- clusterpval::test_hier_clusters_exact(as.matrix(X), link = linkage, K = NC, k1 = clusters[1], k2 = clusters[2], hcl = hcl, dist = dismat)
    S2_pos <- test_gao$trunc # Truncation set for the 2-norm and positive phi's
    stat_2 <- test_gao$stat # Test statistic for the 2-norm (2-norm of diff_means)
    
    NX <- X - 2*nu%*%diff_means/norm2_nu # Transformed dataset to perturb for negative phi's
    test_gao_NX <- clusterpval::test_hier_clusters_exact(as.matrix(NX), link = linkage, K = NC, k1 = clusters[1], k2 = clusters[2], hcl = hcl, dist = stats::dist(as.matrix(NX), method = "euclidean")^2)
    S2_neg <- test_gao_NX$trunc # Truncation set for the 2-norm and negative phi's
    
    S2 <- intervals::interval_union(S2_pos, S2_neg) # Truncation set for all phi's and the 2-norm
    SV <- S2*as.numeric(exp(log(stat_V)-log(stat_2))) # Truncation set for the V-norm
    
    # Computation of the p-value using the chisq distribution
    I1 <- intervals::Intervals(c(stat_V^2, Inf))
    I2 <- suppressWarnings(intervals::interval_intersection(I1, SV^2))
    pv <- clusterpval::TChisqRatioApprox(dim(Sigma)[1], I2, SV^2)
    
    return(list(pvalue = pv, stat = stat_V, hcl = hcl_at_K, Sigma = Sigma))}else{ # Complete linkage
      
      cat('Clustering with complete linkage. Monte-Carlo approximation of the p-value.\n')
      
      # Monte-Carlo approximation of the p-value for complete linkage, without explicit computation of the truncation set.
      # Code adapted from clusterval package (Gao et al. 2022)
      
      prop_k2 <- n2/(n1+n2)
      log_survives <- rep(NA, ndraws)
      phi <- stats::rnorm(ndraws) + stat_V # N(stat, 1)
      
      diff_means <- as.numeric(diff_means)
      k1_constant <- prop_k2*exp(log(abs(diff_means)) - log(stat_V))*sign(diff_means)
      k2_constant <- (prop_k2 - 1)*exp(log(abs(diff_means)) - log(stat_V))*sign(diff_means)
      orig_k1 <- t(X[hcl_at_K == clusters[1], ])
      orig_k2 <- t(X[hcl_at_K == clusters[2], ])
      
      Xphi <- X

      for(j in 1:ndraws) {
        if(phi[j] < 0) next
        
        # Compute perturbed data set for positive phi's
        Xphi <- X
        phi_minus_stat <- phi[j] - stat_V 
        Xphi[hcl_at_K == clusters[1], ] <- t(orig_k1 + sign(k1_constant)*sign(phi_minus_stat)*exp(log(abs(k1_constant)) + log(abs(phi_minus_stat))))
        Xphi[hcl_at_K == clusters[2], ] <- t(orig_k2 + sign(k2_constant)*sign(phi_minus_stat)*exp(log(abs(k2_constant)) + log(abs(phi_minus_stat))))
        
        # Recluster the perturbed data sets
        hcl_Xphi <- fastcluster::hclust(stats::dist(Xphi)^2, method = "complete")
        clusters_Xphi <- stats::cutree(hcl_Xphi, NC)
        
        if((sum(table(hcl_at_K, clusters_Xphi) != 0) == NC)) { # Check for same cluster
          log_survives[j] <- -phi[j]^2/2 + (dim(Sigma)[1]-1)*log(phi[j]) - (dim(Sigma)[1]/2 - 1)*log(2) - lgamma(dim(Sigma)[1]/2) -
            stats::dnorm(phi[j], mean = stat_V, sd = 1, log = TRUE)
        }
      }
      
      # Trim down to only survives
      phi <- phi[!is.na(log_survives)]
      log_survives <- log_survives[!is.na(log_survives)]
      
      survives <- length(log_survives)
      
      # Return nothing if nothing survives
      if(survives == 0) {
        warning("Oops - we didn't generate any samples that preserved the clusters! Try re-running with a larger value of ndraws.")
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
