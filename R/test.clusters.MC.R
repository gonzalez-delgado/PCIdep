#' Test for the difference of two cluster means after any clustering algorithm, for matrix normal model with arbitrary scale matrices.
#' 
#' @param X A \eqn{n \times p} matrix drawn from a \eqn{n \times p} matrix normal distribution \eqn{\mathcal{MN}(}\code{M}, \code{U}, \code{Sigma}\eqn{)}. \code{X} must have \eqn{n} rows and \eqn{p} columns.
#' @param U A \eqn{n \times n} positive-definite matrix describing the dependence structure between the rows in \code{X}. If \code{NULL}, observations are considered independent and \code{U} is set to the \eqn{n \times n} identity matrix.
#' @param Sigma A \eqn{p \times p} positive-definite matrix describing the dependence structure between the columns in \code{X}. If \code{NULL}, \code{Sigma} is over-estimated (in the sense of the Loewner partial order).
#' @param Y If \code{Sigma} is \code{NULL}, an i.i.d. copy of \code{X} allowing its estimation. \code{Y} must have the same number of columns as \code{X}.
#' @param UY If \code{Sigma} is \code{NULL}, a positive-definite matrix describing the dependence structure between the rows in \code{Y}. If \code{NULL} and its inverse is not provided, set to the identity matrix by default.
#' @param precUY The inverse matrix of \code{UY}, that can be provided to increase computational efficiency. If \code{UY} is not \code{NULL} and \code{precUY} is \code{NULL}, \code{precUY} is obtained by inverting \code{UY}.
#' @param clusters A vector of two integers from 1 to \code{NC} indicating the pair of clusters whose means have to be compared.
#' @param cl_fun A function returning assignments to clusters. The function must take as input the data matrix \code{X} and the number of clusters \code{NC}.
#' @param NC The number of clusters to choose, that will be passed as argument to \code{cl_fun}. Must be set to \code{NULL} if not required by \code{cl_fun}.
#' @param cl The result of clustering \code{X} using \code{cl_fun}. It can useful to precompute this quantity before choosing \code{clusters}.
#' @param ndraws The number of Monte Carlo iterations.
#'
#' @return 
#' \itemize{
#'   \item pvalue - The p-value for the difference of cluster means.
#'   \item stat - The test statistic.
#'   \item stdrr - The Monte Carlo standard error.
#'   \item clusters - The partition of the \code{n} observations retrieved by the clustering algorithm.
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
#' # Using HDBSCAN clustering from dbscan package. This algorithm selects 
#' # automatically the number of clusters NC.
#' # Additional clustering parameters must be set as default values
#' # when defining cl_fun.
#' 
#' # install.packages('dbscan')
#'
#' hdbscan.clustering <- function(X, NC = NULL, min.occupancy = 5){
#'  
#'  X.clus <- dbscan::hdbscan(X, minPts = min.occupancy)
#'  return(X.clus$cluster + 1)
#'  
#' }
#'
#' # We start by clustering the data
#' clusters_X <- hdbscan.clustering(X)
#' # We test for the equality of clusters 3 and 1
#' test.clusters.MC(X, U = U, Sigma = Sigma, clusters = c(3,1),
#'  cl = clusters_X, cl_fun = hdbscan.clustering, NC = NULL, ndraws = 500)
#'
#' @references [1] L. L. Gao, J. Bien, and D. Witten. Selective inference for hierarchical clustering. Journal of the American Statistical Association, 0(0):1–11, 2022.
#' 
#' @export

test.clusters.MC <- function(X, U = NULL, Sigma = NULL, Y = NULL, UY = NULL, precUY = NULL, 
                                 clusters, cl_fun, NC = NULL, cl = NULL, ndraws = 2000, ISK = 1){
  
  #### Initial checks and pre-processing #######################################
  
  # Check U
  if(is.null(U)){
    
    U <- Matrix::Diagonal(dim(X)[1]) # Independent observations by default
    cat('U is not provided: observations are considered independent with unit variance.\n')  
    
  }else{ 
    
    if(!matrixNormal::is.positive.definite(U)){ # Check for positive-definiteness of U
      
      stop('U must be positive-definite.')}else{ 
        
        if(!is.CS(U)){warning('U is not Compound Symmetry: selective type I error control might be lost if the deviation from the CS structure is large.\n')}
        U <- Matrix::Matrix(U) # Memory efficiency
        
      }}

  # Check Sigma
  # If Sigma is not provided, estimate it
  if(is.null(Sigma)){
    
    if(is.null(Y)){
      
      stop('Sigma is not provided. An i.i.d. sample Y must be provided to allow its over-estimation.\n')} # Need to provide i.i.d. copy of Y if Sigma is NULL
    
    if(dim(Y)[2] != dim(X)[2]){
      
      stop('Y and X must have the number of variables.\n')} # X and Y must have the same number of features
    
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
      
      stop('Sigma must be positive-definite.\n')}else{
        
        Sigma <- Matrix::Matrix(Sigma) # Memory efficiency
      }
  }
  
  #### Cluster data ############################################################
  
  if(is.null(cl)){cl <- cl_fun(X, NC)}
  NC <- length(unique(cl))

  # Check for correct clustering specification
  if(!all(clusters %in% c(1:NC)) | length(clusters) != 2){stop('clusters must be a vector of two integers between 1 and NC.\n')}
  
  #### Test for the difference of cluster means ################################
  
  # Select individuals from each cluster
  n1 <- sum(cl == clusters[1])
  n2 <- sum(cl == clusters[2])
  nu <- Matrix::Matrix((cl == clusters[1])/n1 - (cl == clusters[2])/n2)
  norm2_nu <- 1/n1 + 1/n2
  
  # Difference of cluster means 
  diff_means <- Matrix::Matrix(Matrix::t(nu)%*%X) 
  
  # Computation of the norm \norm{x}_V = \sqrt{x^T V^{-1} x},
  # where V = nu^T U nu Sigma
  norm2U_nu <- Matrix::crossprod(nu, U)%*%nu
  V_g1g2 <- norm2U_nu[1]*Sigma
  V_g1g2_inv <- Matrix::solve(V_g1g2)
  stat_V <- as.numeric(sqrt(diff_means %*% Matrix::tcrossprod(V_g1g2_inv,diff_means))) # Test statistic (norm_V{diff_means})
  
  # Monte-Carlo approximation of the p-value, without explicit computation of the truncation set.
  # Code adapted from clusterval package (Gao et al. 2022)
     
  prop_k2 <- n2/(n1+n2)
  log_survives <- rep(NA, ndraws)
  phi <- stats::rnorm(ndraws, stat_V, ISK*sqrt(sum(diag(Sigma)))) # N(stat, Tr(Sigma)^0.5)
      
  diff_means <- as.numeric(diff_means)
  k1_constant <- prop_k2*exp(log(abs(diff_means)) - log(stat_V))*sign(diff_means)
  k2_constant <- (prop_k2 - 1)*exp(log(abs(diff_means)) - log(stat_V))*sign(diff_means)
  orig_k1 <- t(X[cl == clusters[1], ])
  orig_k2 <- t(X[cl == clusters[2], ])
      
  Xphi <- X
  
  log_survives <- unlist(future.apply::future_lapply(X = 1:ndraws, FUN = function(j) {
    
    if(phi[j] < 0) return(NA)
    
    # Compute perturbed data set for positive phi's
    Xphi <- X
    phi_minus_stat <- phi[j] - stat_V 
    Xphi[cl == clusters[1], ] <- t(orig_k1 + sign(k1_constant)*sign(phi_minus_stat)*exp(log(abs(k1_constant)) + log(abs(phi_minus_stat))))
    Xphi[cl == clusters[2], ] <- t(orig_k2 + sign(k2_constant)*sign(phi_minus_stat)*exp(log(abs(k2_constant)) + log(abs(phi_minus_stat))))
    
    # Recluster the perturbed data set
    cl_Xphi <- cl_fun(Xphi, NC)

    if(preserve.cl(cl, cl_Xphi, clusters)) {
      log_survives <- -(phi[j])^2/2 + (dim(Sigma)[1]-1)*log(phi[j]) - (dim(Sigma)[1]/2 - 1)*log(2) - lgamma(dim(Sigma)[1]/2) -
        stats::dnorm(phi[j], mean=stat_V, sd=ISK*sqrt(sum(diag(Sigma))), log=TRUE)
      return(log_survives)
    }
    
    return(NA)
    
  }, future.seed=TRUE))
  
  # Trim down to only survives
  phi <- phi[!is.na(log_survives)]
  log_survives <- log_survives[!is.na(log_survives)]
  
  survives <- length(log_survives)
  
  # Return nothing if nothing survives
  if(survives == 0) {
    warning("Oops - we didn't generate any samples that preserved the clusters! Try re-running with a larger value of ndraws.")
    return(list(stat = stat_V, pval=NA, stderr=NA, clusters=cl))
  }
      
  #  Approximate p-values
  log_survives_shift <- log_survives - max(log_survives)
  props <- exp(log_survives_shift)/sum(exp(log_survives_shift))
  pv <- sum(props[phi >= stat_V])
  var_pv <- (1 - pv)^2*sum(props[phi >= stat_V]^2) + pv^2*sum(props[phi < stat_V]^2)
      
  return(list(pvalue = pv, stat = stat_V, stderr = sqrt(var_pv),  clusters=cl))
      
}
