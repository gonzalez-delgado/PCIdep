#' Test for the difference of two cluster means after k-means clustering, for matrix normal model with arbitrary scale matrices.
#' 
#' @param X n x p matrix drawn from a n x p matrix normal distribution MN(M,U,Sigma). X must have n rows and p columns.
#' @param U A n x n positive-definite matrix describing the dependence structure between the rows in X. If NULL, observations are considered independent and U is set to the n x x identity matrix.
#' @param Sigma A p x p positive-definite matrix describing the dependence structure between the columns in X. If NULL, Sigma is over-estimated (in the sens of the Loewner partial order).
#' @param Y If Sigma is NULL, an i.i.d. copy of X allowing its estimation. Y must have the same number of columns as X.
#' @param UY If Sigma is NULL, a positive-definite matrix describing the dependence structure between the rows in Y. If NULL and its inverse is not provided, set to the identity matrix by default.
#' @param precUY The inverse matrix of UY, that can be provided to incredevlase computational efficiency. If UY is not NULL and precUY is NULL, precUY is obtained by inversing UY.
#' @param NC The number of clusters to choose.
#' @param clusters A vector of two integers from 1 to NC indicating the pair of clusters whose means have to be compared.
#' @param itermax The iter.max parameter of the k-means algorithm in kmeans_estimation function of KmeansInference package.
#' @param tol The tol_eps parameter of the k-means algorithm in kmeans_estimation function of KmeansInference package.
#' 
#' @return 
#' \itemize{
#'   \item pvalue - The p-value for the difference between cluster means.
#'   \item stat - The test statistic.
#'   \item km - The partition of the n observations retrieved by the clustering algorithm.
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
#' # k-means under the global null hypothesis
#' test.clusters.km(X, U, Sigma, NC = 3, clusters = sample(1:3, 2))
#' # k-means under the global null hypothesis and over-estimation of Sigma
#' test.clusters.km(X, U, Sigma = NULL, Y = Y, NC = 3, clusters = sample(1:3, 2))
#'
#' @references
#'  [1] L. L. Gao, J. Bien, and D. Witten. Selective inference for hierarchical clustering. Journal of the American Statistical Association, 0(0):1–11, 2022.
#'  [2] Y. T. Chen and D. M. Witten. Selective inference for k-means clustering, 2022. arXiv:2203.15267.
#'  
#' @export
#' 

test.clusters.km <- function(X, U = NULL, Sigma = NULL, Y = NULL, UY = NULL, precUY = NULL, NC, clusters, itermax = 10, tol = 1e-6){
  
  #### Initial checks and pre-processing #######################################

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
  
  # K-means clustering from KmeansInference package 
  test_kmeans <- KmeansInference::kmeans_inference(as.matrix(X), k = NC, cluster_1 = clusters[1], cluster_2 = clusters[2], verbose = FALSE, seed = 1234, sig = 1, tol_eps = tol, iter.max = itermax)
  
  #### Test for the difference of cluster means ################################
  
  # Select individuals from each cluster
  km_at_cl <- as.vector(test_kmeans$final_cluster) # Clustering partition
  n1 <- sum(km_at_cl == clusters[1])
  n2 <- sum(km_at_cl == clusters[2])
  nu <- Matrix::Matrix((km_at_cl == clusters[1])/n1 - (km_at_cl == clusters[2])/n2)
  norm2_nu <- 1/n1 + 1/n2
  
  # Difference of cluster means 
  diff_means <- Matrix::Matrix(Matrix::t(nu)%*%X) 
  
  # Computation of the norm \norm{x}_V = \sqrt{x^T V^{-1} x},
  # where V = nu^T U nu Sigma
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
  
  return(list(pvalue = pv, stat = stat_V, km = km_at_cl))}


