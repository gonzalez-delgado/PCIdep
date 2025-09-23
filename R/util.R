
#' Check whether a matrix has the Compound Symmetry (CS) structure.
#'
#' @keywords internal
#'
#' @param U A square matrix
#'
#' @return Whether \code{U} is CS

is.CS <- function(U){
  
  CS <- TRUE
  if(!matrixNormal::is.square.matrix(U)){CS <- FALSE}
  if(length(unique(U[col(U) == row(U)])) != 1){CS <- FALSE}
  if(length(unique(U[col(U) != row(U)])) != 1){CS <- FALSE}
  return(CS)  
  
}

#' Assess whether a pair of clusters remains invariant across two partitions.
#'
#' @keywords internal
#'
#' @param cl A vector of size \eqn{n} with integers in \eqn{1,\ldots,n} defining the classes of the first partition.
#' @param cl_phi A vector of size  \eqn{n} with integers in \eqn{1,\ldots,n} defining the classes of the second partition.
#' @param clusters A vector of two integers chosen among the coordinates of \code{cl}, denoting the pair of clusters in \code{cl} to be analyzed.
#'
#' @return Whether the cluster \code{clusters}[1] and the cluster \code{clusters}[2] in \code{cl} remain invariant in \code{cl_phi}, independently of label modifications.


preserve.cl <- function(cl, cl_phi, clusters, uncl_label = 0) {
  
  new_cl1 <- unique(cl_phi[which(cl == clusters[1])])
  new_cl2 <- unique(cl_phi[which(cl == clusters[2])])
  
  # Remove unclassified points
  new_cl1 <- setdiff(new_cl1, uncl_label)
  new_cl2 <- setdiff(new_cl2, uncl_label)
   
  in_k1 <- length(new_cl1) == 1 # Individuals in the first clusters stay all in the same cluster
  in_k2 <- length(new_cl2) == 1 # Individuals in the second clusters stay all in the same cluster

  out_k1_k2 <- all(! cl_phi[-which(cl %in% c(clusters, uncl_label))] %in% c(new_cl1, new_cl2)) # New individuals are not assigned to the new clusters after perturbation
  
  in_k1 & in_k2 & out_k1_k2
}

