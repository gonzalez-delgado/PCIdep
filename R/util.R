
#' Check whether a matrix has the Compound Symmetry (CS) structure
#'
#' @keywords internal
#'
#' @param U a square matrix
#'
#' @return Whether U is CS

is.CS <- function(U){
  
  CS <- TRUE
  if(!matrixNormal::is.square.matrix(U)){CS <- FALSE}
  if(length(unique(U[col(U) == row(U)])) != 1){CS <- FALSE}
  if(length(unique(U[col(U) != row(U)])) != 1){CS <- FALSE}
  return(CS)  
  
}