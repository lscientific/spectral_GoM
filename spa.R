#' @title SPA
#'
#' @description A sequential projection algorithm to find the pure subjects
#'
#' @param mat the left singular matrix. 

#' @return The function returns a vector of the pure subject indices.


spa <- function(mat) {
  N <- dim(mat)[1]
  K <- dim(mat)[2]
  indices = rep(0, K)
  
  Y <- mat  # make a copy of the left singular matrix
  for (k in 1:K) {
    row_norms <- apply(Y, 1, function(x) sqrt(sum(x^2)))  # calculate the row norms
    idx <- which.max(row_norms)  # identify the largest norm
    indices[k] = idx
    u = Y[idx, ] / row_norms[idx]
    Y = Y %*% (diag(K) - matrix(u, nc=1) %*% matrix(u, nr=1))  # projection into the orthogonal subspace
  }
  
  return(indices)
}
