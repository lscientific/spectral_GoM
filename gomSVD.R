#' @title gomSVD
#'
#' @description A singular value decomposition algorithm for identifiable GoM model with binary responses
#'
#' @param U the left singular matrix. 
#' @param V the right singular matrix. 
#' @param d the vector containing the singular values. 
#' @param r the number of neighbors. Default value is 10.
#' @param q the probability for the upper quantile for row norms. Default value is 0.4.
#' @param e the probability for the upper quantile for average distance. Default value is 0.2.
#' @param eps the minimum value for item response probabilities. Default value is 0
#'
#' @return The function returns a list with the following components:
#' \describe{
#'   \item{P_hat}{The estimated membership scores.}
#'   \item{T_hat}{The estimated item response parameters.}
#'   \item{R_hat}{The estimated response expectation.}
#'   \item{S}{The estimated indices of pure subjects.}
#'   \item{t}{The computation time.}
#' }

gomSVD <- function(U, V, d, r=10, q=0.4, e=0.2, eps=0) {
  t1 <- Sys.time()
  
  N <- dim(U)[1]
  K <- dim(U)[2]
  S <- pruning(U, r, q, e)  # obtain the subjects to be pruned
  X = U[-S,]
  indices_X <- spa(X)$indices  # SPA to find the pure subjects for the pruned matrix
  indices_U <- ((1:N)[-S])[indices_X]  # the pure subjects for the original matrix
  vertices <- X[indices_X, ]  # the simplex vertices
  
  # estimation for Pi
  P1 <- U %*% solve(vertices)  
  P2 <- t(apply(P1, 1, function(x) ifelse(x<0, 0, x)))  # make every element non-negative
  P_hat <- t(apply(P2, 1, function(x) x / sum(x)))  # re-scale
  
  # estimation for R0
  R_hat = U %*% diag(d) %*% t(V)
  R_hat[R_hat>1-eps] <- 1 - eps
  R_hat[R_hat<eps] <- eps
  
  # estimation for Theta
  T_hat <- t(solve(t(P_hat) %*% P_hat) %*% t(P_hat) %*% U %*% diag(d) %*% t(V))
  T_hat[T_hat > 1-eps] <- 1 - eps
  T_hat[T_hat < eps] <- eps
  
  t2 <- Sys.time()
  
  return(list(P_hat=P_hat, T_hat=T_hat, R_hat=R_hat, 
              S=S, indices_U=indices_U, t=t2-t1))
}

