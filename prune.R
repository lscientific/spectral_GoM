#' @title Pruning 
#'
#' @description Pruning step to remove noisy points from the empirical simplex
#'
#' @param mat the left singular matrix. 
#' @param r the number of neighbors. Default value is 10.
#' @param q the probability for the upper quantile for row norms. Default value is 0.4.
#' @param e the probability for the upper quantile for average distance. Default value is 0.2.
#'
#' @return The function returns the index vector of the rows to be pruned from the left singular matrix.
#'

pruning <- function(mat, r=10, q=0.4, e=0.2) {
  N <- dim(mat)[1]
  row_norms <- apply(mat, 1, function(x) sqrt(sum(x^2)))
  S0 <- which(row_norms > quantile(row_norms, 1-q))

  x = c()  # x records the average distance to the r neighbors
  for (s in S0){
    mat_s <- t(apply(mat, 1, function(x) x - mat[s,]))
    norm_s <- apply(mat_s, 1, function(x) sqrt(sum(x^2)))
    d = norm_s[sort(norm_s, index.return=T)$ix[2:(r+1)]]
    x <- c(x, mean(d))
  }
  S <- S0[which(x > quantile(x, 1-e))]

  return(S)
}

