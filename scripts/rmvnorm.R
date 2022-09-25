#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 3: An Overview of Linear Regression 
#        and Data Generation
#
#' Generate data for multiple regression
#'
#' @param n sample size
#' @param mu mean vector, p x 1 
#' @param Sigma variance-covariance matrix, p x p 
#' @return Data for multiple regression model, n x p 
rmvnorm <- function(n, mu, Sigma) {
  
  p <- length(mu)
  
  if (dim(mu) == c(p, 1) && dim(Sigma) == c(p, p)) {
    
    Z <- matrix(rnorm(p * n, 0, 1), ncol = p)
    X <-  matrix(1, n) %*% t(mu) + Z %*% chol(Sigma)
    return(X)
    
  }else{
    stop("Input dimensions should match!")
  }
}

