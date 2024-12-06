#' Alternating De-biased Sparse PCA
#'
#' This function implements the de-biased sparse principal component analysis (PCA) algorithm, 
#' which alternates between optimizing the principal component vector and estimating 
#' the inverse covariance matrix using nodewise Lasso. It produces a sparse, de-biased 
#' estimate of the principal component and the associated eigenvalue.
#'
#' @param X A numeric matrix of size \(n \\times p\), where \(n\) is the number of samples, 
#'   and \(p\) is the number of features.
#' @param beta_hat A numeric vector of length \(p\), representing the initial guess 
#'   for the principal component vector. The user must provide this value.
#' @param lambda A positive numeric value representing the regularization parameter for 
#'   the sparse optimization.
#' @param T A positive numeric value representing the constraint on the \(L_1\)-norm of 
#'   the principal component vector.
#' @param max_iter An integer specifying the maximum number of iterations for the 
#'   optimization process. Default is 100.
#' @param tol A positive numeric value representing the convergence threshold for 
#'   stopping the optimization. Default is \eqn{10^{-6}}.
#'
#' @return A list containing the following components:
#'   \item{beta_debiased}{A numeric vector representing the de-biased sparse estimate 
#'     of the principal component vector.}
#'   \item{eigenvalue}{A numeric value representing the estimated eigenvalue 
#'     associated with the de-biased principal component.}
#'   \item{Theta}{A numeric matrix representing the estimated inverse covariance matrix.}
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib SA24204138, .registration = TRUE
#'
#' @author SA24204138
#'
#' @export
alternating_debiased_sparse_pca <- function(X, beta_hat, lambda, T, max_iter = 100, tol = 1e-6) {
  # Call the C++ implementation
  .Call("_SA24204138_alternating_debiased_sparse_pca", PACKAGE = "SA24204138", X, beta_hat, lambda, T, max_iter, tol)
}