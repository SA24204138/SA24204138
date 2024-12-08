#' @title De-biased Sparse PCA
#' @description Implements the de-biased sparse PCA algorithm.
#' @param X Numeric matrix of size \(n \times p\).
#' @param beta_hat Initial guess for the principal component vector.
#' @param lambda Regularization parameter.
#' @param T Constraint on \(L_1\)-norm.
#' @param max_iter Maximum number of iterations (default = 100).
#' @param tol Convergence threshold (default = 1e-06).
#' @return A list with `beta_debiased`, `eigenvalue`, and `Theta`.
#' @export
alternating_debiased_sparse_pca <- function(X, beta_hat, lambda, T, max_iter = 100, tol = 1e-6) {
  # Call the C++ implementation
  .Call("_SA24204138_debiased_sparse_pca", PACKAGE = "SA24204138", X, beta_hat, lambda, T, max_iter, tol)
}