#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
//' @useDynLib SA24204138

// Soft-thresholding function
double soft_threshold(double x, double lambda) {
  if (x > lambda) return x - lambda;
  if (x < -lambda) return x + lambda;
  return 0.0;
}

// Optimize principal component vector
arma::vec optimize_beta(const arma::mat& Sigma_hat, arma::vec beta, double lambda, double T, int max_iter = 100, double tol = 1e-6) {
  int p = Sigma_hat.n_cols;
  arma::vec beta_old = beta;
  
  for (int iter = 0; iter < max_iter; ++iter) {
    for (int j = 0; j < p; ++j) {
      double partial_residual = -arma::dot(Sigma_hat.row(j).t(), beta) + Sigma_hat(j, j) * beta(j);
      beta(j) = soft_threshold(partial_residual, lambda) / Sigma_hat(j, j);
    }
    // Enforce L1 norm constraint
    if (arma::norm(beta, 1) > T) {
      beta *= T / arma::norm(beta, 1);
    }
    // Convergence check
    if (arma::norm(beta - beta_old, 2) < tol) break;
    beta_old = beta;
  }
  return beta;
}

// Nodewise Lasso function
arma::mat nodewise_lasso(const arma::mat& A, double lambda, double T) {
  int p = A.n_cols;
  arma::mat Theta = arma::zeros<arma::mat>(p, p);
  
  for (int j = 0; j < p; ++j) {
    arma::vec gamma = arma::zeros<arma::vec>(p - 1);
    arma::vec Aj = A.col(j);
    arma::mat A_minus_j = A;
    A_minus_j.shed_col(j);
    
    double error_tolerance = 1e-6;
    arma::vec gamma_old = gamma;
    double error = 1.0;
    int max_iter = 1000;
    int iter = 0;
    
    while (error > error_tolerance && iter < max_iter) {
      for (int k = 0; k < gamma.n_elem; ++k) {
        double partial_residual = Aj(k) - arma::dot(A_minus_j.row(k).t(), gamma);
        gamma(k) = soft_threshold(partial_residual, lambda) / (A_minus_j(k, k) + 1e-10);
      }
      error = arma::norm(gamma - gamma_old, 2);
      gamma_old = gamma;
      iter++;
    }
    Theta.col(j) = -arma::join_cols(gamma.head(j), arma::vec(1, arma::fill::zeros), gamma.tail(p - j - 1));
  }
  return Theta;
}

// Debiased sparse PCA function
List debiased_sparse_pca(const arma::mat& Sigma_hat, const arma::vec& beta_hat, double lambda, double T, const arma::mat& Theta) {
  double beta_hat_norm_sq = arma::dot(beta_hat, beta_hat);
  arma::vec bias_correction = Theta * (beta_hat_norm_sq * beta_hat - Sigma_hat * beta_hat);
  arma::vec beta_debiased = beta_hat - bias_correction;
  
  double eigenvalue = arma::dot(beta_hat, Sigma_hat * beta_hat);
  
  return List::create(
    Named("beta_debiased") = beta_debiased,
    Named("eigenvalue") = eigenvalue
  );
}

// Alternating optimization main function
// [[Rcpp::export]]
List alternating_debiased_sparse_pca(const arma::mat& X, arma::vec beta_hat, double lambda, double T, int max_iter = 100, double tol = 1e-6) {
  if (X.n_rows < X.n_cols) {
    stop("Input matrix X should have more rows than columns.");
  }
  
  int p = X.n_cols;
  arma::mat Sigma_hat = X.t() * X / X.n_rows;
  beta_hat = optimize_beta(Sigma_hat, beta_hat, lambda, T);
  arma::mat Theta;
  arma::vec beta_debiased;
  double eigenvalue;
  
  for (int iter = 0; iter < max_iter; ++iter) {
    arma::mat Rn_beta = -Sigma_hat + arma::eye<arma::mat>(p, p) * arma::dot(beta_hat, beta_hat) + 2 * beta_hat * beta_hat.t();
    Theta = nodewise_lasso(Rn_beta, lambda, T);
    
    List pca_result = debiased_sparse_pca(Sigma_hat, beta_hat, lambda, T, Theta);
    beta_debiased = as<arma::vec>(pca_result["beta_debiased"]);
    eigenvalue = as<double>(pca_result["eigenvalue"]);
    
    if (arma::norm(beta_debiased - beta_hat, 2) < tol) break;
    beta_hat = beta_debiased;
  }
  
  return List::create(
    Named("beta_debiased") = beta_debiased,
    Named("eigenvalue") = eigenvalue,
    Named("Theta") = Theta
  );
}


