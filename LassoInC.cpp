library(Rcpp); library(RcppArmadillo)
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Soft-thresholding function, returns scalar
// [[Rcpp::export]]
double soft_c(double a, double lambda){
  // Your function code goes here
  if (a >  lambda) return a - lambda;
  if (a < -lambda) return a + lambda;
  return 0.0;
}

// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Your function code goes here
  const double n = static_cast<double>(Xtilde.n_rows);
  arma::colvec r = Ytilde - Xtilde * beta;
  double loss = arma::dot(r, r) / (2.0 * n);
  double pen  = lambda * arma::accu(arma::abs(beta));
  return loss + pen;
}

// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){
  // Your function code goes here
  const double n = static_cast<double>(Xtilde.n_rows);
  arma::colvec r = Ytilde - Xtilde * beta;
  double loss = arma::dot(r, r) / (2.0 * n);
  double pen  = lambda * arma::accu(arma::abs(beta));
  return loss + pen;
}

// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){
  // Your function code goes here
  const int n = static_cast<int>(Xtilde.n_rows);
  const int p = static_cast<int>(Xtilde.n_cols);
  const double nd = static_cast<double>(n);
  const int maxit = 10000;
  
  arma::colvec beta = beta_start;               // no input checks per assignment
  arma::colvec r = Ytilde - Xtilde * beta;      // residuals
  
  double prev_obj = std::numeric_limits<double>::infinity();
  
  for (int it = 0; it < maxit; ++it) {
    for (int j = 0; j < p; ++j) {
      const arma::colvec xj = Xtilde.col(j);
      r += xj * beta[j];                         // add back old contribution
      const double rho = arma::dot(xj, r) / nd;  // (1/n) x_j^T r  (standardized)
      const double bj = soft_c(rho, lambda);     // unit denom
      r -= xj * bj;                              // subtract new contribution
      beta[j] = bj;
    }
    const double obj = arma::dot(r, r) / (2.0 * nd) + lambda * arma::accu(arma::abs(beta));
    if (std::abs(prev_obj - obj) < eps) break;  // stop at first time < eps
    prev_obj = obj;
  }
  return beta;
}  

// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  // Your function code goes here
  const int p = static_cast<int>(Xtilde.n_cols);
  const int L = static_cast<int>(lambda_seq.n_elem);
  const int n = static_cast<int>(Xtilde.n_rows);
  const double nd = static_cast<double>(n);
  const int maxit = 10000;
  
  arma::mat B(p, L, arma::fill::zeros);
  arma::colvec beta(p, arma::fill::zeros);   // warm start across lambdas
  arma::colvec r = Ytilde;                   // residual for beta=0
  
  for (int k = 0; k < L; ++k) {
    const double lambda = lambda_seq[k];
    r = Ytilde - Xtilde * beta;              // ensure consistency
    double prev_obj = std::numeric_limits<double>::infinity();
    
    for (int it = 0; it < maxit; ++it) {
      for (int j = 0; j < p; ++j) {
        const arma::colvec xj = Xtilde.col(j);
        r += xj * beta[j];
        const double rho = arma::dot(xj, r) / nd;
        const double bj = soft_c(rho, lambda);
        r -= xj * bj;
        beta[j] = bj;
      }
      const double obj = arma::dot(r, r) / (2.0 * nd) + lambda * arma::accu(arma::abs(beta));
      if (std::abs(prev_obj - obj) < eps) break;
      prev_obj = obj;
    }
    B.col(k) = beta;
  }
  return B;
}