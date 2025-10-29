#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#ifdef lambda
#undef lambda
#endif

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
  double pen  = lambda * arma::norm(beta, 1);
  return loss + pen;
}

// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps){
  // Your function code goes here
  const int n = static_cast<int>(Xtilde.n_rows);
  const int p = static_cast<int>(Xtilde.n_cols);
  const double nd = static_cast<double>(n);
  const int maxit = 10000;
  
  arma::colvec beta = beta_start;
  arma::colvec r = Ytilde - Xtilde * beta;
  double prev_obj = 1e300;
  
  for (int it = 0; it < maxit; ++it) {
    for (int j = 0; j < p; ++j) {
      arma::colvec xj = Xtilde.col(j);
      r += xj * beta[j];
      double rho = arma::dot(xj, r) / nd;
      double bj  = soft_c(rho, lambda);
      r -= xj * bj;
      beta[j] = bj;
    }
    double obj = arma::dot(r, r) / (2.0 * nd) + lambda * arma::norm(beta, 1);
    if (std::abs(prev_obj - obj) < eps) break;
    prev_obj = obj;
  }
  return beta;
}  

// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps){
  // Your function code goes here
  const int p = static_cast<int>(Xtilde.n_cols);
  const int L = static_cast<int>(lambda_seq.n_elem);
  const int n = static_cast<int>(Xtilde.n_rows);
  const double nd = static_cast<double>(n);
  const int maxit = 10000;
  
  arma::mat B(p, L, arma::fill::zeros);
  arma::colvec beta(p, arma::fill::zeros);
  arma::colvec r = Ytilde;
  
  for (int k = 0; k < L; ++k) {
    double lam_here = lambda_seq[k];
    r = Ytilde - Xtilde * beta;
    double prev_obj = 1e300;
    
    for (int it = 0; it < maxit; ++it) {
      for (int j = 0; j < p; ++j) {
        arma::colvec xj = Xtilde.col(j);
        r += xj * beta[j];
        double rho = arma::dot(xj, r) / nd;
        double bj  = soft_c(rho, lam_here);
        r -= xj * bj;
        beta[j] = bj;
      }
      double obj = arma::dot(r, r) / (2.0 * nd) + lam_here * arma::norm(beta, 1);
      if (std::abs(prev_obj - obj) < eps) break;
      prev_obj = obj;
    }
    B.col(k) = beta;
  }
  return B;
}

