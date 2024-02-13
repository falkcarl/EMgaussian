// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// First part of E step for an entire data matrix
// [[Rcpp::export]]
arma::mat imp1matcov(Rcpp::NumericMatrix D, const arma::colvec muest, const arma::mat sigest){  

  arma::mat d(D.begin(), D.rows(), D.cols(), true); // I think false means no copy; here we want one to avoid replacing original data
  
  int N = d.n_rows;

  arma::mat siginv = inv(sigest);
  
  for(int i = 0; i<N; i++){
    arma::uvec io = find_finite(d.row(i));
    arma::uvec im = find_nonfinite(d.row(i));
    arma::uvec idx = {static_cast<arma::uword>(i)};
    
    d.submat(idx,im) = (muest(im) - inv(siginv.submat(im,im)) * siginv.submat(im,io) * (d.submat(idx,io).t() - muest(io))).t();        
  }
  
  return d;

}
      
// Second part of E step for an entire data matrix
// operates directly on T2 to avoid making a copy
// [[Rcpp::export]]
void imp2matcov(Rcpp::NumericMatrix D, const arma::mat sigest, arma::mat& t2){

  arma::mat d(D.begin(), D.rows(), D.cols(), true);
  
  int N = d.n_rows;

  arma::mat siginv = inv(sigest);
  
  for(int i = 0; i<N; i++){
    arma::uvec im = find_nonfinite(d.row(i));
    t2.submat(im,im) += inv(siginv.submat(im,im));
  }

}

// EM cycle all in one shot, covariance matrix parameterization
// [[Rcpp::export]]
Rcpp::List EMcyclecov(const Rcpp::NumericMatrix D, const arma::colvec muest, const arma::mat sigest){
  
  // First part of imputation
  arma::mat dimp = imp1matcov(D, muest, sigest);

  int N = dimp.n_rows;
  
  // Vector of 1's  
  arma::colvec ones = arma::ones<arma::vec>(N);
  
  // Compute T1 matrix
  arma::mat t1 = dimp.t() * ones;

  // Compute T2 matrix, then add stuff to T2
  arma::mat t2 = dimp.t() * dimp;

  // Then, add stuff to T2
  imp2matcov(D, sigest, t2);

  // Now, compute mu
  arma::colvec newmu = t1/N;

  // update S
  arma::mat newS = t2/N - newmu*newmu.t();

  Rcpp::List out;
  out["mu"] = newmu;
  out["S"] = newS;

  return out;
}

// Negative log-likelihood; covariance matrix parameterization
// [[Rcpp::export]]
double nllcov(const arma::mat d, const arma::colvec muest, const arma::mat sigest, int np){
  
  int N = d.n_rows;
  int J = d.n_cols;
  
  double nll = 0;
  
  for(int i = 0; i<N; i++){
    arma::uvec io = find_finite(d.row(i));
    arma::uvec idx = {static_cast<arma::uword>(i)};       
    nll += (0.5*(log(det(sigest.submat(io,io))) + (d.submat(idx,io).t() - muest(io)).t() * inv(sigest.submat(io,io)) * (d.submat(idx,io).t() - muest(io)) + J*log(2*arma::datum::pi))).eval()(0,0);
  }
  
  return nll;
}