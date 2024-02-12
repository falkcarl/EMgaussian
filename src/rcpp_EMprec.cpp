// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat imp1mat(Rcpp::NumericMatrix D, const arma::colvec muest, const arma::mat kest){
    
  arma::mat d(D.begin(), D.rows(), D.cols(), true); // I think false means no copy; here we want one to avoid replacing original data
  
  int N = d.n_rows;

  for(int i = 0; i<N; i++){
    arma::uvec io = find_finite(d.row(i));
    arma::uvec im = find_nonfinite(d.row(i));
    arma::uvec idx = {static_cast<arma::uword>(i)};

    d.submat(idx,im) = (muest(im) - inv(kest.submat(im,im)) * kest.submat(im,io) * (d.submat(idx,io).t() - muest(io))).t();        
  }
  return(d);
}
      
// Second part of E step for an entire data matrix
// operates directly on T2 to avoid making a copy
// [[Rcpp::export]]
void imp2mat(const arma::mat d, const arma::mat kest, arma::mat t2){
  
  int N = d.n_rows;
  //arma::mat t2(T2.begin(), T2.rows(), T2.cols(), false);
  
  for(int i = 0; i<N; i++){
    arma::uvec im = find_nonfinite(d.row(i));
    t2.submat(im,im) += inv(kest.submat(im,im));
  }
  
}
        
// Negative log-likelihood; precision matrix parameterization
// [[Rcpp::export]]
double nllprec(const arma::mat d, const arma::colvec muest, const arma::mat kest, int np){
    
    int N = d.n_rows;
    int J = d.n_cols;
    
    arma::mat kinv = inv(kest);
    
    double nll = 0;
    
    for(int i = 0; i<N; i++){
      arma::uvec io = find_finite(d.row(i));
      arma::uvec idx = {static_cast<arma::uword>(i)};
      nll += (0.5*(log(det(kinv.submat(io,io))) + (d.submat(idx,io).t() - muest(io)).t() * inv(kinv.submat(io,io)) * (d.submat(idx,io).t() - muest(io)) + J*log(2*arma::datum::pi))).eval()(0,0);
    }
    
    return(nll);
}