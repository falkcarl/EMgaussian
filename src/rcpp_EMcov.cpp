// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// First part of E step for an entire data matrix
// [[Rcpp::export]]
arma::mat imp1matsigma(Rcpp::NumericMatrix D, Rcpp::NumericVector Muest, Rcpp::NumericMatrix Sigest){

  arma::mat d(D.begin(), D.rows(), D.cols(), true); // I think false means no copy
  
  int N = d.n_rows;
  //int J = d.n_cols;
  
  arma::colvec muest(Muest.begin(),Muest.length(),false);
  arma::mat sigest(Sigest.begin(), Sigest.rows(), Sigest.cols(), false);
  
  arma::mat siginv = inv(sigest);
  
  for(int i = 0; i<N; i++){
    arma::uvec io = find_finite(d.row(i));
    arma::uvec im = find_nonfinite(d.row(i));
    //arma::uvec idx = {i};
    arma::uvec idx = {static_cast<arma::uword>(i)};    
    
    d.submat(idx,im) = (muest(im) - inv(siginv.submat(im,im)) * siginv.submat(im,io) * (d.submat(idx,io).t() - muest(io))).t();        
  }
  
  return(d);

}
      
// Second part of E step for an entire data matrix
// operates directly on T2 to avoid making a copy
// [[Rcpp::export]]
void imp2matsigma(Rcpp::NumericMatrix D, Rcpp::NumericMatrix Sigest, Rcpp::NumericMatrix T2){

  arma::mat d(D.begin(), D.rows(), D.cols(), false);

  int N = d.n_rows;
  //int J = d.n_cols;
  
  arma::mat sigest(Sigest.begin(), Sigest.rows(), Sigest.cols(), false);
  arma::mat t2(T2.begin(), T2.rows(), T2.cols(), false);
  
  arma::mat siginv = inv(sigest);
  
  for(int i = 0; i<N; i++){
    //arma::uvec io = find_finite(d.row(i));
    arma::uvec im = find_nonfinite(d.row(i));
    //arma::uvec idx = {i};
    
    t2.submat(im,im) += inv(siginv.submat(im,im));
  }

}
        
// Negative log-likelihood; covariance matrix parameterization
// [[Rcpp::export]]
double nllmvn(Rcpp::NumericMatrix D, Rcpp::NumericVector Muest, Rcpp::NumericMatrix Sigest, int np){
  
  arma::mat d(D.begin(), D.rows(), D.cols(), false);
  
  int N = d.n_rows;
  int J = d.n_cols;
  
  arma::colvec muest(Muest.begin(),Muest.length(),false);
  arma::mat sigest(Sigest.begin(), Sigest.rows(), Sigest.cols(), false);
  
  double nll = 0;
  
  for(int i = 0; i<N; i++){
    arma::uvec io = find_finite(d.row(i));
    //arma::uvec idx = {i};
    arma::uvec idx = {static_cast<arma::uword>(i)};       
    
    nll += (0.5*(log(det(sigest.submat(io,io))) + (d.submat(idx,io).t() - muest(io)).t() * inv(sigest.submat(io,io)) * (d.submat(idx,io).t() - muest(io)) + J*log(2*arma::datum::pi))).eval()(0,0);
  }
  
  return(nll);
}