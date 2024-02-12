######################################################################
##
## Copyright 2019-2024 Carl F. Falk
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
## <http://www.gnu.org/licenses/>

#' EM algorithm originally in Stadler & Buhlmann (2012)
#' 
#' @param dat blah
#' @param max.iter blah
#' @param tol blah
#' @param start blah
#' @param debug blah
#' @param ... blah
#' @details
#' Additional details...
#' @importFrom stats cov
#' @importFrom lavaan lavCor lav_matrix_vechr lav_matrix_vechr_reverse
#' @importFrom matrixcalc is.positive.definite is.symmetric.matrix
#' @importFrom Matrix nearPD forceSymmetric
#' @importFrom Rcpp evalCpp 
#' @useDynLib EMgaussian 
#' @export
#' @examples
#' \dontrun{
#'   library(psych)
#'   data(bfi)
#'   test <- em.prec(bfi[,1:25])
#' }
##########################################################################################
# EM algorithm originally in Stadler & Buhlmann (2012), which is on missing data w/ Gaussian graphical model
# Try implementing their algorithm. Section 2.3.2
# Just remove the glasso part and we have a saturated mean and covariance matrix
em.prec <- function(dat, max.iter = 500, tol=1e-5, start=c("diag","pairwise","listwise","full"), debug=0,
                     ...){
  
  start <- match.arg(start)

  dat <- as.matrix(dat)
  
  # obtain starting values: note that lavaan is kind of cheating
  mustart.sb <- colMeans(dat,na.rm=T)
  if(start=="listwise"){
    covstart.sb <- cov(dat, use="complete.obs")
  } else if (start=="pairwise"){
    covstart.sb <- cov(dat,use="pairwise.complete.obs")
  } else if (start=="full"){
    covstart.sb <- lavCor(as.data.frame(dat), missing="ml", output="cov")
  } else if (start=="diag"){
    covstart.sb <- diag(diag(cov(dat,use="pairwise.complete.obs")))
  }
  
  if(any(is.na(covstart.sb))){
    covstart.sb[is.na(covstart.sb)]<-0
    warning("Some NAs in starting vals for covariance matrix. Replacing with zeros.")
  }
  
  # force positive-definite covstart.sb
  if(!is.positive.definite(covstart.sb)){
    covstart.sb<-as.matrix(nearPD(covstart.sb)$mat)
  }

  Kstart.sb.mat<-solve(covstart.sb)
  Kstart.sb <- lav_matrix_vechr(Kstart.sb.mat)
  
  mu.est <- as.matrix(mustart.sb)
  
  K.est <- Kstart.sb
  K.est.mat <- lav_matrix_vechr_reverse(K.est)

  pkstart.sb<-c(mustart.sb,Kstart.sb)
  p.est<-pkstart.sb
  
  ones <- as.matrix(rep(1,nrow(dat)))
  
  # Estep
  # Is this just the same as described in section 2.3
  # Conditional expectation for just xji is the same!
  N<-nrow(dat)
  
  # loop until convergence
  conv<-FALSE
  for(it in 1:max.iter){
    if(debug>1){
      print(paste0('iter ',it))
    }

    d.imp<-imp1mat(dat,mu.est,K.est.mat) # here, K.est.mat is used instead of S.est.mat

    T1 <- t(d.imp)%*%ones
    
    # What about for E[xijxij']?
    # It looks like the imputations can be used directly, except for cases where both values are missing.
    # We need to add K^m whatever.
    # 1. Just do cross product with imputed dataset
    T2 <- t(d.imp)%*%d.imp

    # 2. Then, add stuff to T2. Can loop over individuals
    imp2mat(dat,K.est.mat,T2) # again here

    # Now, compute mu
    mu.est <- T1/N
    
    # update S, then feed to glasso
    S<- (1/N)*T2 - mu.est%*%t(mu.est)
    
    if(!is.symmetric.matrix(S)){
      S<-as.matrix(forceSymmetric(S))
      warnings("Non-symmetric S found after E step. Could indicate identification or estimation problem.")
    }
    
    # force positive-definite S
    if(!is.positive.definite(S)){
      S<-as.matrix(nearPD(S)$mat)
      warnings("Non-positive definite S found after E step. Could indicate identification or estimation problem.")
      
    }

    K.est.mat <- solve(S)
    K.est <- lav_matrix_vechr(K.est.mat)
    
    S.est.mat <- S

    p.old <- p.est
    p.est <- c(mu.est,K.est)
    
    # check for convergence
    if(debug>2){print(max(abs(p.old-p.est)))}
    if(max(abs(p.old-p.est)) < tol){
      conv<-TRUE
      break
    }
  }
  
  out<-list(p.est=p.est, mu=mu.est, S=S.est.mat, K=K.est.mat, it=it, conv=conv)
  return(out)
  
}