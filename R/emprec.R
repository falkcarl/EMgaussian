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

  # loop until convergence
  conv<-FALSE
  for(it in 1:max.iter){
    if(debug>1){
      print(paste0('iter ',it))
    }

    # Do one EM cycle
    EMres <- EMcycleprec(dat, mu.est, K.est.mat)
    mu.est <- EMres$mu
    S.est.mat <- EMres$S
    
    # Check S
    if(!is.symmetric.matrix(S.est.mat)){
      S.est.mat<-as.matrix(forceSymmetric(S.est.mat))
      warnings("Non-symmetric S found after E step. Could indicate identification or estimation problem.")
    }
    
    # force positive-definite S
    if(!is.positive.definite(S.est.mat)){
      S.est.mat<-as.matrix(nearPD(S.est.mat)$mat)
      warnings("Non-positive definite S found after E step. Could indicate identification or estimation problem.")
    }

    # Save estimates
    K.est.mat <- solve(S.est.mat)
    K.est <- lav_matrix_vechr(K.est.mat)
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