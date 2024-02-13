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

#' EM algorithm using covariance matrix parameterization
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
#'   test <- em.cov(bfi[,1:25])
#' }
#' 
##########################################################################################
# EM algorithm originally in Stadler & Buhlmann (2012), which is on missing data w/ Gaussian graphical model
# Try implementing their algorithm. Section 2.3.2
# Just remove the glasso part and we have a saturated mean and covariance matrix
em.cov <- function(dat, max.iter = 500, tol=1e-5, start=c("diag","pairwise","listwise","full"), debug=0,
                         ...){
  
  start <- match.arg(start)
  
  dat <- as.matrix(dat)
  
  # obtain starting values
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
  
  # track change in means and covariances
  mu.est <- as.matrix(mustart.sb)
  S.est.mat <- covstart.sb
  S.est <- lav_matrix_vechr(S.est.mat)
  p.est <- c(mustart.sb, S.est)

  # loop until convergence
  conv<-FALSE
  for(it in 1:max.iter){
    if(debug>1){
      print(paste0('iter ',it))
    }
    
    # Do one EM cycle
    EMres <- EMcyclecov(dat, mu.est, S.est.mat)
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
    S.est <- lav_matrix_vechr(S.est.mat)
    p.old <- p.est
    p.est <- c(mu.est,S.est)
    
    # check for convergence
    if(debug>2){print(max(abs(p.old-p.est)))}
    if(max(abs(p.old-p.est)) < tol){
      conv<-TRUE
      break
    }
  }
  
  out<-list(p.est=p.est, mu=mu.est, S=S.est.mat, it=it, conv=conv)
  return(out)
  
}


# likelihood contribution for single respondent
lli <- function(p,y,ni){
  
  # obtain model-implied stuff
  mu <- p[1:ni]
  sig <- lav_matrix_vechr_reverse(p[(ni+1):length(p)])
  
  # based on whatever subject has available
  nnay<-!is.na(y)
  mu<-mu[nnay]
  sig<-sig[nnay,nnay]
  y<-y[nnay]
  
  siginv<-solve(sig)
  
  # log likelihood
  out<-.5*log(det(siginv))-.5*t(y-mu)%*%siginv%*%(y-mu) - .5*length(p)*log(2*pi)
  return(out)
}

# negative log-likelihood
# lambda allows penalization of correlations
nll <- function(p, dat, lambda=NULL, pendiag=FALSE){
  ni<-ncol(dat)
  
  out<-sapply(1:nrow(dat),function(i){
    lli(p,dat[i,,drop=FALSE],ni)
  })
  
  if(!is.null(lambda)){
    sig<-p[(ni+1):length(p)]
    if(!pendiag){
      sig <- lav_matrix_vechr_reverse(sig)
      sig <- sig[lower.tri(sig)]      
    }
    pen <- lambda*sum(abs(sig))
    out<-out-pen
  }
  
  -sum(out)
}

# hack to fix a particular element to a fixed value
# use for profile likelihood/constrained estimation
nll.fix<-function(p,dat,fix.idx,fix.val,...){
  p[fix.idx] <- fix.val
  return(nll(p,dat,...))
}