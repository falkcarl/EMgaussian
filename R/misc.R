# Not sure where to put some of these functions. Some may not be used currently,
# but may see use in the near future.

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

# log-likelihood for individual participant
# using precision matrix parameterization
#' @importFrom lavaan lav_matrix_vechr_reverse
lli7 <- function(p,y,ni){
  
  # obtain model-implied stuff
  mu <- p[1:ni]
  K <- lav_matrix_vechr_reverse(p[(ni+1):length(p)])
  Kinv <- solve(K) # need to take inverse first
  
  # based on whatever subject has available
  nnay<-!is.na(y)
  mu<-mu[nnay]
  K<-K[nnay,nnay,drop=FALSE]
  Kinv<-Kinv[nnay,nnay,drop=FALSE]
  y<-y[nnay]
  
  # log likelihood - last part is a constant and therefore omitted from Stadler & Buhlmann
  #out<- - .5*(log(det(sig))+t(y-mu)%*%siginv%*%(y-mu))
  # Stadler & Buhlmann Eq 7
  #out<- - .5*(log(det(Kinv))+t(y-mu)%*%solve(Kinv)%*%(y-mu) + length(p)*log(2*pi))
  out<- - .5*(log(det(Kinv))+t(y-mu)%*%solve(Kinv)%*%(y-mu) + ncol(K)*log(2*pi))  
  return(out)
}

# negative log-likelihood, using precision matrix parameterization
# lambda allows penalization
#' @importFrom lavaan lav_matrix_vechr_reverse
nll7<-function(p, dat, lambda=NULL, pendiag=FALSE){
  ni<-ncol(dat)
  
  out<-sapply(1:nrow(dat),function(i){
    lli7(p,dat[i,,drop=FALSE],ni)
  })
  
  if(!is.null(lambda)){
    K<-p[(ni+1):length(p)]
    if(!pendiag){
      K <- lav_matrix_vechr_reverse(K)
      K <- K[lower.tri(K)]
    }
    pen <- lambda*sum(abs(K))
    out<-out-pen
  }
  
  -sum(out)
  
}

# hack to fix a particular element to a fixed value
# use for profile likelihood/constrained estimation
nll7.fix<-function(p,dat,fix.idx,fix.val,...){
  p[fix.idx] <- fix.val
  return(nll7(p,dat,...))
}