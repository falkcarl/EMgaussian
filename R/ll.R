# Functions for computing (negative) log-likelihood and related stuff

#' @importFrom lavaan lav_matrix_vechr_reverse
#' @useDynLib EMgaussian
nllggm.wrap<-function(p, dat){
  
  ni <- ncol(dat)
  np <- length(p)
  mu <- p[1:ni]
  K <- lav_matrix_vechr_reverse(p[(ni+1):np])
  
  nll <- nllprec(as.matrix(dat), mu, K)
  return(nll)
  
}

#TODO: functions that fix particular elements of the precision matrix to a value,
# with optional tuning parameter, and compute log-likelihood. Use same nomenclature as
# zero argument in glasso package

#TODO: change below so that it calls nllggm.wrap or it just calls nll prec directly
# hack to fix a particular element to a fixed value
# use for profile likelihood/constrained estimation
nll7.fix<-function(p,dat,fix.idx,fix.val,...){
  p[fix.idx] <- fix.val
  return(nll7(p,dat,...))
}