# Functions for computing (negative) log-likelihood and related stuff

#' @importFrom lavaan lav_matrix_vechr_reverse
#' @useDynLib EMgaussian
nllggm.wrap<-function(p, dat, param = c("prec","cov")){
  
  param <- match.arg(param)
  
  ni <- ncol(dat)
  np <- length(p)
  mu <- p[1:ni]

  mat <- lav_matrix_vechr_reverse(p[(ni+1):np])
  
  if(param == "prec"){
    nll <- nllprec(as.matrix(dat), mu, mat)
  } else if (param == "cov"){
    nll <- nllcov(as.matrix(dat), mu, mat)    
  }

  return(nll)
  
}

#' Information criteria for GGM under missing data
#' 
#' @param dat Data frame or matrix that contains the raw data.
#' @param index Which information criteria to compute? See details.
#' @param mu Numeric vector of means.
#' @param K Precision matrix. If not provided, \code{sig} is required.
#' @param sig Covariance matrix. If not provided, \code{K} is required.
#' @param N Optional sample size for use in information criteria computations.
#' @param gam Single numeric value for \eqn{\gamma} in EBIC computations.
#' @param tol Tolerance for declaring edges to be zero. Anything in absolute value
#'   below \code{tol} will be considered zero and will not count as a parameter.
#' @details Computes information criterion (EBIC, BIC, AIC, or negative log-likelihood)
#'   for the given dataset, which may contain missing values, under the provided
#'   means (\code{mu}) and either precision matrix (\code{K}) or covariance matrix
#'   (\code{sig}) and assuming multivariate normality. Provide only \code{K} or
#'   \code{sig}, not both. Sample size computations are based on average pair-wise
#'   observations, unless a value for \code{N} is provided. Note that \code{gam}
#'   is the value for \eqn{\gamma} in EBIC computations. The typical choice is .5,
#'   but will be set to 0 if BIC is desired. Equations used are as follows:
#'   
#'   \deqn{EBIC = -2ll + E\log(N) + 4\gamma E \log(P)}
#'   \deqn{AIC = -2ll + 2E}
#'   
#'   Where \eqn{ll} is the log-likelihood, \eqn{E} is the number of edges (i.e., number
#'   of estimated parameters, though mean structure is ignored), \eqn{N} is the sample
#'   size, and \eqn{P} is the number of nodes (i.e., columns in the dataset).
#'   
#' @return Single numeric value corresponding to the selected information criterion.
#' @export
ICggm <- function(dat, index = c("EBIC","BIC","AIC","nll"), mu, K=NULL, sig=NULL, N=NULL, gam=.5, tol=1e-32){
  
  # set up index
  index <- match.arg(index)
  
  if(index=="BIC"){
    gam <- 0
  }
  
  # which matrices we have
  if(is.null(K) & is.null(sig)){
    stop("K or sig must be specified")
  }
  if(!is.null(K) & !is.null(sig)){
    stop("Only specify K or sig, not both")
  }
  
  # Compute -2ll and number of off-diagonal elements that are non-zero for either matrix
  if(!is.null(sig)){
    neg2l <- 2*nllcov(as.matrix(dat), mu, sig)
    E<-sum(abs(sig[lower.tri(sig)])>tol)    
  } else {
    neg2l <- 2*nllprec(as.matrix(dat), mu, K)
    E<-sum(abs(K[lower.tri(K)])>tol)
  }
  
  # set up N
  if(is.null(N)){
    # next 3 lines copied and modified from bootnet, then modified
    xmat <- as.matrix(!is.na(dat))
    misMatrix <- t(xmat) %*% xmat
    N<-mean(misMatrix[lower.tri(misMatrix)])
  }

  
  # number of "nodes"
  P <- ncol(dat)
  
  # how do we know how many non-zero edges there are?
  # actually look at off-diagonal elements of whatever matrix
  #K <- lav_matrix_vechr_reverse(p[(P+1):length(p)])
  #E<-sum(abs(K[lower.tri(K)])>tol)
  
  if(index %in% c("EBIC","BIC")){
    out <- neg2l + E*log(N) + 4*gam*E*log(P)
  } else if (index == "AIC"){
    out <- neg2l + 2*E
  } else if (index == "nll"){
    out <- neg2l
  }
  
  return(out)
}

#' Penalized negative log-likelihood, parameter vector as input
#' 
#' @param p Parameter vector for means, then precision or covariance matrix
#'   (rows of lower triangular portion).
#' @param dat Data frame or matrix that contains the raw data.
#' @param param Parameterization to use: Precision matrix or covariance matrix.
#' @param fix.idx Optional vector of parameter index numbers that will be fixed
#'   to particular values.
#' @param fix.val Optional vector of parameter values; should match length of
#'   \code{fix.idx}.
#' @param lambda Optional numeric value for tuning parameter for L1 penalty
#'   term.
#' @param pendiag Boolean value indicating whether to penalize the diagonal of
#'   the precision matrix.
#' @details Assuming multivariate normality and accommodating missing data by
#'   using all available information (e.g., as would be used in direct maximum
#'   likelihood), computes value of L1 penalized log-likelihood.
#'   Off-diagonal elements of the matrix (precision or covariance)
#'   are penalized, with the option to also penalize the diagonal elements.
#'   
#'   Let \eqn{\theta} contain all model parameters, and \eqn{l()} be the
#'   log-likelihood, then the negative penalized log-likelihood is:
#'   \deqn{-l(\theta) + \lambda \sum |k_{jj'} |}
#'   
#'   where \eqn{\lambda} is the tuning parameter and the sum is over all unique
#'   elements in the matrix (precision or covariance), which, depending on the
#'   options chosen may or may not include the diagonal.
#'   
#'   This function uses a parameter vector as input, thereby possibly
#'   facilitating optimization. It also supports fixing particular parameters
#'   indexed by \code{fix.idx} to particular values in \code{fix.val}.
#'   
#'   
#' @return Single numeric value corresponding to the negative penalized
#'   log-likelihood.
#' @importFrom lavaan lav_matrix_vechr_reverse
#' @export
nll.parpen <- function(p, dat, param = c("prec","cov"), fix.idx=NULL, fix.val=NULL,
                    lambda = NULL, pendiag = FALSE){
  param <- match.arg(param)
  
  p[fix.idx] <- fix.val
  
  out <- nllggm.wrap(p=p, dat=dat, param=param)
  
  if(!is.null(lambda)){
    ni <- ncol(dat)
    K<-p[(ni+1):length(p)]
    if(!pendiag){
      K <- lav_matrix_vechr_reverse(K)
      K <- K[lower.tri(K)]
    }
    pen <- lambda*sum(abs(K))
    out<-out-pen
  }
  return(out)
}