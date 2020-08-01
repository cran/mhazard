#' Estimates the survival function for a bivariate outcome
#'
#' Computes the survival function for a bivariate outcome using one of
#' three possible estimators. The survival function for a bivariate outcome
#' is analogous to the Kaplan-Meier estimator for a univariate outcome.
#' Optionally (bootstrap) confidence intervals for the survival function
#' may also be computed.
#'
#' @param Y1,Y2 Vectors of event times (continuous).
#' @param Delta1,Delta2 Vectors of censoring indicators (1=event,
#' 0=censored).
#' @param newT1,newT2 Optional vectors of new times at which to estimate
#' the survival function. Defaults to the unique values in Y1/Y2 if not
#' specified.
#' @param estimator Which estimator of the survival function should be used.
#' Possible values include "dabrowska", "volterra", and "prentice-cai".
#' Defaults to "dabrowska".
#' @param conf.int Should bootstrap confidence intervals be computed?
#' @param R Number of bootstrap replicates. This argument is passed to the
#' boot function. Defaults to 1000. Ignored if conf.int is FALSE.
#' @param ... Additional arguments to the boot function.
#' @section Details:
#' If conf.int is TRUE, confidence intervals will be computed using the
#' boot function in the boot package. Currently only 95% confidence
#' intervals computed using the percentile method are implemented. If
#' conf.int is FALSE, confidence intervals will not be computed, and
#' confidence bounds will not be returned in the output.
#' @return A list containing the following elements:
#' \describe{
#' \item{T1:}{Unique values of Y1 at which Fhat was computed}
#' \item{T2:}{Unique values of Y2 at which Fhat was computed}
#' \item{Fhat:}{Estimated survival function (computed at T1, T2)}
#' \item{Fhat.lci:}{Lower 95% confidence bounds for Fhat}
#' \item{Fhat.uci:}{Upper 95% confidence bounds for Fhat}
#' \item{Fmarg1:}{Estimated marginal survival function for variable 1
#' (computed at newT1)}
#' \item{Fmarg1.lci:}{Lower 95% confidence bounds for Fmarg1}
#' \item{Fmarg1.uci:}{Upper 95% confidence bounds for Fmarg1}
#' \item{Fmarg2:}{Estimated marginal survival function for variable 2
#' (computed at newT2)}
#' \item{Fmarg2.lci:}{Lower 95% confidence bounds for Fmarg2}
#' \item{Fmarg2.uci:}{Upper 95% confidence bounds for Fmarg2}
#' \item{Fhat_est:}{Estimated survival function (computed at newT1, newT2)}
#' \item{Fhat_est.lci:}{Lower 95% confidence bounds for Fhat_est}
#' \item{Fhat_est.uci:}{Upper 95% confidence bounds for Fhat_est}
#' \item{CR:}{Estimated cross ratio (computed at T1, T2)}
#' \item{KT:}{Estimated Kendall\'s tau (computed at T1, T2)}
#' \item{CR_est:}{Estimated cross ratio (computed at newT1, newT2)}
#' \item{KT_est:}{Estimated Kendall\'s tau (computed at newT1, newT2)}
#' }
#' @seealso \code{\link[boot]{boot}}
#' @references
#' Prentice, R., Zhao, S. "Nonparametric estimation of the multivariate
#' survivor function: the multivariate Kaplan–Meier estimator", Lifetime
#' Data Analysis (2018) 24:3-27.
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
#' @useDynLib mhazard
#' @importFrom Rcpp sourceCpp
#' @importFrom boot boot
#' @export
#' @examples
#' x <- genClayton2(1000, 0, 1, 1, 2, 2)
#' x.km2 <- KM2(x$Y1, x$Y2, x$Delta1, x$Delta2)
#' \donttest{x.km2.ci <- KM2(x$Y1, x$Y2, x$Delta1, x$Delta2, conf.int=TRUE)}
KM2 <- function(Y1, Y2, Delta1, Delta2, newT1=NULL, newT2=NULL,
	        estimator=c("dabrowska", "volterra", "prentice-cai"),
		conf.int=FALSE, R=1000, ...){
    if (min(c(Y1, Y2))<0) {
        stop("Y1, Y2 must be nonnegative")
    }
    if (sum(!(c(Delta1, Delta2) %in% c(0, 1)))>0) {
        stop("Delta1, Delta2 must be 0/1 vectors")
    }

	  # The estimator seems to have problems if event times are tied.
	  # We can avoid this issue by discarding duplicated event times.
	  # This fixes the issue and will not affect the final estimates.
          T1 <- unique(sort(Y1[Delta1==1]))
          T2 <- unique(sort(Y2[Delta2==1]))

	  if (is.null(newT1)) {
	    newT1 <- T1
	  }
	  if (is.null(newT2)) {
	    newT2 <- T2
	  }

    if (min(c(newT1, newT2))<0) {
        stop("newT1, newT2 must be nonnegative")
    }

	  estimator <- match.arg(estimator)
	  if (estimator=="dabrowska") {
	    method <- 1
	  }
	  else if (estimator=="volterra") {
	    method <- 2
	  }
	  else if (estimator=="prentice-cai") {
	    method <- 3
	  }

	  junk <- jointHazLambda(Y1, Y2, T1, T2, Delta1, Delta2, method)
	  Fhat <- junk[,,1]
	  temp <- calcTemp2(junk[,,1], junk[,,2], junk[,,3], junk[,,4])
	  CR <- temp[,,1]/temp[,,2]
	  KT <- temp[,,3]/temp[,,4]

	  # Calculate confidence intervals. I hard coded percentile
	  # intervals for performance reasons. I examined BCa intervals
	  # (which have better theoretical properties), but it was
	  # extremely slow. Should we consider giving users that as an
	  # option?
	  if (conf.int) {
	    Fhat.boot <- boot::boot(cbind(Y1, Y2, Delta1, Delta2), KM2.boot,
                                    R=R, method=method, ...)
	    Fhat.ci <- getBootCI(Fhat.boot$t)
	    Fhat.lci <- matrix(Fhat.ci[1,], ncol=ncol(Fhat), nrow=nrow(Fhat))
	    Fhat.uci <- matrix(Fhat.ci[2,], ncol=ncol(Fhat), nrow=nrow(Fhat))
	  }

          ##########summary results at specified time points
	  index.i <- sapply(newT1, function(x1, x2) {sum(x2<=x1)},
	                    x2=T1)
	  index.j <- sapply(newT2, function(x1, x2) {sum(x2<=x1)},
	                    x2=T2)

          ###########estimated marginal survival probability#########
          Fmarg1 <- Fhat[index.i+1, 1]
          Fmarg2 <- Fhat[1, index.j+1]

	  if (!conf.int) {
          return(list(Fhat=Fhat, T1=T1, T2=T2, CR=CR, KT=KT, Fmarg1=Fmarg1,
	  	      Fmarg2=Fmarg2, Fhat_est=Fhat[index.i+1, index.j+1],
		      CR_est=CR[index.i+1, index.j+1],
		      KT_est=KT[index.i+1, index.j+1]))
	  }
	  else {
	  Fmarg1.lci <- Fhat.lci[index.i+1, 1]
	  Fmarg1.uci <- Fhat.uci[index.i+1, 1]
	  Fmarg2.lci <- Fhat.lci[1, index.j+1]
	  Fmarg2.uci <- Fhat.uci[1, index.j+1]
          return(list(Fhat=Fhat, Fhat.lci=Fhat.lci, Fhat.uci=Fhat.uci, T1=T1,
	  	      T2=T2, CR=CR, KT=KT, Fmarg1=Fmarg1, Fmarg2=Fmarg2,
		      Fmarg1.lci=Fmarg1.lci, Fmarg1.uci=Fmarg1.uci,
		      Fmarg2.lci=Fmarg2.lci, Fmarg2.uci=Fmarg2.uci,
		      Fhat_est=Fhat[index.i+1, index.j+1],
		      Fhat_est.lci=Fhat.lci[index.i+1, index.j+1],
		      Fhat_est.uci=Fhat.uci[index.i+1, index.j+1],
		      CR_est=CR[index.i+1, index.j+1],
		      KT_est=KT[index.i+1, index.j+1]))
	  }
}

KM2.boot <- function(x, ndx, method) {
  # replace the original data with the bootstrap sample
  Y1 <- x[ndx,1]
  Y2 <- x[ndx,2]
  Delta1 <- x[ndx,3]
  Delta2 <- x[ndx,4]
  # but estimate the survival function at the original T's (rather than
  # at the bootstrapped version of the T's)
  newT1 <-  unique(sort(x[x[,3]==1,1]))
  newT2 <-  unique(sort(x[x[,4]==1,2]))

          T1 <- unique(sort(Y1[Delta1==1]))
          T2 <- unique(sort(Y2[Delta2==1]))

	  Fhat <- jointHazLambda(Y1, Y2, T1, T2, Delta1, Delta2, method)[,,1]

	  # find the indices of the bootstrapped T's that correspond to
	  # the original T's
	  index.i <- sapply(newT1, function(x1, x2) {sum(x2<=x1)},
	                    x2=T1)
	  index.j <- sapply(newT2, function(x1, x2) {sum(x2<=x1)},
	                    x2=T2)
	  return(as.vector(Fhat[c(1,index.i+1),c(1,index.j+1)]))
}
