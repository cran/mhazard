#' Nonparametric estimates of the survival function for bivariate failure
#' time data
#'
#' Computes the survival function for bivariate failure time data using one
#' of three possible estimators, including Dabrowska, Volterra and
#' Prentice-Cai estimators. Optionally (bootstrap) confidence intervals for
#' the survival function may also be computed. 
#'
#' @param Y1,Y2 Vectors of event times (continuous).
#' @param Delta1,Delta2 Vectors of censoring indicators (1=event,
#' 0=censored).
#' @param newT1,newT2 Optional vectors of times at which to estimate
#' the survival function (which do not need to be subsets of Y1/Y2).
#' Defaults to the unique values in Y1/Y2 if not specified.
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
#' \item{T1:}{Unique uncensored Y1 values}
#' \item{T2:}{Unique uncensored Y2 values}
#' \item{Fhat:}{Estimated bivariate survival function (computed at T1, T2)}
#' \item{Fhat.lci:}{Lower 95% confidence bounds for Fhat}
#' \item{Fhat.uci:}{Upper 95% confidence bounds for Fhat}
#' \item{Fmarg1.est:}{Estimated marginal survival function for T1
#' (computed at newT1)}
#' \item{Fmarg1.lci:}{Lower 95% confidence bounds for Fmarg1}
#' \item{Fmarg1.uci:}{Upper 95% confidence bounds for Fmarg1}
#' \item{Fmarg2.est:}{Estimated marginal survival function for T2
#' (computed at newT2)}
#' \item{Fmarg2.lci:}{Lower 95% confidence bounds for Fmarg2}
#' \item{Fmarg2.uci:}{Upper 95% confidence bounds for Fmarg2}
#' \item{F.est:}{Estimated survival function (computed at newT1, newT2)}
#' \item{F.est.lci:}{Lower 95% confidence bounds for F.est}
#' \item{F.est.uci:}{Upper 95% confidence bounds for F.est}
#' \item{CR:}{Estimated cross ratio (computed at T1, T2)}
#' \item{KT:}{Estimated Kendall\'s tau (computed at T1, T2)}
#' \item{CR.est:}{Estimated cross ratio (computed at newT1, newT2)}
#' \item{KT.est:}{Estimated Kendall\'s tau (computed at newT1, newT2)}
#' }
#' @seealso \code{\link[boot]{boot}}, \code{\link{mhazard-deprecated}}
#' @references
#' Prentice, R., Zhao, S. "Nonparametric estimation of the multivariate
#' survivor function: the multivariate Kaplan–Meier estimator", Lifetime
#' Data Analysis (2018) 24:3-27.
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
#' @examples
#' \donttest{x <- genClayton2(1000, 0, 1, 1, 2, 2)}
#' \donttest{x.km2 <- KM2(x$Y1, x$Y2, x$Delta1, x$Delta2)}
#' \donttest{x.km2.ci <- KM2(x$Y1, x$Y2, x$Delta1, x$Delta2, conf.int=TRUE)}
#' @name KM2-deprecated
#' @usage KM2(Y1, Y2, Delta1, Delta2, newT1=NULL, newT2=NULL,
#'	        estimator=c("dabrowska", "volterra", "prentice-cai"),
#'		conf.int=FALSE, R=1000, ...)
#' @keywords internal
NULL

#' @rdname mhazard-deprecated
#' @section \code{KM2}:
#' For \code{KM2}, use \code{\link{npSurv2}}.
#'
#' @export
KM2 <- function(Y1, Y2, Delta1, Delta2, newT1=NULL, newT2=NULL,
	        estimator=c("dabrowska", "volterra", "prentice-cai"),
		conf.int=FALSE, R=1000, ...){
    .Deprecated("npSurv2", package="mhazard")
    npSurv2(Y1, Y2, Delta1, Delta2, newT1, newT2, estimator, conf.int, R, ...)
}
