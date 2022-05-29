#' Nonparametric estimates of the survival function for trivariate failure
#' time data
#'
#' Computes the survival function for a trivariate failure time data. The
#' survival function for trivariate failure time data is analogous to the
#' Kaplan-Meier estimator for a univariate failure time data and Dabrowska
#' estimator for bivariate failure time data. Optionally (bootstrap)
#' confidence intervals for the survival function may also be computed.
#'
#' @param Y1,Y2,Y3 Vectors of event times (continuous).
#' @param Delta1,Delta2,Delta3 Vectors of censoring indicators (1=event,
#' 0=censored).
#' @param newT1,newT2,newT3 Optional vectors of times at which to estimate
#' the survival function (which do not need to be subsets of Y1/Y2/Y3).
#' Defaults to the unique values in Y1/Y2/Y3 if not specified.
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
#' \item{T3:}{Unique values of Y3 at which Fhat was computed}
#' \item{Fhat:}{Estimated survival function (computed at T1, T2, T3)}
#' \item{Fhat.lci:}{Lower 95% confidence bounds for Fhat}
#' \item{Fhat.uci:}{Upper 95% confidence bounds for Fhat}
#' \item{Fmarg1.est:}{Estimated marginal survival function for T11
#' (computed at newT1)}
#' \item{Fmarg1.lci:}{Lower 95% confidence bounds for Fmarg1}
#' \item{Fmarg1.uci:}{Upper 95% confidence bounds for Fmarg1}
#' \item{Fmarg2.est:}{Estimated marginal survival function for T2
#' (computed at newT2)}
#' \item{Fmarg2.lci:}{Lower 95% confidence bounds for Fmarg2}
#' \item{Fmarg2.uci:}{Upper 95% confidence bounds for Fmarg2}
#' \item{Fmarg3.est:}{Estimated marginal survival function for T3
#' (computed at newT3)}
#' \item{Fmarg3.lci:}{Lower 95% confidence bounds for Fmarg3}
#' \item{Fmarg3.uci:}{Upper 95% confidence bounds for Fmarg3}
#' \item{F.est:}{Estimated survival function (computed at newT1, newT2,
#' newT3)}
#' \item{F.est.lci:}{Lower 95% confidence bounds for F.est}
#' \item{F.est.uci:}{Upper 95% confidence bounds for F.est}
#' \item{C110:}{Pairwise marginal cross ratio estimator C110 (computed at
#' newT1, newT2, newT3)}
#' \item{C101:}{Pairwise marginal cross ratio estimator C101 (computed at
#' newT1, newT2, newT3)}
#' \item{C011:}{Pairwise marginal cross ratio estimator C011 (computed at
#' new T1, newT2, newT3)}
#' \item{C111:}{Trivariate dependency estimator C111 (computed at newT1,
#' newT2, newT3)}
#' }
#' @seealso \code{\link[boot]{boot}}, \code{\link{mhazard-deprecated}}
#' @references
#' Prentice, R., Zhao, S. "Nonparametric estimation of the multivariate
#' survivor function: the multivariate Kaplan–Meier estimator", Lifetime
#' Data Analysis (2018) 24:3-27.
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
#' @examples
#' \donttest{x <- genClayton3(200, 0, 0.5, 0.5, 0.5)}
#' \donttest{x.km3 <- KM3(x$Y1, x$Y2, x$Y3, x$Delta1, x$Delta2, x$Delta3)}
#' \donttest{x.km3.ci <- KM3(x$Y1, x$Y2, x$Y3, x$Delta1, x$Delta2,
#' x$Delta3, conf.int=TRUE, R=500)}
#' @name KM3-deprecated
#' @usage KM3(Y1, Y2, Y3, Delta1, Delta2, Delta3, newT1=NULL, newT2=NULL,
#'            newT3=NULL, conf.int=FALSE, R=1000, ...)
#' @keywords internal
NULL

#' @rdname mhazard-deprecated
#' @section \code{KM3}:
#' For \code{KM3}, use \code{\link{npSurv3}}.
#'
#' @export
KM3 <- function(Y1, Y2, Y3, Delta1, Delta2, Delta3, newT1=NULL, newT2=NULL,
       		newT3=NULL, conf.int=FALSE, R=1000, ...){
    .Deprecated("npSurv3", package="mhazard")
    npSurv3(Y1, Y2, Y3, Delta1, Delta2, Delta3, newT1, newT2, newT3, conf.int,
            R, ...)
}
