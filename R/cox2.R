#' Cox regression for a bivariate outcome
#'
#' Fits a semiparametric Cox regression model for a bivariate
#' outcome. This function computes the regression coefficients,
#' baseline hazards, and sandwich estimates of the standard
#' deviation of the regression coefficients. If desired, estimates
#' of the survival function F and marginal hazard rates Lambda11
#' can be computed using the cox2.LF function.
#'
#' @param Y1,Y2 Vectors of event times (continuous).
#' @param Delta1,Delta2 Vectors of censoring indicators (1=event,
#' 0=censored).
#' @param X Matrix of covariates (continuous or binary).
#' @return A list containing the following elements:
#' \describe{
#' \item{Y1, Y2:}{Original vectors of event times}
#' \item{Delta1, Delta2:}{Original vectors of censoring indicators}
#' \item{X:}{Original covariate matrix}
#' \item{n10, n01:}{Total number of events for the first/second outcome}
#' \item{n11:}{Total number of double events}
#' \item{beta10, beta01, beta11:}{Regression coefficient estimates}
#' \item{lambda10, lambda01, lambda11:}{Baseline hazard estimates}
#' \item{SD.beta10, SD.beta01, SD.beta11:}{Sandwich estimates of the
#' standard deviation of the regression coefficients}
#' \item{SD.beta10.cox, SD.beta01.cox:}{Standard deviation estimates
#' for the regression coefficients based on a univariate Cox model}
#' }
#' @seealso \code{\link{mhazard-deprecated}}
#' @references
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
#' Prentice, R., Zhao, S. "Regression models and multivariate life tables",
#' Journal of the American Statistical Association (2020) In press.
#' @examples
#' \donttest{x <- genClaytonReg(1000, 2, 0.5, 1, 1, log(2), log(2),
#' log(8/3), 2, 2)}
#' \donttest{x.cox2 <- cox2(x$Y1, x$Y2, x$Delta1, x$Delta2, x$X)}
#' @name cox2-deprecated
#' @usage cox2(Y1, Y2, Delta1, Delta2, X)
#' @keywords internal
NULL

#' @rdname mhazard-deprecated
#' @section \code{cox2}:
#' For \code{cox2}, use \code{\link{mHR2}}.
#'
#' @export
cox2 <- function(Y1, Y2, Delta1, Delta2, X){
    .Deprecated("mHR2", package="mhazard")
    mHR2(Y1, Y2, Delta1, Delta2, X)
}

#' Bivariate regression survival function and marginal hazards estimation
#'
#' Estimates the survival function F and the marginal hazards Lambda11
#' for a bivariate Cox regression model. F and Lambda11 are estimated
#' at two specified values of the covariates. If desired, (bootstrap)
#' confidence intervals or confidence bounds for F and Lambda11 may also
#' be computed.
#'
#' @param cox2.obj Output from the cox2 function.
#' @param X0_out,X1_out Two possible sets of values for the covariates.
#' F and Lambda will be estimated at X=X0_out and X=X1_out.
#' @param T1_out,T2_out Vector of time points at which F and Lambda11
#' should be estimated. If confidence="CB", then both vectors must
#' have length 3.
#' @param confidence Type of confidence estimate to be computed.
#' Possible values include "none", "CI" (to compute confidence
#' intervals), and "CB" (to compute confidence bands). Defaults to
#' "none".
#' @param n.boot Number of bootstrap iterations for computing the
#' confidence intervals/bands. Defaults to 100. Ignored if
#' confidence="none".
#' @section Details:
#' If confidence="CI" or confidence="CB", then 95% bootstrap confidence
#' bounds are computed by estimating the standard errors of F/Lambda11
#' based on n.boot bootstrap iterations. Currently confidence bounds
#' can only be computed at three specified T1out/T2out combinations
#' (meaning that T1out and T2out must both have length 3 if
#' confidence="CB"). No confidence measures will be returned if
#' confidence="none".
#' @return A list containing the following elements:
#' \describe{
#' \item{n10, n01:}{Total number of events for the first/second outcome}
#' \item{n11:}{Total number of double events}
#' \item{beta10, beta01, beta11:}{Regression coefficient estimates}
#' \item{lambda10, lambda01, lambda11:}{Baseline hazard estimates}
#' \item{Lambda11_out_Z0, Lambda11_out_Z1:}{Estimates of Lambda11 at
#' T1_out, T2_out for X=X0_out and X=X1_out}
#' \item{F_out_X0, F_out_X1:}{Estimates of F at T1_out, T2_out for
#' X=X0_out and X=X1_out}
#' \item{CI_Lambda11_X0.lb, CI_Lambda11_X0.ub:}{Lower and upper bounds
#' for Lambda11 at X=X0_out}
#' \item{CI_Lambda11_X1.lb, CI_Lambda11_X1.ub:}{Lower and upper bounds
#' for Lambda11 at X=X1_out}
#' \item{CI_F_X0.lb, CI_F_X0.ub:}{Lower and upper bounds for F at
#' X=X0_out}
#' \item{CI_F_X1.lb, CI_F_X1.ub:}{Lower and upper bounds for F at
#' X=X1_out}
#' \item{CB1_Lambda11_X0.lb, CB1_Lambda11_X0.ub, CB2_Lambda11_X0.lb,
#' CB2_Lambda11_X0.ub, CB3_Lambda11_X0.lb, CB3_Lambda11_X0.ub:}{Lower
#' and upper bounds for Lambda11 at X=X0_out, at three T1_out, T2_out
#' combinations}
#' \item{CB1_Lambda11_X1.lb, CB1_Lambda11_X1.ub, CB2_Lambda11_X1.lb,
#' CB2_Lambda11_X1.ub, CB3_Lambda11_X1.lb, CB3_Lambda11_X1.ub:}{Lower
#' and upper bounds for Lambda11 at X=X1_out, at three T1_out, T2_out
#' combinations}
#' \item{CB1_F_X0.lb, CB1_F_X0.ub, CB2_F_X0.lb, CB2_F_X0.ub, CB3_F_X0.lb,
#' CB3_F_X0.ub:}{Lower and upper bounds for F at X=X0_out, at three
#' T1_out, T2_out combinations}
#' \item{CB1_F_X1.lb, CB1_F_X1.ub, CB2_F_X1.lb, CB2_F_X1.ub, CB3_F_X1.lb,
#' CB3_F_X1.ub:}{Lower and upper bounds for F at X=X1_out, at three
#' T1_out, T2_out combinations}
#' }
#' @seealso \code{\link{mhazard-deprecated}}
#' @references
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
#' Prentice, R., Zhao, S. "Regression models and multivariate life tables",
#' Journal of the American Statistical Association (2020) In press.
#' @examples
#' \donttest{x <- genClaytonReg(1000, 2, 0.5, 1, 1, log(2), log(2),
#' log(8/3), 2, 2)}
#' \donttest{x.cox2 <- cox2(x$Y1, x$Y2, x$Delta1, x$Delta2, x$X)}
#' \donttest{x.LF <- cox2.LF(x.cox2, 0, 1, c(0.25, 0.5, 1), c(0.25, 0.5, 1))}
#' \donttest{x.LF.CI <- cox2.LF(x.cox2, 0, 1, c(0.25, 0.5, 1),
#' c(0.25, 0.5, 1), confidence="CI")}
#' \donttest{x.LF.CB <- cox2.LF(x.cox2, 0, 1, c(0.25, 0.5, 1),
#' c(0.25, 0.5, 1), confidence="CB")}
#' @name cox2.LF-deprecated
#' @usage cox2.LF(cox2.obj, X0_out, X1_out, T1_out, T2_out,
#'                confidence=c("none", "CI", "CB"), n.boot=100)
#' @keywords internal
NULL

#' @rdname mhazard-deprecated
#' @section \code{cox2.LF}:
#' For \code{cox2.LF}, use \code{\link{mHR2.LF}}.
#'
#' @export
cox2.LF <- function(cox2.obj, X0_out, X1_out, T1_out, T2_out,
                    confidence=c("none", "CI", "CB"), n.boot=100){
    .Deprecated("mHR2.LF", package="mhazard")
    mHR2.LF(cox2.obj, X0_out, X1_out, T1_out, T2_out, confidence, n.boot)
}
