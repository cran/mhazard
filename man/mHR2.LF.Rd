% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mHR2.R
\name{mHR2.LF}
\alias{mHR2.LF}
\title{Bivariate regression survival function and marginal hazards estimation}
\usage{
mHR2.LF(
  mHR2.obj,
  X0_out,
  X1_out,
  T1_out,
  T2_out,
  confidence = c("none", "CI", "CB"),
  n.boot = 100
)
}
\arguments{
\item{mHR2.obj}{Output from the mHR2 function.}

\item{X0_out, X1_out}{Two possible sets of values for the covariates.
F and Lambda will be estimated at X=X0_out and X=X1_out.}

\item{T1_out, T2_out}{Vector of time points at which F and Lambda11
should be estimated. If confidence="CB", then both vectors must
have length 3.}

\item{confidence}{Type of confidence estimate to be computed.
Possible values include "none", "CI" (to compute confidence
intervals), and "CB" (to compute confidence bands). Defaults to
"none".}

\item{n.boot}{Number of bootstrap iterations for computing the
confidence intervals/bands. Defaults to 100. Ignored if
confidence="none".}
}
\value{
A list containing the following elements:
\describe{
\item{n10, n01:}{Total number of events for the first/second outcome}
\item{n11:}{Total number of double events}
\item{beta10, beta01, beta11:}{Regression coefficient estimates}
\item{lambda10, lambda01, lambda11:}{Baseline hazard estimates}
\item{Lambda11_out_Z0, Lambda11_out_Z1:}{Estimates of Lambda11 at
T1_out, T2_out for X=X0_out and X=X1_out}
\item{F_out_X0, F_out_X1:}{Estimates of F at T1_out, T2_out for
X=X0_out and X=X1_out}
\item{CI_Lambda11_X0.lb, CI_Lambda11_X0.ub:}{Lower and upper bounds
for Lambda11 at X=X0_out}
\item{CI_Lambda11_X1.lb, CI_Lambda11_X1.ub:}{Lower and upper bounds
for Lambda11 at X=X1_out}
\item{CI_F_X0.lb, CI_F_X0.ub:}{Lower and upper bounds for F at
X=X0_out}
\item{CI_F_X1.lb, CI_F_X1.ub:}{Lower and upper bounds for F at
X=X1_out}
\item{CB1_Lambda11_X0.lb, CB1_Lambda11_X0.ub, CB2_Lambda11_X0.lb,
CB2_Lambda11_X0.ub, CB3_Lambda11_X0.lb, CB3_Lambda11_X0.ub:}{Lower
and upper bounds for Lambda11 at X=X0_out, at three T1_out, T2_out
combinations}
\item{CB1_Lambda11_X1.lb, CB1_Lambda11_X1.ub, CB2_Lambda11_X1.lb,
CB2_Lambda11_X1.ub, CB3_Lambda11_X1.lb, CB3_Lambda11_X1.ub:}{Lower
and upper bounds for Lambda11 at X=X1_out, at three T1_out, T2_out
combinations}
\item{CB1_F_X0.lb, CB1_F_X0.ub, CB2_F_X0.lb, CB2_F_X0.ub, CB3_F_X0.lb,
CB3_F_X0.ub:}{Lower and upper bounds for F at X=X0_out, at three
T1_out, T2_out combinations}
\item{CB1_F_X1.lb, CB1_F_X1.ub, CB2_F_X1.lb, CB2_F_X1.ub, CB3_F_X1.lb,
CB3_F_X1.ub:}{Lower and upper bounds for F at X=X1_out, at three
T1_out, T2_out combinations}
}
}
\description{
Estimates the survival function F and the marginal hazards Lambda11
for a bivariate Cox regression model. F and Lambda11 are estimated
at two specified values of the covariates. If desired, (bootstrap)
confidence intervals or confidence bounds for F and Lambda11 may also
be computed.
}
\section{Details}{

If confidence="CI" or confidence="CB", then 95\% bootstrap confidence
bounds are computed by estimating the standard errors of F/Lambda11
based on n.boot bootstrap iterations. Currently confidence bounds
can only be computed at three specified T1out/T2out combinations
(meaning that T1out and T2out must both have length 3 if
confidence="CB"). No confidence measures will be returned if
confidence="none".
}

\examples{
x <- genClaytonReg(50, 2, 0.5, 1, 1, log(2), log(2), log(8/3), 2, 2)
x.mHR2 <- mHR2(x$Y1, x$Y2, x$Delta1, x$Delta2, x$X)
x.LF <- mHR2.LF(x.mHR2, 0, 1, c(0.25, 0.5, 1), c(0.25, 0.5, 1))
x.LF.CI <- mHR2.LF(x.mHR2, 0, 1, c(0.25, 0.5, 1),
c(0.25, 0.5, 1), confidence="CI")
x.LF.CB <- mHR2.LF(x.mHR2, 0, 1, c(0.25, 0.5, 1),
c(0.25, 0.5, 1), confidence="CB")
}
\references{
Prentice, R., Zhao, S. "The statistical analysis of multivariate
failure time data: A marginal modeling approach", CRC Press (2019).
Prentice, R., Zhao, S. "Regression models and multivariate life tables",
Journal of the American Statistical Association (2020) In press.
}
\seealso{
\code{\link{mHR2}}
}
