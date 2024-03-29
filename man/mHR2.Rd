% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mHR2.R
\name{mHR2}
\alias{mHR2}
\title{Cox regression for a bivariate outcome}
\usage{
mHR2(Y1, Y2, Delta1, Delta2, X)
}
\arguments{
\item{Y1, Y2}{Vectors of event times (continuous).}

\item{Delta1, Delta2}{Vectors of censoring indicators (1=event,
0=censored).}

\item{X}{Matrix of covariates (continuous or binary).}
}
\value{
A list containing the following elements:
\describe{
\item{Y1, Y2:}{Original vectors of event times}
\item{Delta1, Delta2:}{Original vectors of censoring indicators}
\item{X:}{Original covariate matrix}
\item{n10, n01:}{Total number of events for the first/second outcome}
\item{n11:}{Total number of double events}
\item{beta10, beta01, beta11:}{Regression coefficient estimates}
\item{lambda10, lambda01, lambda11:}{Baseline hazard estimates}
\item{SD.beta10, SD.beta01, SD.beta11:}{Sandwich estimates of the
standard deviation of the regression coefficients}
\item{SD.beta10.cox, SD.beta01.cox:}{Standard deviation estimates
for the regression coefficients based on a univariate Cox model}
}
}
\description{
Fits a semiparametric Cox regression model for a bivariate
outcome. This function computes the regression coefficients,
baseline hazards, and sandwich estimates of the standard
deviation of the regression coefficients. If desired, estimates
of the survival function F and marginal hazard rates Lambda11
can be computed using the mHR2.LF function.
}
\examples{
x <- genClaytonReg(1000, 2, 0.5, 1, 1, log(2), log(2), log(8/3), 2, 2)
x.mHR2 <- mHR2(x$Y1, x$Y2, x$Delta1, x$Delta2, x$X)
}
\references{
Prentice, R., Zhao, S. "The statistical analysis of multivariate
failure time data: A marginal modeling approach", CRC Press (2019).
Prentice, R., Zhao, S. "Regression models and multivariate life tables",
Journal of the American Statistical Association (2021) 116(535):
1330-1345. https://doi.org/10.1080/01621459.2020.1713792
}
\seealso{
\code{\link{mHR2.LF}}
}
