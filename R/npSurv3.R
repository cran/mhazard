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
#' @seealso \code{\link[boot]{boot}}
#' @references
#' Prentice, R., Zhao, S. "Nonparametric estimation of the multivariate
#' survivor function: the multivariate Kaplan–Meier estimator", Lifetime
#' Data Analysis (2018) 24:3-27.
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
#' pp. 120-123.
#' @useDynLib mhazard
#' @importFrom Rcpp sourceCpp
#' @importFrom boot boot
#' @export
#' @examples
#' x <- genClayton3(200, 0, 0.5, 0.5, 0.5)
#' x.npSurv3 <- npSurv3(x$Y1, x$Y2, x$Y3, x$Delta1, x$Delta2, x$Delta3)
#' \donttest{x.npSurv3.ci <- npSurv3(x$Y1, x$Y2, x$Y3, x$Delta1, x$Delta2,
#' x$Delta3, conf.int=TRUE, R=500)}
npSurv3 <- function(Y1, Y2, Y3, Delta1, Delta2, Delta3, newT1=NULL, newT2=NULL,
                    newT3=NULL, conf.int=FALSE, R=1000, ...){
    if (min(c(Y1, Y2, Y3))<0) {
        stop("Y1, Y2, Y3 must be nonnegative")
    }
    if (sum(!(c(Delta1, Delta2, Delta3) %in% c(0, 1)))>0) {
        stop("Delta1, Delta2, Delta3 must be 0/1 vectors")
    }
          I <- sum(Delta1)
          J <- sum(Delta2)
          K <- sum(Delta3)

          T1 <- unique(sort(Y1[Delta1==1]))
          T2 <- unique(sort(Y2[Delta2==1]))
          T3 <- unique(sort(Y3[Delta3==1]))

	  if (is.null(newT1)) {
	    newT1 <- T1
	  }
	  if (is.null(newT2)) {
	    newT2 <- T2
	  }
	  if (is.null(newT3)) {
	    newT3 <- T3
	  }

    if (min(c(newT1, newT2, newT3))<0) {
        stop("newT1, newT2, newT3 must be nonnegative")
    }

	  Fhat <- triHaz(Y1, Y2, Y3, T1, T2, T3, Delta1, Delta2, Delta3)

	  if (conf.int) {
	    Fhat.boot <- boot::boot(cbind(Y1, Y2, Y3, Delta1, Delta2, Delta3),
                                    KM3.boot, R=R, ...)
	    Fhat.ci <- getBootCI(Fhat.boot$t)
#	    Fhat.ci <- apply(Fhat.boot$t, 2, quantile, probs=c(0.025, 0.975),
#	    	     	     type=8)
	    Fhat.lci <- array(Fhat.ci[1,], dim=dim(Fhat))
	    Fhat.uci <- array(Fhat.ci[2,], dim=dim(Fhat))
	  }

          ##########summary results at
	  index.i <- sapply(newT1, function(x1, x2) {sum(x2<=x1)},
	                    x2=T1)
	  index.j <- sapply(newT2, function(x1, x2) {sum(x2<=x1)},
	                    x2=T2)
	  index.k <- sapply(newT3, function(x1, x2) {sum(x2<=x1)},
	                    x2=T3)

          ###########estimated marginal survival probability#########
          Fmarg1 <- Fhat[index.i+1, 1, 1]
          Fmarg2 <- Fhat[1, index.j+1, 1]
          Fmarg3 <- Fhat[1, 1, index.k+1]

          ##########estimate average cross-ratio##########
          Fhat_d1 <- Fhat[1:I, , , drop=FALSE]-Fhat[2:(I+1), , , drop=FALSE]
          Fhat_d2 <- Fhat[, 1:J, , drop=FALSE]-Fhat[, 2:(J+1), , drop=FALSE]
          Fhat_d3 <- Fhat[, , 1:K, drop=FALSE]-Fhat[, , 2:(K+1), drop=FALSE]

          Fhat_d12 <- Fhat_d1[, 1:J, , drop=FALSE]-Fhat_d1[, 2:(J+1), , drop=FALSE]
          Fhat_d13 <- Fhat_d1[, , 1:K, drop=FALSE]-Fhat_d1[, , 2:(K+1), drop=FALSE]
          Fhat_d23 <- Fhat_d2[, , 1:K, drop=FALSE]-Fhat_d2[, , 2:(K+1), drop=FALSE]

          Fhat_d123 <- Fhat_d12[, , 1:K, drop=FALSE]-Fhat_d12[, , 2:(K+1), drop=FALSE]

          temp <- calcTemp3(Fhat)
          Lambda_111_0 <- 1-(Fhat_d1[, 1:J, 1:K, drop=FALSE]+Fhat_d2[1:I, , 1:K, drop=FALSE]+Fhat_d3[1:I, 1:J, , drop=FALSE]-Fhat_d12[, , 1:K, drop=FALSE]-Fhat_d13[, 1:J, , drop=FALSE]-Fhat_d23[1:I, , , drop=FALSE])/Fhat[1:I, 1:J, 1:K, drop=FALSE]-temp

          C110 <- C101 <- C011 <- matrix(rep(NA, 9), ncol=3)
          for (m in 1:3){
               for (n in 1:3){
                    C110[m, n] <- sum(Fhat_d12[1:index.i[m], 1:index.j[n], 1])/sum(Fhat_d1[1:index.i[m], 1:index.j[n], 1]*Fhat_d2[1:index.i[m], 1:index.j[n], 1]/Fhat[1:index.i[m], 1:index.j[n], 1])
                    C101[m, n] <- sum(Fhat_d13[1:index.i[m], 1, 1:index.k[n]])/sum(Fhat_d1[1:index.i[m], 1, 1:index.k[n]]*Fhat_d3[1:index.i[m], 1, 1:index.k[n]]/Fhat[1:index.i[m], 1, 1:index.k[n]])
                    C011[m, n] <- sum(Fhat_d23[1, 1:index.j[m], 1:index.k[n]])/sum(Fhat_d2[1, 1:index.j[m], 1:index.k[n]]*Fhat_d3[1, 1:index.j[m], 1:index.k[n]]/Fhat[1, 1:index.j[m], 1:index.k[n]])
               }
          }

          ##########estimate trivariate dependency##########
          C111 <- array(rep(NA, 27), dim=c(3, 3, 3))
          for (m in 1:3){
               for (n in 1:3){
                    for (l in 1:3){
                         C111[m, n, l] <- sum(Fhat_d123[1:index.i[m], 1:index.j[n], 1:index.k[l]], na.rm=TRUE)/sum(Fhat[1:index.i[m], 1:index.j[n], 1:index.k[l]]*Lambda_111_0[1:index.i[m], 1:index.j[n], 1:index.k[l]], na.rm=TRUE)
                    }
              }
          }

	  if (!conf.int) {
          return(list(Fhat=Fhat, T1=T1, T2=T2, T3=T3, Fmarg1=Fmarg1, Fmarg2=Fmarg2, Fmarg3=Fmarg3, C110=C110, C101=C101, C011=C011, C111=C111, Fhat_est=Fhat[index.i+1, index.j+1, index.k+1]))
	  }
	  else {
	     Fmarg1.lci <- Fhat.lci[index.i+1, 1, 1]
             Fmarg2.lci <- Fhat.lci[1, index.j+1, 1]
             Fmarg3.lci <- Fhat.lci[1, 1, index.k+1]
	     Fmarg1.uci <- Fhat.uci[index.i+1, 1, 1]
             Fmarg2.uci <- Fhat.uci[1, index.j+1, 1]
             Fmarg3.uci <- Fhat.uci[1, 1, index.k+1]
	     return(list(Fhat=Fhat, Fhat.lci=Fhat.lci, Fhat.uci=Fhat.uci, T1=T1, T2=T2, T3=T3, Fmarg1=Fmarg1, Fmarg2=Fmarg2, Fmarg3=Fmarg3, Fmarg1.lci=Fmarg1.lci, Fmarg1.uci=Fmarg1.uci, Fmarg2.lci=Fmarg2.lci, Fmarg2.uci=Fmarg2.uci, Fmarg3.lci=Fmarg3.lci, Fmarg3.uci=Fmarg3.uci, C110=C110, C101=C101, C011=C011, C111=C111, Fhat_est=Fhat[index.i+1, index.j+1, index.k+1], Fhat_est.lci=Fhat.lci[index.i+1, index.j+1, index.k+1], Fhat_est.uci=Fhat.uci[index.i+1, index.j+1, index.k+1]))
	  }
}

KM3.boot <- function(x, ndx) {
  Y1 <- x[ndx,1]
  Y2 <- x[ndx,2]
  Y3 <- x[ndx,3]
  Delta1 <- x[ndx,4]
  Delta2 <- x[ndx,5]
  Delta3 <- x[ndx,6]

  newT1 <-  unique(sort(x[x[,4]==1,1]))
  newT2 <-  unique(sort(x[x[,5]==1,2]))
  newT3 <-  unique(sort(x[x[,6]==1,3]))
  T1 <- unique(sort(Y1[Delta1==1]))
  T2 <- unique(sort(Y2[Delta2==1]))
  T3 <- unique(sort(Y3[Delta3==1]))

  Fhat <- triHaz(Y1, Y2, Y3, T1, T2, T3, Delta1, Delta2, Delta3)

  index.i <- sapply(newT1, function(x1, x2) {sum(x2<=x1)}, x2=T1)
  index.j <- sapply(newT2, function(x1, x2) {sum(x2<=x1)}, x2=T2)
  index.k <- sapply(newT3, function(x1, x2) {sum(x2<=x1)}, x2=T3)

  return(as.vector(Fhat[c(1,index.i+1),c(1,index.j+1),c(1,index.k+1)]))
}
