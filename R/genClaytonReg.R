#' Generates regression data from a bivariate Clayton-Oakes model
#'
#' Generates simulated survival data from a bivariate Clayton-Oakes model
#' where the hazard depends on a binary coefficient X. This can be used to
#' create example data for bivariate Cox regression. The marginal
#' distributions are exponential with given rate parameters. The joint
#' distribution is defined using a Clayton copula. The censoring times are
#' also exponentially distributed with given rate parameters.
#'
#' @param n Sample size for the simulated data set.
#' @param theta Parameter for the Clayton copula. Must be -1 or larger.
#' @param Xp Probability that the covariate is equal to 1. Must satisfy
#' 0<Xp<1.
#' @param lambda10,lambda01 Rate parameters for the (marginal) exponential
#' distributions when X=0.
#' @param b10,b01,b11 Regression coefficient values.
#' @param lambdaC1,lambdaC2 Rate parameters for the censoring times. No
#' censoring occurs if this parameter is equal to 0.
#' @section Details:
#' This function simulates data with the following survival function:
#' F(t1,t2) = \[F(t1,0)^(-eta) + F(0,t2)^(-eta) - 1\]^(-1/eta)
#' (The survival function is defined to be equal to 0 if this
#' quantity is negative.) Here eta=theta*exp(X*b11). The marginal
#' survival functions F(t1,0) and F(0,t2) are exponentially distributed
#' with rate parameters lambda10*exp(X*b10) and lambda01*exp(X*b01),
#' respectively. After generating survival times Y1 and Y2 (of
#' length n) under this distribution, censoring times C1 and C2 (also of
#' length n) are generated. C1/C2 are generated under an exponential
#' distribution with rate parameters lambdaC1 and lambdaC2. If
#' C1\[i\]<Y1\[i\] for a given observation i, then observation i is
#' considered to be censored (i.e., Delta1\[i\]=0). Delta2 is defined in
#' a similar manner. If lambdaC1 or lambdaC2 is equal to 0, then the
#' corresponding variable is uncensored (meaning that Delta\[i\]=1 for
#' all i).
#' @return A data frame containing the following elements:
#' \describe{
#' \item{Y1, Y2:}{Survival times for the simulated data}
#' \item{Delta1, Delta2:}{Censoring indicators for the simulated data}
#' }
#' \item{X}{Covariate matrix (of dimension n x 1).}
#' @references
#' Clayton, D. "Model for association in bivariate life tables and its
#' application in epidemiological studies of familial tendency in chronic
#' disease incidence.", Biometrika (1978) 65:141-151.
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
#' @importFrom stats runif rexp rbinom
#' @export
#' @examples
#' x <- genClaytonReg(1000, 2, 0.5, 1, 1, log(2), log(2), log(8/3), 2, 2)
genClaytonReg <- function(n, theta, Xp, lambda10, lambda01, b10, b01, b11, lambdaC1, lambdaC2){
    if (theta<(-1)) {
        stop("theta must be >= -1")
    }
    if (Xp<=0 | Xp>=1) {
        stop("Xp must satisfy 0<Xp<1")
    }
    if (min(c(lambda01, lambda10))<=0) {
        stop("lambda10, lambda01 must be positive")
    }
    if (min(c(lambdaC1, lambdaC2))<0) {
        stop("lambdaC1, lambdaC2 must be nonnegative")
    }
               X <- rbinom(n, 1, Xp)
               X <- as.matrix(X)
               U1 <- runif(n, 0, 1)
               U2 <- runif(n, 0, 1)

               eta <- theta*exp(-X%*%b11)   #note: not true for every simulation

               T2 <- -log(1-U2)/lambda01/as.vector(exp(X%*%b01))
               T1 <- rep(NA, n)
               for (i in 1:n){
                    if (eta[i]==0){
                        T1[i] <- -log(1-U1[i])/lambda10/exp(X[i, , drop=FALSE]%*%b10)
                    }else{
                        T1[i] <- log(1-(1-U2[i])^(-eta[i])+(1-U1[i])^(-eta[i]/(1+eta[i]))*(1-U2[i])^(-eta[i]))/eta[i]/lambda10/exp(X[i, , drop=FALSE]%*%b10)
                    }
               }

               if (lambdaC1!=0){
                   C1 <- rexp(n, lambdaC1)
               }else{
                   C1 <- max(T1)+1  #this generates no censoring on T1
               }

               if (lambdaC2!=0){
                   C2 <- rexp(n, lambdaC2)
               }else{
                   C2 <- max(T2)+1  #this generates no censoring on T2
               }

               Y1 <- pmin(T1, C1)
               Y2 <- pmin(T2, C2)
               Delta1 <- (T1<=C1)
               Delta2 <- (T2<=C2)

    return(data.frame(Y1=Y1, Y2=Y2, Delta1=Delta1, Delta2=Delta2,
                      X=X))
}

#' Creates an example of a matrix of time-varying covariates
#'
#' Given a set of (non-time-varying) covariates, creates a simple example
#' of a matrix of time-varying covariates that can be used as input data
#' for the mHR2.tvc function.
#'
#' @param Y1,Y2 Vectors of event times (continuous).
#' @param Delta1,Delta2 Vectors of censoring indicators (1=event,
#' 0=censored).
#' @param X Matrix of covariates (continuous or binary).
#' @section Details:
#' For each (non-time-varying) covariate in X, two time-varying
#' covariates are created. The first time-varying covariate is equal to
#' X*log(T1), and the second is equal to X*log(T2). (If T=0, then the
#' time-varying covariate is set to be 0.) A vector of ID numbers and a
#' matrix of time-varying covariates are created in a format that can be
#' passed to the mHR2.tvc function.
#' @return A list containing the following elements:
#' \describe{
#' \item{ids:}{A vector of ids}
#' \item{X.tv:}{Time-varying covariate matrix}
#' }
#' @seealso \code{\link{mHR2}}, \code{\link{genClaytonReg}}
#' @export
#' @examples
#' x <- genClaytonReg(250, 2, 0.5, 1, 1, 0, log(2), 0, 5, 5)
#' x.tv <- tvc.example(x$Y1, x$Y2, x$Delta1, x$Delta2, x$X)
tvc.example <- function(Y1, Y2, Delta1, Delta2, X)
{
    T1 <- c(0,unique(sort(Y1[Delta1==1])))
    T2 <- c(0,unique(sort(Y2[Delta2==1])))

    n <- length(Y1)
    X <- as.matrix(X)
    np <- ncol(X)
    I <- length(T1)
    J <- length(T2)

    X.tv <- matrix(nrow=n*I*J, ncol=3*(np+1))

    ids <- seq(from=1001, to=(1000+n), by=1)

    for (k in 1:n) {
        for (i in 1:I) {
            for (j in 1:J) {
                X.tv[(k-1)*I*J+(i-1)*J+j,1] <- 1000+k
                X.tv[(k-1)*I*J+(i-1)*J+j,2] <- T1[i]
                X.tv[(k-1)*I*J+(i-1)*J+j,3] <- T2[j]
                X.tv[(k-1)*I*J+(i-1)*J+j,(4:(np+3))] <- X[k,]
                if (T1[i]>0) {
                    X.tv[(k-1)*I*J+(i-1)*J+j,((np+4):(2*np+3))] <-
                        X[k,]*log(T1[i])
                }
                else {
                    X.tv[(k-1)*I*J+(i-1)*J+j,((np+4):(2*np+3))] <- 0
                }
                if (T2[j]>0) {
                    X.tv[(k-1)*I*J+(i-1)*J+j,((2*np+4):(3*np+3))] <-
                        X[k,]*log(T2[j])
                }
                else {
                    X.tv[(k-1)*I*J+(i-1)*J+j,((2*np+4):(3*np+3))] <- 0
                }
            }
        }
    }

    return(list(ids=ids, X.tv=X.tv))
}
