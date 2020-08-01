#' Generates survival data from a trivariate Clayton-Oakes model
#'
#' Generates simulated survival data from a trivariate Clayton-Oakes model,
#' which can be used to create example data for trivariate survival
#' function estimation. The marginal distributions are exponential with
#' rate parameter 1. The joint distribution is defined using a Clayton
#' copula. The censoring times are also exponentially distributed with
#' given rate parameters.
#'
#' @param n Sample size for the simulated data set.
#' @param theta Parameter for the Clayton copula. Must be -1 or larger.
#' @param lambdaC1,lambdaC2,lambdaC3 Rate parameters for the censoring
#' times. No censoring occurs if this paramter is equal to 0.
#' @section Details:
#' This function simulates data with the following survival function:
#' F(t1,t2,t3) = \[F(t1,0,0)^(-theta) + F(0,t2,0)^(-theta) +
#' F(0,0,t3) - 2\]^(-1/theta)
#' (The survival function is defined to be equal to 0 if this
#' quantity is negative.) The marginal survival functions F(t1,0,0),
#' F(0,t2,0), and F(0,0,t3) are exponentially distributed with rate
#' parameter 1. After generating survival times Y1, Y2, and Y3 (of
#' length n) under this distribution, censoring times C1, C2, and C3 (also
#' of length n) are generated. C1/C2/C3 are generated under an exponential
#' distribution with rate parameters lambdaC1, lambdaC2, lambdaC3,
#' respectively. If C1\[i\]<Y1\[i\] for a given observation i, then
#' observation i is considered to be censored (i.e., Delta1\[i\]=0).
#' Delta2 and Delta3 are defined in a similar manner. If lambdaC1,
#' lambdaC2, and/or lambdaC3 is equal to 0, then the corresponding
#' variable is uncensored (meaning that Delta\[i\]=1 for all i).
#' @return A data frame containing the following elements:
#' \describe{
#' \item{Y1, Y2, Y3:}{Survival times for the simulated data}
#' \item{Delta1, Delta2, Delta3:}{Censoring indicators for the simulated
#' data}
#' }
#' @references
#' Clayton, D. "Model for association in bivariate life tables and its
#' application in epidemiological studies of familial tendency in chronic
#' disease incidence.", Biometrika (1978) 65:141-151.
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
#' @importFrom stats runif rexp
#' @export
#' @examples
#' x <- genClayton3(200, 0, 0.5, 0.5, 0.5)
genClayton3 <- function(n, theta, lambdaC1, lambdaC2, lambdaC3) {
    if (theta<(-1)) {
        stop("theta must be >= -1")
    }
    if (min(c(lambdaC1, lambdaC2, lambdaC3))<0) {
        stop("lambdaC1, lambdaC2, lambdaC3 must be nonnegative")
    }

  u1 <- runif(n)
  u2 <- runif(n)
  u3 <- runif(n)
  t1 <- -1*log(u1)
  if (theta==0) {
    t2 <- -1*log(u2)
    t3 <- -1*log(u3)
  }
  else {
    t2 <- theta^(-1) * log(1 - u1^(-1*theta) * (1 - u2^(-1*theta*
       	    (1+theta)^(-1))))
    t3 <- theta^(-1) * log(1 - u1^(-1*theta)*u2^(-1*theta*(1+theta)^(-1))*
       	    (1 - u3^(-1*theta*(1+2*theta)^(-1))))
  }
    failt <- cbind(t1, t2, t3)
    centt <- matrix(nrow=n, ncol=3)
    if (lambdaC1 != 0) {
        centt[,1] <- rexp(n, lambdaC1)
    }
    else {
        centt[,1] <- max(t1)+1
    }
    if (lambdaC2 != 0) {
        centt[,2] <- rexp(n, lambdaC2)
    }
    else {
        centt[,2] <- max(t2)+1
    }
    if (lambdaC3 != 0) {
        centt[,3] <- rexp(n, lambdaC3)
    }
    else {
        centt[,3] <- max(t3)+1
    }

    Delta <- 1*(failt<centt)
    Y1 <- apply(cbind(failt[,1], centt[,1]), 1, min)
    Y2 <- apply(cbind(failt[,2], centt[,2]), 1, min)
    Y3 <- apply(cbind(failt[,3], centt[,3]), 1, min)
    return(data.frame(Y1=Y1, Y2=Y2, Y3=Y3, Delta1=Delta[,1], Delta2=Delta[,2],
                      Delta3=Delta[,3]))
}
