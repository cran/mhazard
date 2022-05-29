#' Generates survival data from a bivariate Clayton-Oakes model
#'
#' Generates simulated survival data from a bivariate Clayton-Oakes model,
#' which can be used to create example data for bivariate survival
#' function estimation. The marginal distributions are exponential with
#' given rate parameters. The joint distribution is defined using a Clayton
#' copula. The censoring times are also exponentially distributed with
#' given rate parameters.
#'
#' @param n Sample size for the simulated data set.
#' @param theta Parameter for the Clayton copula. Must be -1 or larger.
#' @param lambda10,lambda01 Rate parameters for the (marginal) exponential
#' distributions.
#' @param lambdaC1,lambdaC2 Rate parameters for the censoring times. No
#' censoring occurs if this parameter is equal to 0.
#' @section Details:
#' This function simulates data with the following survival function:
#' F(t1,t2) = \[F(t1,0)^(-theta) + F(0,t2)^(-theta) - 1\]^(-1/theta)
#' (The survival function is defined to be equal to 0 if this
#' quantity is negative.) The marginal survival functions F(t1,0) and
#' F(0,t2) are exponentially distributed with rate parameters lambda10 and
#' lambda01, respectively. After generating survival times Y1 and Y2 (of
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
#' @references
#' Clayton, D. "Model for association in bivariate life tables and its
#' application in epidemiological studies of familial tendency in chronic
#' disease incidence.", Biometrika (1978) 65:141-151.
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
#' @importFrom stats runif rexp
#' @export
#' @examples
#' x <- genClayton2(1000, 0, 1, 1, 2, 2)
genClayton2 <- function(n, theta, lambda10, lambda01, lambdaC1, lambdaC2){
    if (theta<(-1)) {
        stop("theta must be >= -1")
    }
    if (min(c(lambda01, lambda10))<=0) {
        stop("lambda10, lambda01 must be positive")
    }
    if (min(c(lambdaC1, lambdaC2))<0) {
        stop("lambdaC1, lambdaC2 must be nonnegative")
    }

               U1 <- runif(n, 0, 1)
               U2 <- runif(n, 0, 1)

               eta <- theta   #note: not true for every simulation

               T2 <- -log(1-U2)/lambda01
               T1 <- rep(NA, n)

               if (eta==0){
                   T1 <- -log(1-U1)/lambda10
               }else{
                   T1 <- log(1-(1-U2)^(-eta)+(1-U1)^(-eta/(1+eta))*(1-U2)^(-eta))/eta/lambda10
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

               return(data.frame(Y1=Y1, Y2=Y2, Delta1=Delta1, Delta2=Delta2))
}
