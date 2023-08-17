#' Cox regression for a bivariate outcome with time-varying covariates
#'
#' Fits a semiparametric Cox regression model for a bivariate
#' outcome with time-varying covariates. Currently only the regression
#' coefficients are computed.
#'
#' @param Y1,Y2 Vectors of event times (continuous).
#' @param Delta1,Delta2 Vectors of censoring indicators (1=event,
#' 0=censored).
#' @param ids Vector of ID numbers. It is used to map the values of
#' the time-varying covariates back to the original
#' Y1/Y2/Delta1/Delta2 values. See Details.
#' @param X Matrix of covariates (continuous or binary). See Details
#' for the proper format of this matrix.
#' @section Details:
#' X must be a matrix with at least four columns. The first column
#' contains the ID numbers. Each ID number in this column must map
#' to a unique element of the ids vector. The second and third columns
#' consists of time points for T1 and T2, respectively. They specify
#' the time points at which the covariates take on the specified
#' value(s). The remaining columns represent the values
#' of the covariates on the specified time interval. For example,
#' if we define
#' X.tv <- matrix(c(1001, 1001, 0, 0, 0, 5, 1, 2), nrow=2)
#' then, for the observation with ID number 1001, then when T1=0, the
#' time-varying covariate has a value of 1 on when T2 is in [0,5) and a
#' value of 2 when T2 is in [5,Inf). Note that the values of the
#' time-varying covariates must be specified for when T1=0 (or T2=0)
#' in order to compute beta10 and beta01. If a value of a covariate is
#' constant when T1=0 or T2=0, that covariate will be dropped when
#' computing beta10 or beta01.
#'
#' Support for time-varying covariates is experimental and has not
#' been tested extensively. Use this function at your own risk.
#' @return A list containing the following elements:
#' \describe{
#' \item{beta10, beta01, beta11:}{Regression coefficient estimates}
#' }
#' @seealso \code{\link{mHR2}}
#' @references
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
#' Prentice, R., Zhao, S. "Regression models and multivariate life tables",
#' Journal of the American Statistical Association (2021) 116(535):
#' 1330-1345. https://doi.org/10.1080/01621459.2020.1713792
#' @useDynLib mhazard
#' @importFrom rootSolve multiroot
#' @export
#' @examples
#' x <- genClaytonReg(50, 2, 0.5, 1, 1, 0, log(2), 0, 5, 5)
#' x.tv <- tvc.example(x$Y1, x$Y2, x$Delta1, x$Delta2, x$X)
#' x.mHR2 <- mHR2.tvc(x$Y1, x$Y2, x$Delta1, x$Delta2,
#' x.tv$ids, x.tv$X.tv)
mHR2.tvc <- function(Y1, Y2, Delta1, Delta2, ids, X) {
    n <- length(Y1)
    if (length(Y2)!=n | length(Delta1)!=n | length(Delta2)!=n |
        length(ids)!=n) {
        stop("Y1, Y2, Delta1, Delta2, ids must have the same length")
    }
    if (sum(duplicated(ids))>0) {
        stop("ids must be unique")
    }
    if (ncol(X)<4) {
        stop("X must have at least 4 columns")
    }
    if (!all(sort(ids)==sort(unique(X[,1])))) {
            stop("Invalid ids for X")
    }
    X <- as.matrix(X)
    X <- X[order(X[,1], X[,2], X[,3]),]
    p <- ncol(X)-3

          T1 <- unique(sort(Y1[Delta1==1]))
          T2 <- unique(sort(Y2[Delta2==1]))
          I <- length(T1)
          J <- length(T2)

          #####counts#####
          dij_11 <- matrix(rep(0, (I+1)*(J+1)), ncol=J+1)
          dij_10 <- rep(0, I)
          dij_01 <- rep(0, J)

          id11 <- which(Delta1*Delta2==1)
          for (k in id11){
               Ik <- sum(T1<=Y1[k])
               Jk <- sum(T2<=Y2[k])
               dij_11[Ik+1, Jk+1] <- dij_11[Ik+1, Jk+1]+1
          }

          for (i in 1:I){
               dij_10[i] <- sum(Y1==T1[i] & Delta1==1)
          }
          for (j in 1:J){
               dij_01[j] <- sum(Y2==T2[j] & Delta2==1)
          }

          #####Marginal Cox Models: time-varying#####
          EE10 <- function(beta10){
                  EE10.temp <- matrix(rep(NA, p.cur*I), ncol=I)
                  for (i in 1:I){
                       t1 <- T1[i]
                       i.D <- which(Y1==t1 & Delta1==1)
                       i.D.ids <- ids[i.D]
                       i.R <- which(Y1>=t1)
                       i.R.ids <- ids[i.R]

                       Xt.D <- matrix(nrow=length(i.D), ncol=p.cur)
                       Xt.R <- matrix(nrow=length(i.R), ncol=p.cur)
                       for (k in 1:length(i.D)) {
                           cur.tv <- X.cur[X.cur[,1]==i.D.ids[k],]
                           cur.t <- max(which(cur.tv[,2]<=t1))
                           if (cur.t<0) {
                               stop("Time values must be nonnegative")
                           }
                           Xt.D[k,] <- cur.tv[cur.t,-(1:3)]
                       }
                       for (k in 1:length(i.R)) {
                           cur.tv <- X.cur[X.cur[,1]==i.R.ids[k],]
                           cur.t <- max(which(cur.tv[,2]<=t1))
                           if (cur.t<0) {
                               stop("Time values must be nonnegative")
                           }
                           Xt.R[k,] <- cur.tv[cur.t,-(1:3)]
                       }

                       EE10.temp[, i] <- colSums(Xt.D)-dij_10[i]*colSums(Xt.R*as.vector(exp(Xt.R%*%beta10)))/colSums(exp(Xt.R%*%beta10))
                  }
                  return(rowSums(EE10.temp))
          }
          X.cur <- X[X[,3]==0,]
          if (nrow(X.cur)==0) {
              stop("X(t,0) undefined")
          }
          X.keep <- apply(X.cur, 2, function(x) length(unique(x)))>1
          X.keep[3] <- TRUE
          X.cur <- X.cur[,X.keep]
          p.cur <- ncol(X.cur)-3

          find.root10 <- multiroot(EE10, rep(0, p.cur))
          beta10 <- find.root10$root

          EE01 <- function(beta01){
                  EE01.temp <- matrix(rep(NA, p.cur*J), ncol=J)
                  for (j in 1:J){
                       t2 <- T2[j]
                       i.D <- which(Y2==t2 & Delta2==1)
                       i.D.ids <- ids[i.D]
                       i.R <- which(Y2>=t2)
                       i.R.ids <- ids[i.R]

                       Xt.D <- matrix(nrow=length(i.D), ncol=p.cur)
                       Xt.R <- matrix(nrow=length(i.R), ncol=p.cur)
                       for (k in 1:length(i.D)) {
                           cur.tv <- X.cur[X.cur[,1]==i.D.ids[k],]
                           cur.t <- max(which(cur.tv[,3]<=t2))
                           if (cur.t<0) {
                               stop("Time values must be nonnegative")
                           }
                           Xt.D[k,] <- cur.tv[cur.t,-(1:3)]
                       }
                       for (k in 1:length(i.R)) {
                           cur.tv <- X.cur[X.cur[,1]==i.R.ids[k],]
                           cur.t <- max(which(cur.tv[,3]<=t2))
                           if (cur.t<0) {
                               stop("Time values must be nonnegative")
                           }
                           Xt.R[k,] <- cur.tv[cur.t,-(1:3)]
                       }

                       EE01.temp[, j] <- colSums(Xt.D)-dij_01[j]*colSums(Xt.R*as.vector(exp(Xt.R%*%beta01)))/colSums(exp(Xt.R%*%beta01))
               }
                  return(rowSums(EE01.temp))
          }
          X.cur <- X[X[,2]==0,]
          if (nrow(X.cur)==0) {
              stop("X(0,t) undefined")
          }
          X.keep <- apply(X.cur, 2, function(x) length(unique(x)))>1
          X.keep[2] <- TRUE
          X.cur <- X.cur[,X.keep]
          p.cur <- ncol(X.cur)-3

          find.root01 <- multiroot(EE01, rep(0, p.cur))
          beta01 <- find.root01$root

          #####Double Failure Hazard#####
          EE11 <- function(beta11){
              EE11.temp <- array(rep(0, p*I*J), dim=c(p, I, J))
                  for (i in 1:I){
                       for (j in 1:J){
                            if (dij_11[i+1, j+1]>0){
                                t1 <- T1[i]
                                t2 <- T2[j]
                                i.D <- which(Y1==t1 & Y2==t2 & Delta1==1 & Delta2==1)
                                i.D.ids <- ids[i.D]
                                i.R <- which(Y1>=t1 & Y2>=t2)
                                i.R.ids <- ids[i.R]
                                Xt.D <- matrix(nrow=length(i.D), ncol=p)
                                Xt.R <- matrix(nrow=length(i.R), ncol=p)

                                for (k in 1:length(i.D)) {
                                    cur.tv <- X[X[,1]==i.D.ids[k],]
                                    cur.t <- max(which(cur.tv[,2]<=t1&
                                                       cur.tv[,3]<=t2))
                                    if (cur.t<0) {
                                        stop("Time values must be nonnegative")
                                    }
                                    Xt.D[k,] <- cur.tv[cur.t,-(1:3)]
                                }
                                for (k in 1:length(i.R)) {
                                    cur.tv <- X[X[,1]==i.R.ids[k],]
                                    cur.t <- max(which(cur.tv[,2]<=t1&
                                                       cur.tv[,3]<=t2))
                                    if (cur.t<0) {
                                        stop("Time values must be nonnegative")
                                    }
                                    Xt.R[k,] <- cur.tv[cur.t,-(1:3)]
                                }


                                EE11.temp[, i, j] <- colSums(Xt.D)-dij_11[i+1, j+1]*colSums(Xt.R*as.vector(exp(Xt.R%*%beta11)))/colSums(exp(Xt.R%*%beta11))
                            }
                       }
                  }
                  EE11 <- apply(EE11.temp, 1, sum)
                  return(EE11)
          }
    find.root <- try(multiroot(EE11, rep(0, p)), silent=TRUE)
    beta11 <- find.root$root

    return(list(beta10=beta10, beta01=beta01, beta11=beta11))
}
