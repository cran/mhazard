#' Cox regression for a bivariate outcome
#'
#' Fits a semiparametric Cox regression model for a bivariate
#' outcome. This function computes the regression coefficients,
#' baseline hazards, and sandwich estimates of the standard
#' deviation of the regression coefficients. If desired, estimates
#' of the survival function F and marginal hazard rates Lambda11
#' can be computed using the mHR2.LF function.
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
#' @seealso \code{\link{mHR2.LF}}
#' @references
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
#' Prentice, R., Zhao, S. "Regression models and multivariate life tables",
#' Journal of the American Statistical Association (2021) 116(535):
#' 1330-1345. https://doi.org/10.1080/01621459.2020.1713792
#' @useDynLib mhazard
#' @importFrom Rcpp sourceCpp
#' @importFrom survival coxph basehaz
#' @importFrom rootSolve multiroot
#' @export
#' @examples
#' x <- genClaytonReg(1000, 2, 0.5, 1, 1, log(2), log(2), log(8/3), 2, 2)
#' x.mHR2 <- mHR2(x$Y1, x$Y2, x$Delta1, x$Delta2, x$X)
mHR2 <- function(Y1, Y2, Delta1, Delta2, X){
          X <- as.matrix(X)
          n <- length(Y1)
          p <- dim(X)[2]

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

          #####Parameter Estimates#####
          mod10 <- survival::coxph(survival::Surv(Y1, Delta1)~X, ties="breslow")
          beta10 <- mod10$coef
          mod01 <- survival::coxph(survival::Surv(Y2, Delta2)~X, ties="breslow")
          beta01 <- mod01$coef

          EE11 <- function(beta11, Y1, Y2, T1, T2, Delta1, Delta2,
                                X, dij_11){
                   return(calcEE11(beta11, Y1, Y2, T1, T2, Delta1, Delta2,
                                   X, dij_11))
          }
          beta11 <- rootSolve::multiroot(EE11, rep(0, p), Y1=Y1, Y2=Y2,
                                         T1=T1, T2=T2, Delta1=Delta1,
                                         Delta2=Delta2, X=X,
                                         dij_11=dij_11)$root


          #####variance estimates#####
          U10.k <- function(beta10, Y1, Delta1, X){
                   temp3 <- rep(NA, I)
                   X_ave <- matrix(rep(NA, p*I), nrow=p)
                   for (i in 1:I){
                        k.R <- which(Y1>=T1[i])
                        X.R <- X[k.R, , drop=FALSE]
                        temp1 <- as.vector(exp(X.R%*%beta10))
                        temp2 <- X.R*temp1
                        X_ave[, i] <- colSums(temp2)/sum(temp1)
                        temp3[i] <- dij_10[i]/sum(temp1)
                   }

                   U10.temp <- matrix(rep(0, p*n), nrow=p)
                   for (k in 1:n){
                        Ik <- sum(T1<=Y1[k])

                        temp4 <- rep(0, p)
                        if (Delta1[k]>0){
                            temp4 <- X[k, ]-X_ave[, Ik]
                        }

                        temp5 <- rep(0, p)
                        if (Ik>0){
                            temp5 <- rowSums((X[k, ]-X_ave[, 1:Ik, drop=FALSE])*temp3[1:Ik])
                        }

                        U10.temp[, k] <- temp4-temp5*as.numeric(exp(X[k, ]%*%beta10))
                   }
                   return(U10.temp)
          }

          U01.k <- function(beta01, Y2, Delta2, X){
                   temp3 <- rep(NA, J)
                   X_ave <- matrix(rep(NA, p*J), nrow=p)
                   for (j in 1:J){
                        k.R <- which(Y2>=T2[j])
                        X.R <- X[k.R, , drop=FALSE]
                        temp1 <- as.vector(exp(X.R%*%beta01))
                        temp2 <- X.R*temp1
                        X_ave[, j] <- colSums(temp2)/sum(temp1)
                        temp3[j] <- dij_01[j]/sum(temp1)
                   }

                   U01.temp <- matrix(rep(0, p*n), nrow=p)
                   for (k in 1:n){
                        Jk <- sum(T2<=Y2[k])

                        temp4 <- rep(0, p)
                        if (Delta2[k]>0){
                            temp4 <- X[k, ]-X_ave[, Jk]
                        }

                        temp5 <- rep(0, p)
                        if (Jk>0){
                            temp5 <- rowSums((X[k, ]-X_ave[, 1:Jk, drop=FALSE])*temp3[1:Jk])
                        }

                        U01.temp[, k] <- temp4-temp5*as.numeric(exp(X[k, ]%*%beta01))
                   }
                   return(U01.temp)
          }

          U11.k <- function(beta11, Y1, Y2, Delta1, Delta2, X){
                   temp3 <- matrix(rep(0, I*J), ncol=J)
                   X_ave <- Xtemp3 <- array(rep(0, p*I*J), dim=c(p, I, J))
                   for (i in 1:I){
                        for (j in 1:J){
                             if (dij_11[i+1, j+1]>0){
                                 k.R <- which(Y1>=T1[i] & Y2>=T2[j])
                                 X.R <- X[k.R, , drop=FALSE]
                                 temp1 <- as.vector(exp(X.R%*%beta11))
                                 temp2 <- X.R*temp1
                                 X_ave[, i, j] <- colSums(temp2)/sum(temp1)
                                 temp3[i, j] <- dij_11[i+1, j+1]/sum(temp1)
                                 Xtemp3[, i, j] <- X_ave[, i, j]*temp3[i, j]
                             }
                        }
                   }

                   U11.temp <- matrix(rep(0, p*n), nrow=p)
                   for (k in 1:n){
                        Ik <- sum(T1<=Y1[k])
                        Jk <- sum(T2<=Y2[k])

                        temp4 <- rep(0, p)
                        if (Delta1[k]*Delta2[k]>0){
                            temp4 <- X[k, ]-X_ave[, Ik, Jk]
                        }

                        temp5 <- temp6 <- rep(0, p)
                        if (Ik>0 & Jk>0){
                            temp5 <- X[k, ]*sum(temp3[1:Ik, 1:Jk])
                            temp6 <- apply(Xtemp3[, 1:Ik, 1:Jk, drop=FALSE], 1, sum)
                        }

                        U11.temp[, k] <- temp4-(temp5-temp6)*as.numeric(exp(X[k, ]%*%beta11))
                   }

                   return(U11.temp)
          }

          ##########Sandwich SD##########
          I10.temp <- array(rep(NA, p*p*I), dim=c(p, p, I))
          for (i in 1:I){
               t1 <- T1[i]
               k.R <- which(Y1>=t1)
               X.R <- X[k.R, , drop=FALSE]

               a <- as.vector(exp(X.R%*%beta10))
               b <- X.R*a
               temp0 <- sum(a)
               temp1 <- colSums(b)
               temp2 <- t(X.R)%*%b
               I10.temp[, , i] <- dij_10[i]*(temp0*temp2-outer(temp1, temp1))/(temp0^2)
          }
          I10 <- apply(I10.temp, c(1, 2), sum)

          I01.temp <- array(rep(NA, p*p*J), dim=c(p, p, J))
          for (j in 1:J){
               t2 <- T2[j]
               k.R <- which(Y2>=t2)
               X.R <- X[k.R, , drop=FALSE]

               a <- as.vector(exp(X.R%*%beta01))
               b <- X.R*a
               temp0 <- sum(a)
               temp1 <- colSums(b)
               temp2 <- t(X.R)%*%b
               I01.temp[, , j] <- dij_01[j]*(temp0*temp2-outer(temp1, temp1))/(temp0^2)
          }
          I01 <- apply(I01.temp, c(1, 2), sum)

          I11.temp <- array(rep(0, p*p*I*J), dim=c(p, p, I, J))
          for (i in 1:I){
               for (j in 1:J){
                    if (dij_11[i+1, j+1]>0){
                        t1 <- T1[i]
                        t2 <- T2[j]
                        k.R <- which(Y1>=t1 & Y2>=t2)

                        if (length(k.R)>0){
                            X.R <- X[k.R, , drop=FALSE]

                            a <- as.vector(exp(X.R%*%beta11))
                            b <- X.R*a
                            temp0 <- sum(a)
                            temp1 <- colSums(b)
                            temp2 <- t(X.R)%*%b
                            I11.temp[, , i, j] <- dij_11[i+1, j+1]*(temp0*temp2-outer(temp1, temp1))/(temp0^2)
                        }
                    }
               }
          }
          I11 <- apply(I11.temp, c(1, 2), sum)

          I_hat <- matrix(rep(0, 9*p*p), ncol=3*p)
          I_hat[1:p, 1:p] <- I10
          I_hat[(p+1):(2*p), (p+1):(2*p)] <- I01
          I_hat[(2*p+1):(3*p), (2*p+1):(3*p)] <- I11
          U_hat <- rbind(U10.k(beta10, Y1, Delta1, X), U01.k(beta01, Y2, Delta2, X), U11.k(beta11, Y1, Y2, Delta1, Delta2, X))
          A <- U_hat%*%t(U_hat)
          SD.sand <- sqrt(diag(solve(I_hat)%*%A%*%solve(I_hat)))

          ##########Baseline hazards##########
          cumhaz10 <- unique(survival::basehaz(mod10, centered=FALSE)$hazard)
          cumhaz01 <- unique(survival::basehaz(mod01, centered=FALSE)$hazard)
          lambda10 <- cumhaz10[2:(I+1)]-cumhaz10[1:I]
          lambda01 <- cumhaz01[2:(J+1)]-cumhaz01[1:J]
          lambda11 <- calc_lambda11(T1, T2, Y1, Y2, X, dij_11, beta11)

          return(list(Y1=Y1, Y2=Y2, Delta1=Delta1, Delta2=Delta2, X=X, n10=sum(Delta1), n01=sum(Delta2), n11=sum(Delta1*Delta2), beta10=beta10, beta01=beta01, beta11=beta11, SD.beta10=SD.sand[1:p], SD.beta01=SD.sand[(p+1):(2*p)], SD.beta11=SD.sand[(2*p+1):(3*p)], SD.beta10.cox=sqrt(mod10$var), SD.beta10.cox=sqrt(mod01$var),lambda10=lambda10, lambda01=lambda01, lambda11=lambda11))
}

#' Bivariate regression survival function and marginal hazards estimation
#'
#' Estimates the survival function F and the marginal hazards Lambda11
#' for a bivariate Cox regression model. F and Lambda11 are estimated
#' at two specified values of the covariates. If desired, (bootstrap)
#' confidence intervals or confidence bounds for F and Lambda11 may also
#' be computed.
#'
#' @param mHR2.obj Output from the mHR2 function.
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
#' @seealso \code{\link{mHR2}}
#' @references
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
#' Prentice, R., Zhao, S. "Regression models and multivariate life tables",
#' Journal of the American Statistical Association (2020) In press.
#' @useDynLib mhazard
#' @importFrom Rcpp sourceCpp
#' @importFrom survival coxph basehaz Surv
#' @importFrom rootSolve multiroot
#' @importFrom stats quantile sd
#' @export
#' @examples
#' x <- genClaytonReg(1000, 2, 0.5, 1, 1, log(2), log(2), log(8/3), 2, 2)
#' x.mHR2 <- mHR2(x$Y1, x$Y2, x$Delta1, x$Delta2, x$X)
#' x.LF <- mHR2.LF(x.mHR2, 0, 1, c(0.25, 0.5, 1), c(0.25, 0.5, 1))
#' \donttest{x.LF.CI <- mHR2.LF(x.mHR2, 0, 1, c(0.25, 0.5, 1),
#' c(0.25, 0.5, 1), confidence="CI")}
#' \donttest{x.LF.CB <- mHR2.LF(x.mHR2, 0, 1, c(0.25, 0.5, 1),
#' c(0.25, 0.5, 1), confidence="CB")}
mHR2.LF <- function(mHR2.obj, X0_out, X1_out, T1_out, T2_out,
                    confidence=c("none", "CI", "CB"), n.boot=100){
    Y1 <- mHR2.obj$Y1
    Y2 <- mHR2.obj$Y2
    Delta1 <- mHR2.obj$Delta1
    Delta2 <- mHR2.obj$Delta2
    X <- as.matrix(mHR2.obj$X)
    n.sub <- length(Y1)

    confidence <- match.arg(confidence)
    if (confidence=="CI") {
        result <- CI.boot.one(Y1, Y2, Delta1, Delta2, X, X0_out, X1_out, T1_out, T2_out)
        Lambda11.boot_Z0 <- Lambda11.boot_Z1 <- F.boot_Z0 <- F.boot_Z1 <- array(rep(NA, length(T1_out)*length(T2_out)*n.boot), dim=c(length(T1_out), length(T2_out), n.boot))

        for (i.boot in 1:n.boot){
            index <- sample(1:n.sub, replace=TRUE)
            Y1.boot <- Y1[index]
            Y2.boot <- Y2[index]
            Delta1.boot <- Delta1[index]
            Delta2.boot <- Delta2[index]
            X.boot <- as.matrix(X[index, ])
            result.boot <- CI.boot.one(Y1.boot, Y2.boot, Delta1.boot, Delta2.boot, X.boot, X0_out, X1_out, T1_out, T2_out)
            Lambda11.boot_Z0[, , i.boot] <- result.boot$Lambda11_out_Z0
            Lambda11.boot_Z1[, , i.boot] <- result.boot$Lambda11_out_Z1
            F.boot_Z0[, , i.boot] <- result.boot$F_out_Z0
            F.boot_Z1[, , i.boot] <- result.boot$F_out_Z1
        }

        SE_Lambda11_Z0 <- apply(Lambda11.boot_Z0, c(1, 2), function(x){sd(x, na.rm=TRUE)})
        SE_Lambda11_Z1 <- apply(Lambda11.boot_Z1, c(1, 2), function(x){sd(x, na.rm=TRUE)})
        SE_F_Z0 <- apply(F.boot_Z0, c(1, 2), function(x){sd(x, na.rm=TRUE)})
        SE_F_Z1 <- apply(F.boot_Z1, c(1, 2), function(x){sd(x, na.rm=TRUE)})

        CI_Lambda11_Z0.lb <- result$Lambda11_out_Z0-1.96*SE_Lambda11_Z0
        CI_Lambda11_Z0.ub <- result$Lambda11_out_Z0+1.96*SE_Lambda11_Z0
        CI_Lambda11_Z1.lb <- result$Lambda11_out_Z1-1.96*SE_Lambda11_Z1
        CI_Lambda11_Z1.ub <- result$Lambda11_out_Z1+1.96*SE_Lambda11_Z1

        CI_F_Z0.lb <- result$F_out_Z0-1.96*SE_F_Z0
        CI_F_Z0.ub <- result$F_out_Z0+1.96*SE_F_Z0
        CI_F_Z1.lb <- result$F_out_Z1-1.96*SE_F_Z1
        CI_F_Z1.ub <- result$F_out_Z1+1.96*SE_F_Z1

        return(list(n10=result$n10, n01=result$n01, n11=result$n11,
                    beta10=result$beta10, beta01=result$beta01, beta11=result$beta11,
                    lambda10=result$lambda10, lambda01=result$lambda01, lambda11=result$lambda11,
                    Lambda11_out_X0=result$Lambda11_out_Z0, Lambda11_out_X1=result$Lambda11_out_Z1,
                    CI_Lambda11_X0.lb=CI_Lambda11_Z0.lb, CI_Lambda11_X0.ub=CI_Lambda11_Z0.ub,
                    CI_Lambda11_X1.lb=CI_Lambda11_Z1.lb, CI_Lambda11_X1.ub=CI_Lambda11_Z1.ub,
                    F_out_X0=result$F_out_Z0, F_out_X1=result$F_out_Z1,
                    CI_F_X0.lb=CI_F_Z0.lb, CI_F_X0.ub=CI_F_Z0.ub,
                    CI_F_X1.lb=CI_F_Z1.lb, CI_F_X1.ub=CI_F_Z1.ub))
    }
    else if (confidence=="CB") {
        if (length(T1_out)!=3 | length(T2_out)!=3) {
            stop("T1_out, T2_out must have length 3 for confidence='CB'")
        }
        T1 <- sort(unique(Y1[Delta1==1]))
        T2 <- sort(unique(Y2[Delta2==1]))
        result <- CI.boot.one(Y1, Y2, Delta1, Delta2, X, X0_out, X1_out, T1, T2)

        Lambda11_Z0_hat <- result$Lambda11_out_Z0
        Lambda11_Z1_hat <- result$Lambda11_out_Z1
        F_Z0_hat <- result$F_out_Z0
        F_Z1_hat <- result$F_out_Z1

        ind.T1 <- unlist(lapply(T1_out, function(x){sum(x>=T1)}))
        ind.T2 <- unlist(lapply(T2_out, function(x){sum(x>=T2)}))

        mW_Lambda11.boot_Z0 <- mW_Lambda11.boot_Z1 <- mW_F.boot_Z0 <- mW_F.boot_Z1 <- matrix(rep(NA, 3*n.boot), ncol=3)

        for (i.boot in 1:n.boot){
            index <- sample(1:n.sub, replace=TRUE)
            Y1.boot <- Y1[index]
            Y2.boot <- Y2[index]
            Delta1.boot <- Delta1[index]
            Delta2.boot <- Delta2[index]
            X.boot <- as.matrix(X[index, ])
            result.boot <- CI.boot.one(Y1.boot, Y2.boot, Delta1.boot, Delta2.boot, X.boot, X0_out, X1_out, T1, T2, calc_rij=TRUE)


            W_Lambda11.boot_Z0 <- sqrt(n.sub)*(result.boot$Lambda11_out_Z0-Lambda11_Z0_hat)*I(result.boot$rij_out>0)
            W_Lambda11.boot_Z1 <- sqrt(n.sub)*(result.boot$Lambda11_out_Z1-Lambda11_Z1_hat)*I(result.boot$rij_out>0)
            W_F.boot_Z0 <- sqrt(n.sub)*(result.boot$F_out_Z0-F_Z0_hat)*I(result.boot$rij_out>0)
            W_F.boot_Z1 <- sqrt(n.sub)*(result.boot$F_out_Z1-F_Z1_hat)*I(result.boot$rij_out>0)

            r1.boot.Z0 <- sum((Y1.boot[X.boot==0]>=T1_out[1]) & (Y2.boot[X.boot==0]>T2_out[1]))
            r2.boot.Z0 <- sum((Y1.boot[X.boot==0]>=T1_out[2]) & (Y2.boot[X.boot==0]>T2_out[2]))
            r3.boot.Z0 <- sum((Y1.boot[X.boot==0]>=T1_out[3]) & (Y2.boot[X.boot==0]>T2_out[3]))
            r1.boot.Z1 <- sum((Y1.boot[X.boot==1]>=T1_out[1]) & (Y2.boot[X.boot==1]>T2_out[1]))
            r2.boot.Z1 <- sum((Y1.boot[X.boot==1]>=T1_out[2]) & (Y2.boot[X.boot==1]>T2_out[2]))
            r3.boot.Z1 <- sum((Y1.boot[X.boot==1]>=T1_out[3]) & (Y2.boot[X.boot==1]>T2_out[3]))

            if (r1.boot.Z0>0){
                mW_Lambda11.boot_Z0[i.boot, 1] <- max(abs(W_Lambda11.boot_Z0[1:ind.T1[1], 1:ind.T2[1]]))
                mW_F.boot_Z0[i.boot, 1] <- max(abs(W_F.boot_Z0[1:ind.T1[1], 1:ind.T2[1]]))
            }
            if (r2.boot.Z0>0){
                mW_Lambda11.boot_Z0[i.boot, 2] <- max(abs(W_Lambda11.boot_Z0[1:ind.T1[2], 1:ind.T2[2]]))
                mW_F.boot_Z0[i.boot, 2] <- max(abs(W_F.boot_Z0[1:ind.T1[2], 1:ind.T2[2]]))
            }
            if (r3.boot.Z0>0){
                mW_Lambda11.boot_Z0[i.boot, 3] <- max(abs(W_Lambda11.boot_Z0[1:ind.T1[3], 1:ind.T2[3]]))
                mW_F.boot_Z0[i.boot, 3] <- max(abs(W_F.boot_Z0[1:ind.T1[3], 1:ind.T2[3]]))
            }

            if (r1.boot.Z1>0){
                mW_Lambda11.boot_Z1[i.boot, 1] <- max(abs(W_Lambda11.boot_Z1[1:ind.T1[1], 1:ind.T2[1]]))
                mW_F.boot_Z1[i.boot, 1] <- max(abs(W_F.boot_Z1[1:ind.T1[1], 1:ind.T2[1]]))
            }
            if (r2.boot.Z1>0){
                mW_Lambda11.boot_Z1[i.boot, 2] <- max(abs(W_Lambda11.boot_Z1[1:ind.T1[2], 1:ind.T2[2]]))
                mW_F.boot_Z1[i.boot, 2] <- max(abs(W_F.boot_Z1[1:ind.T1[2], 1:ind.T2[2]]))
            }
            if (r3.boot.Z1>0){
                mW_Lambda11.boot_Z1[i.boot, 3] <- max(abs(W_Lambda11.boot_Z1[1:ind.T1[3], 1:ind.T2[3]]))
                mW_F.boot_Z1[i.boot, 3] <- max(abs(W_F.boot_Z1[1:ind.T1[3], 1:ind.T2[3]]))
            }
        }

        C.Lambda11_Z0 <- apply(mW_Lambda11.boot_Z0, 2, function(x){quantile(x, 0.95, na.rm=TRUE)})
        C.Lambda11_Z1 <- apply(mW_Lambda11.boot_Z1, 2, function(x){quantile(x, 0.95, na.rm=TRUE)})
        C.F_Z0 <- apply(mW_F.boot_Z0, 2, function(x){quantile(x, 0.95, na.rm=TRUE)})
        C.F_Z1 <- apply(mW_F.boot_Z1, 2, function(x){quantile(x, 0.95, na.rm=TRUE)})

        CB1_Lambda11_Z0.lb <- Lambda11_Z0_hat[1:ind.T1[1], 1:ind.T2[1]]-C.Lambda11_Z0[1]/sqrt(n.sub)
        CB1_Lambda11_Z0.ub <- Lambda11_Z0_hat[1:ind.T1[1], 1:ind.T2[1]]+C.Lambda11_Z0[1]/sqrt(n.sub)
        CB1_Lambda11_Z1.lb <- Lambda11_Z1_hat[1:ind.T1[1], 1:ind.T2[1]]-C.Lambda11_Z1[1]/sqrt(n.sub)
        CB1_Lambda11_Z1.ub <- Lambda11_Z1_hat[1:ind.T1[1], 1:ind.T2[1]]+C.Lambda11_Z1[1]/sqrt(n.sub)

        CB1_F_Z0.lb <- F_Z0_hat[1:ind.T1[1], 1:ind.T2[1]]-C.F_Z0[1]/sqrt(n.sub)
        CB1_F_Z0.ub <- F_Z0_hat[1:ind.T1[1], 1:ind.T2[1]]+C.F_Z0[1]/sqrt(n.sub)
        CB1_F_Z1.lb <- F_Z1_hat[1:ind.T1[1], 1:ind.T2[1]]-C.F_Z1[1]/sqrt(n.sub)
        CB1_F_Z1.ub <- F_Z1_hat[1:ind.T1[1], 1:ind.T2[1]]+C.F_Z1[1]/sqrt(n.sub)

        CB2_Lambda11_Z0.lb <- Lambda11_Z0_hat[1:ind.T1[2], 1:ind.T2[2]]-C.Lambda11_Z0[2]/sqrt(n.sub)
        CB2_Lambda11_Z0.ub <- Lambda11_Z0_hat[1:ind.T1[2], 1:ind.T2[2]]+C.Lambda11_Z0[2]/sqrt(n.sub)
        CB2_Lambda11_Z1.lb <- Lambda11_Z1_hat[1:ind.T1[2], 1:ind.T2[2]]-C.Lambda11_Z1[2]/sqrt(n.sub)
        CB2_Lambda11_Z1.ub <- Lambda11_Z1_hat[1:ind.T1[2], 1:ind.T2[2]]+C.Lambda11_Z1[2]/sqrt(n.sub)

        CB2_F_Z0.lb <- F_Z0_hat[1:ind.T1[2], 1:ind.T2[2]]-C.F_Z0[2]/sqrt(n.sub)
        CB2_F_Z0.ub <- F_Z0_hat[1:ind.T1[2], 1:ind.T2[2]]+C.F_Z0[2]/sqrt(n.sub)
        CB2_F_Z1.lb <- F_Z1_hat[1:ind.T1[2], 1:ind.T2[2]]-C.F_Z1[2]/sqrt(n.sub)
        CB2_F_Z1.ub <- F_Z1_hat[1:ind.T1[2], 1:ind.T2[2]]+C.F_Z1[2]/sqrt(n.sub)

        CB3_Lambda11_Z0.lb <- Lambda11_Z0_hat[1:ind.T1[3], 1:ind.T2[3]]-C.Lambda11_Z0[3]/sqrt(n.sub)
        CB3_Lambda11_Z0.ub <- Lambda11_Z0_hat[1:ind.T1[3], 1:ind.T2[3]]+C.Lambda11_Z0[3]/sqrt(n.sub)
        CB3_Lambda11_Z1.lb <- Lambda11_Z1_hat[1:ind.T1[3], 1:ind.T2[3]]-C.Lambda11_Z1[3]/sqrt(n.sub)
        CB3_Lambda11_Z1.ub <- Lambda11_Z1_hat[1:ind.T1[3], 1:ind.T2[3]]+C.Lambda11_Z1[3]/sqrt(n.sub)

        CB3_F_Z0.lb <- F_Z0_hat[1:ind.T1[3], 1:ind.T2[3]]-C.F_Z0[3]/sqrt(n.sub)
        CB3_F_Z0.ub <- F_Z0_hat[1:ind.T1[3], 1:ind.T2[3]]+C.F_Z0[3]/sqrt(n.sub)
        CB3_F_Z1.lb <- F_Z1_hat[1:ind.T1[3], 1:ind.T2[3]]-C.F_Z1[3]/sqrt(n.sub)
        CB3_F_Z1.ub <- F_Z1_hat[1:ind.T1[3], 1:ind.T2[3]]+C.F_Z1[3]/sqrt(n.sub)

        return(list(n=n.sub, n10=result$n10, n01=result$n01, n11=result$n11, T1=T2, T2=T2,
                    beta10=result$beta10, beta01=result$beta01, beta11=result$beta11, lambda10=result$lambda10, lambda01=result$lambda01, lambda11=result$lambda11,
                    Lambda11_out_X0=result$Lambda11_out_Z0, Lambda11_out_X1=result$Lambda11_out_Z1,
                    CB1_Lambda11_X0.lb=CB1_Lambda11_Z0.lb, CB1_Lambda11_X0.ub=CB1_Lambda11_Z0.ub,
                    CB1_Lambda11_X1.lb=CB1_Lambda11_Z1.lb, CB1_Lambda11_X1.ub=CB1_Lambda11_Z1.ub,
                    CB2_Lambda11_X0.lb=CB2_Lambda11_Z0.lb, CB2_Lambda11_X0.ub=CB2_Lambda11_Z0.ub,
                    CB2_Lambda11_X1.lb=CB2_Lambda11_Z1.lb, CB2_Lambda11_X1.ub=CB2_Lambda11_Z1.ub,
                    CB3_Lambda11_X0.lb=CB3_Lambda11_Z0.lb, CB3_Lambda11_X0.ub=CB3_Lambda11_Z0.ub,
                    CB3_Lambda11_X1.lb=CB3_Lambda11_Z1.lb, CB3_Lambda11_X1.ub=CB3_Lambda11_Z1.ub,
                    F_out_X0=result$F_out_Z0, F_out_X1=result$F_out_Z1,
                    CB1_F_X0.lb=CB1_F_Z0.lb, CB1_F_X0.ub=CB1_F_Z0.ub,
                    CB1_F_X1.lb=CB1_F_Z1.lb, CB1_F_X1.ub=CB1_F_Z1.ub,
                    CB2_F_X0.lb=CB2_F_Z0.lb, CB2_F_X0.ub=CB2_F_Z0.ub,
                    CB2_F_X1.lb=CB2_F_Z1.lb, CB2_F_X1.ub=CB2_F_Z1.ub,
                    CB3_F_X0.lb=CB3_F_Z0.lb, CB3_F_X0.ub=CB3_F_Z0.ub,
                    CB3_F_X1.lb=CB3_F_Z1.lb, CB3_F_X1.ub=CB3_F_Z1.ub))
    }
    else {
        result <- CI.boot.one(Y1, Y2, Delta1, Delta2, X, X0_out, X1_out, T1_out, T2_out)
        return(list(n10=result$n10, n01=result$n01, n11=result$n11,
                    beta10=result$beta10, beta01=result$beta01, beta11=result$beta11,
                    lambda10=result$lambda10, lambda01=result$lambda01, lambda11=result$lambda11,
                    Lambda11_out_X0=result$Lambda11_out_Z0, Lambda11_out_X1=result$Lambda11_out_Z1,
                    F_out_X0=result$F_out_Z0, F_out_X1=result$F_out_Z1))
    }
}

do.F <- function(T1, T2, beta10, beta01, beta11, lambda10, lambda01, lambda11, x){
        I <- length(T1)
        J <- length(T2)

        lambda10.x <- rep(NA, I+1)
        lambda01.x <- rep(NA, J+1)
        lambda11.x <- matrix(rep(NA, (I+1)*(J+1)), ncol=(J+1))
        lambda10.x[1] <- 0
        lambda01.x[1] <- 0
        lambda10.x[2:(I+1)] <- lambda10*exp(as.numeric(x%*%beta10))
        lambda01.x[2:(J+1)] <- lambda01*exp(as.numeric(x%*%beta01))
        lambda11.x <- lambda11*exp(as.numeric(x%*%beta11))

        F.x <- calcF(lambda10.x, lambda01.x, lambda11.x)
        return(F.x)
}

CI.boot.one <- function(Y1, Y2, Delta1, Delta2, X, X0_out, X1_out, T1_out, T2_out, calc_rij=FALSE){
               X <- as.matrix(X)
               n <- length(Y1)
               p <- dim(X)[2]

               T1 <- unique(sort(Y1[Delta1==1]))
               T2 <- unique(sort(Y2[Delta2==1]))
               I <- length(T1)
               J <- length(T2)

               #####counts#####
               if (calc_rij) {
                   junk <- calc_drij(Y1, Y2, c(0, T1), c(0, T2), Delta1, Delta2)
                   dij_11 <- junk[,,2]
                   rij <- junk[,,1]
               }
               else {
                   dij_11 <- matrix(rep(0, (I+1)*(J+1)), ncol=J+1)
                   id11 <- which(Delta1*Delta2==1)
                   for (k in id11){
                       Ik <- sum(T1<=Y1[k])
                       Jk <- sum(T2<=Y2[k])
                       dij_11[Ik+1, Jk+1] <- dij_11[Ik+1, Jk+1]+1
                   }
               }

               #####Parameter Estimates#####
               mod10 <- survival::coxph(survival::Surv(Y1, Delta1)~X, ties="breslow")
               beta10 <- mod10$coef
               mod01 <- survival::coxph(survival::Surv(Y2, Delta2)~X, ties="breslow")
               beta01 <- mod01$coef

               EE11 <- function(beta11, Y1, Y2, T1, T2, Delta1, Delta2,
                                X, dij_11){
                   return(calcEE11(beta11, Y1, Y2, T1, T2, Delta1, Delta2,
                                   X, dij_11))
               }
               beta11 <- rootSolve::multiroot(EE11, rep(0, p), Y1=Y1, Y2=Y2,
                                              T1=T1, T2=T2, Delta1=Delta1,
                                              Delta2=Delta2, X=X,
                                              dij_11=dij_11)$root

               cumhaz10 <- unique(survival::basehaz(mod10,
                                                    centered=FALSE)$hazard)
               cumhaz01 <- unique(survival::basehaz(mod01,
                                                    centered=FALSE)$hazard)
               lambda10 <- cumhaz10[2:(I+1)]-cumhaz10[1:I]
               lambda01 <- cumhaz01[2:(J+1)]-cumhaz01[1:J]
               lambda11 <- calc_lambda11(T1, T2, Y1, Y2, X, dij_11, beta11)

               lambda11_Z0 <- lambda11*exp(as.numeric(X0_out%*%beta11))
               lambda10_Z0 <- lambda10*exp(as.numeric(X0_out%*%beta10))
               lambda01_Z0 <- lambda01*exp(as.numeric(X0_out%*%beta01))

               lambda11_Z1 <- lambda11*exp(as.numeric(X1_out%*%beta11))
               lambda10_Z1 <- lambda10*exp(as.numeric(X1_out%*%beta10))
               lambda01_Z1 <- lambda01*exp(as.numeric(X1_out%*%beta01))

               Lambda11_Z0 <- t(apply(apply(lambda11_Z0, 2, cumsum), 1, cumsum))
               Lambda11_Z1 <- t(apply(apply(lambda11_Z1, 2, cumsum), 1, cumsum))

               F_Z0 <- do.F(T1, T2, beta10, beta01, beta11, lambda10, lambda01, lambda11, X0_out)
               F_Z1 <- do.F(T1, T2, beta10, beta01, beta11, lambda10, lambda01, lambda11, X1_out)

               index.i <- unlist(lapply(T1_out, function(x){sum(T1<=x)}))
               index.j <- unlist(lapply(T2_out, function(x){sum(T2<=x)}))

               if (calc_rij) {
                   rij_out <- rij[index.i+1, index.j+1]
               }
               Lambda11_out_Z0 <- Lambda11_Z0[index.i+1, index.j+1]
               Lambda11_out_Z1 <- Lambda11_Z1[index.i+1, index.j+1]
               F_out_Z0 <- F_Z0[index.i+1, index.j+1]
               F_out_Z1 <- F_Z1[index.i+1, index.j+1]

               if (calc_rij) {
                   return(list(n10=sum(Delta1), n01=sum(Delta2), n11=sum(Delta1*Delta2),
                               T1=T1, T2=T2,
                               beta10=beta10, beta01=beta01, beta11=beta11,
                               lambda10=lambda10, lambda01=lambda01, lambda11=lambda11,
                               rij_out=rij_out,
                               Lambda11_out_Z0=Lambda11_out_Z0, Lambda11_out_Z1=Lambda11_out_Z1,
                               F_out_Z0=F_out_Z0, F_out_Z1=F_out_Z1))
               }
               else {
                   return(list(n10=sum(Delta1), n01=sum(Delta2), n11=sum(Delta1*Delta2),
                               T1=T1, T2=T2,
                               beta10=beta10, beta01=beta01, beta11=beta11,
                               lambda10=lambda10, lambda01=lambda01, lambda11=lambda11,
                               Lambda11_out_Z0=Lambda11_out_Z0, Lambda11_out_Z1=Lambda11_out_Z1,
                               F_out_Z0=F_out_Z0, F_out_Z1=F_out_Z1))
               }
}
