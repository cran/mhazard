#' Uses a 3D perspective plot to visualize a nonparametric bivariate
#' survival function
#'
#' Plots a 3D perspective plot of an estimated nonparametric bivariate
#' survival function. This function is a wrapper for the persp3D function
#' from the plot3D package with default parameters chosen to make the
#' data easier to visualize.
#'
#' @param npSurv2.obj Output of the npSurv2 function.
#' @param col Color palette to be used for the plot. Defaults to "grey".
#' See persp3D.
#' @param shade The degree of shading of the surface facets. Defaults to
#' 0.25. See persp.
#' @param theta The azimuthal viewing direction. See persp.
#' @param xlab The x-axis label. Defaults to "T1".
#' @param ylab The y-axis label. Defaults to "T2".
#' @param zlab The z-axis label. Defaults to "Fhat".
#' @param ... Additional parameters to the persp3D function.
#' @seealso \code{\link{npSurv2}}, \code{\link[plot3D]{persp3D}}
#' @importFrom plot3D persp3D
#' @export
#' @examples
#' x <- genClayton2(100, 0, 1, 1, 2, 2)
#' x.npSurv2 <- npSurv2(x$Y1, x$Y2, x$Delta1, x$Delta2)
#' plotnpSurv2.3D(x.npSurv2)
#'
#' x2 <- genClayton2(100, 2, 1, 1, 2, 2)
#' x2.npSurv2 <- npSurv2(x2$Y1, x2$Y2, x2$Delta1, x2$Delta2)
#' plotnpSurv2.3D(x2.npSurv2)
plotnpSurv2.3D <- function(npSurv2.obj, col="grey", shade=0.25, theta=120,
                           xlab="T1", ylab="T2", zlab="Fhat", ...) {
    T1 <- npSurv2.obj$T1
    T2 <- npSurv2.obj$T2
    Fhat <- npSurv2.obj$Fhat
    plot3D::persp3D(c(0, T1), c(0, T2), Fhat, col=col, shade=shade,
                    theta=theta, xlab=xlab, ylab=ylab, zlab=zlab, ...)
}

#' Uses a heat map to visualize a nonparametric bivariate survival function
#'
#' Plots a heat map of an estimated nonparametric bivariate survival
#' function. This function is a wrapper for the image function with default
#' parameters chosen to make the data easier to visualize.
#'
#' @param npSurv2.obj Output of the npSurv2 function.
#' @param contour Should contour lines be added to the plot? Defaults to
#' TRUE.
#' @param col List of colors for the heat map. Defaults to
#' terrain.colors(100).
#' @param xlab The x-axis label. Defaults to "T1".
#' @param ylab The y-axis label. Defaults to "T2".
#' @param ... Additional parameters to the image function.
#' @seealso \code{\link{npSurv2}}, \code{\link[graphics]{image}}
#' @importFrom graphics image
#' @importFrom grDevices terrain.colors
#' @export
#' @examples
#' x <- genClayton2(1000, 0, 1, 1, 2, 2)
#' x.npSurv2 <- npSurv2(x$Y1, x$Y2, x$Delta1, x$Delta2)
#' plotnpSurv2.HM(x.npSurv2)
#'
#' x2 <- genClayton2(1000, 2, 1, 1, 2, 2)
#' x2.npSurv2 <- npSurv2(x2$Y1, x2$Y2, x2$Delta1, x2$Delta2)
#' plotnpSurv2.HM(x2.npSurv2)
plotnpSurv2.HM <- function(npSurv2.obj, contour=TRUE, col=terrain.colors(100),
                           xlab="T1", ylab="T2", ...) {
    T1 <- npSurv2.obj$T1
    T2 <- npSurv2.obj$T2
    Fhat <- npSurv2.obj$Fhat
    image(c(0, T1), c(0, T2), Fhat, col=col, xlab=xlab, ylab=ylab, ...)
    if (contour) {
        contour(c(0, T1), c(0, T2), Fhat, add=TRUE)
    }
}
