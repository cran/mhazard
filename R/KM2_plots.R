#' Uses a 3D perspective plot to visualize a bivariate survival function
#'
#' Plots a 3D perspective plot of an estimated bivariate survival
#' function. This function is a wrapper for the persp3D function
#' from the plot3D package with default parameters chosen to make the
#' data easier to visualize.
#'
#' @param km2.obj Output of the KM2 function.
#' @param col Color palette to be used for the plot. Defaults to "grey".
#' See persp3D.
#' @param shade The degree of shading of the surface facets. Defaults to
#' 0.25. See persp.
#' @param theta The azimuthal viewing direction. See persp.
#' @param xlab The x-axis label. Defaults to "T1".
#' @param ylab The y-axis label. Defaults to "T2".
#' @param zlab The z-axis label. Defaults to "Fhat".
#' @param ... Additional parameters to the persp3D function.
#' @seealso \code{\link{mhazard-deprecated}}, \code{\link[plot3D]{persp3D}}
#' @examples
#' \donttest{x <- genClayton2(1000, 0, 1, 1, 2, 2)}
#' \donttest{x.km2 <- KM2(x$Y1, x$Y2, x$Delta1, x$Delta2)}
#' \donttest{plotKM2.3D(x.km2)}
#'
#' \donttest{x2 <- genClayton2(1000, 2, 1, 1, 2, 2)}
#' \donttest{x2.km2 <- KM2(x2$Y1, x2$Y2, x2$Delta1, x2$Delta2)}
#' \donttest{plotKM2.3D(x2.km2)}
#' @name plotKM2.3D-deprecated
#' @usage plotKM2.3D(km2.obj, col="grey", shade=0.25, theta=120,
#'                   xlab="T1", ylab="T2", zlab="Fhat", ...)
#' @keywords internal
NULL

#' @rdname mhazard-deprecated
#' @section \code{plotKM2.3D}:
#' For \code{plotKM2.3D}, use \code{\link{plotnpSurv2.3D}}.
#'
#' @export
plotKM2.3D <- function(km2.obj, col="grey", shade=0.25, theta=120,
                        xlab="T1", ylab="T2", zlab="Fhat", ...) {
    .Deprecated("plotnpSurv2.3D", package="mhazard")
    plotnpSurv2.3D(km2.obj, col, shade, theta, xlab, ylab, zlab, ...)
}

#' Uses a heat map to visualize a bivariate survival function
#'
#' Plots a heat map of an estimated bivariate survival function. This
#' function is a wrapper for the image function with default parameters
#' chosen to make the data easier to visualize.
#'
#' @param km2.obj Output of the KM2 function.
#' @param contour Should contour lines be added to the plot? Defaults to
#' TRUE.
#' @param col List of colors for the heat map. Defaults to
#' terrain.colors(100).
#' @param xlab The x-axis label. Defaults to "T1".
#' @param ylab The y-axis label. Defaults to "T2".
#' @param ... Additional parameters to the image function.
#' @seealso \code{\link{mhazard-deprecated}}, \code{\link[graphics]{image}}
#' @examples
#' \donttest{x <- genClayton2(1000, 0, 1, 1, 2, 2)}
#' \donttest{x.km2 <- KM2(x$Y1, x$Y2, x$Delta1, x$Delta2)}
#' \donttest{plotKM2.HM(x.km2)}
#'
#' \donttest{x2 <- genClayton2(1000, 2, 1, 1, 2, 2)}
#' \donttest{x2.km2 <- KM2(x2$Y1, x2$Y2, x2$Delta1, x2$Delta2)}
#' \donttest{plotKM2.HM(x2.km2)}
#' @name plotKM2.HM-deprecated
#' @usage plotKM2.HM(km2.obj, contour=TRUE, col=terrain.colors(100),
#'                   xlab="T1", ylab="T2", ...)
#' @keywords internal
NULL

#' @rdname mhazard-deprecated
#' @section \code{plotKM2.HM}:
#' For \code{plotKM2.HM}, use \code{\link{plotnpSurv2.HM}}.
#'
#' @export
plotKM2.HM <- function(km2.obj, contour=TRUE, col=terrain.colors(100),
                        xlab="T1", ylab="T2", ...) {
    .Deprecated("plotnpSurv2.HM", package="mhazard")
    plotnpSurv2.HM(km2.obj, contour, col, xlab, ylab, ...)
}
