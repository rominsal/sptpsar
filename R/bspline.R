#' Build design matrix for B-spline basis.
#'
#' @param x vector of data.
#' @param xl lower bound for knots.
#' @param xr upper bound for knots.
#' @param ndx number of inner knots.
#' @param bdeg degree of B-Spline basis (3 by default).
#' @return B design matrix of B-Spline basis.
#' @return knots vector of inner and outer knots.

bspline <- function(x, xl, xr, ndx, bdeg){
	dx <- (xr - xl)/ndx
	knots <- seq(xl - bdeg*dx, xr + bdeg*dx, by=dx)
	B <- splines::spline.des(knots, x, bdeg+1, 0*x)$design
	res <- list(B = B, knots = knots)
	res
}
