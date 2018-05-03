#' Build mixed-model matrices for P-Spline models.
#'
#' @param x vector of data.
#' @param xl lower bound for knots.
#' @param xr upper bound for knots.
#' @param ndx number of inner knots.
#' @param bdeg degree of B-Spline basis (3 by default).
#' @param pord degree of penalty differences (2 by default).
#' @param decom type of fixed effect matrix \code{X}.
#' @param \itemize{ \item 1 (default) \code{X=P*U} \item 2 \code{X=(1,x,...,x^{pord-1})} \item 3 FILL }
#' @return X design matrix of fixed-effects in mixed models.
#' @return Z design Matrix of random-effects in mixed models.
#' @return d vector of eigenvalues of SVD decomposition.
#' @return B design matrix of B-Spline basis.
#' @return m number of columns of B matrix.
#' @return D penalty matrix of differences.
#' @return U matrix including all vectors of null-space of B.
#' @return \code{knots} number of knots.

MM_basis <- function (x, xl, xr, ndx, bdeg, pord, decom = 1) {
	Bb <- bspline(x,xl,xr,ndx,bdeg)
	knots <- Bb$knots
	B <- Bb$B
	m <- ncol(B)
	n <- nrow(B)
	D <- diff(diag(m), differences=pord)
	P.svd <- svd(crossprod(D))
	U <- (P.svd$u)[,1:(m-pord)] # eigenvectors
	d <- (P.svd$d)[1:(m-pord)]  # eigenvalues
	Z <- B %*% U
	if(decom == 1) {
		X <- B %*% ((P.svd$u)[,-(1:(m-pord))])
	} else if (decom == 2){
		X <- NULL
		for(i in 0:(pord-1)){ X <- cbind(X,x^i) }
	} else if(decom == 3) {
		Xf <- NULL
		for(i in 0:(pord-1)){
	    	Xf <- cbind(Xf,knots[-c((1:pord),
	    	                  (length(knots)- pord + 1):length(knots))]^i)
		}
		X <- B %*% Xf
	}
	list(X = X, Z = Z, d = d, B = B, m = m, D = D, U = P.svd$u, knots = knots)
}
