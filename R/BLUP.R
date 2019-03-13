#' @title Best Linear Unbiased Predictor
#'
#' @description
#' \code{BLUP(X, F, V, var, f, v)} calculates the Best Linear Unbiased Predictor (BLUP) of FDSLRM time series in time \eqn{t_{n+d}}.
#'
#' @param X time series.
#' @param F design matrix for fixed effects.
#' @param V design matrix for random effects.
#' @param var variance parameters.
#' @param f values of functions forming the columns of design matrix \eqn{F} computed in time \eqn{t_{n+d}}.
#' @param v values of functions forming the columns of design matrix \eqn{V} computed in time \eqn{t_{n+d}}.
#'
#' @return BLUP of \eqn{X_{t_{n+d}}}.
#'
#' @note Ver.: 13-Mar-2019 18:00:08.
#'
#' @example R/Examples/example_BLUP.R
#'
#' @export
#'
BLUP <- function(X, F, V, var, f, v){

        D <- diag(var[-1])
        H <- t(V) %*% V %*% D + var[1] * diag(length(var)-1)
        blup <- t(c(f, D %*% v) %*% solve(rbind(cbind(t(F) %*% F, t(F) %*% V %*% D), cbind(t(V) %*% F, H))) %*% c(t(F) %*% X, t(V) %*% X))

        return(blup)

}
