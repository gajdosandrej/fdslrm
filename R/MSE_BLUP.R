#' @title Mean Square Error of the Best Linear Unbiased Predictor
#'
#' @description
#' \code{MSE_BLUP(F, V, var, f, v)} calculates the Mean Square Error (MSE) of the Best Linear Unbiased Predictor (BLUP) of FDSLRM time series in time \eqn{t_{n+d}}.
#'
#' @param F design matrix for fixed effects.
#' @param V design matrix for random effects.
#' @param var variance parameters.
#' @param f values of functions forming the columns of design matrix \eqn{F} computed in time \eqn{t_{n+d}}.
#' @param v values of functions forming the columns of design matrix \eqn{V} computed in time \eqn{t_{n+d}}.
#'
#' @return MSE of the BLUP of \eqn{X_{t_{n+d}}}.
#'
#' @note Ver.: 13-Mar-2019 18:50:41.
#'
#' @example R/Examples/example_MSE_BLUP.R
#'
#' @export
#'
MSE_BLUP <- function(F, V, var, f, v){

        D <- diag(var[-1])
        H <- t(V) %*% V %*% D + var[1] * diag(length(var)-1)
        mse <- var[1] * t(c(f, D %*% v) %*% solve(rbind(cbind(t(F) %*% F, t(F) %*% V %*% D), cbind(t(V) %*% F, H))) %*% c(f, v))

        return(mse)

}
