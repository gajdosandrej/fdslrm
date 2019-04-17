#' @title Non-negative double ordinary least squares estimates
#'
#' @description
#' \code{NNDOOLSE(X, F, V)} calculates non-negative double ordinary least squares estimates (NNDOOLSE) of variance parameters in FDSLRM.
#'
#' @param X time series realization.
#' @param F design matrix for fixed effects.
#' @param V design matrix for random effects.
#'
#' @return NNDOOLSE of variance parameters \eqn{\sigma_1^2, ..., \sigma_l^2, \sigma^2}.
#'
#' @note Ver.: 06-Feb-2019 08:51:47.
#'
#' @import CVXR
#'
#' @example R/Examples/example_NNDOOLSE.R
#'
#' @export
#'
NNDOOLSE <- function(X, F, V){

        n <- length(X)
        k <- ncol(F)
        l <- ncol(V)

        # GF <- t(F) %*% F
        # InvGF <- solve(GF)
        I <- diag(n)
        # PF <- F %*% InvGF %*% t(F)
        #MF <- I - PF
        MF <- makeM_F(F)
        MFV <- MF %*% V

        SX <- MF %*% X %*% t(X) %*% MF
        s <- Variable(l+1)
        p_obj <- Minimize(sum_squares(SX - (s[1] %*% I) - (V %*% diag(s[2:(l+1)]) %*% t(V))))
        constr <- list(s >= 0)
        prob <- Problem(p_obj, constr)

        sol <- solve(prob)

        return(c(sol$getValue(s)))

}
