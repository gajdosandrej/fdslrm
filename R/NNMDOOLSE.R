#' @title Non-negative modified double ordinary least squares estimates
#'
#' @description
#' \code{NNMDOOLSE(X, F, V)} calculates non-negative modified double ordinary least squares estimates (NNMDOOLSE) of variance parameters in FDSLRM.
#'
#' @param X time series realization.
#' @param F design matrix for fixed effects.
#' @param V design matrix for random effects.
#'
#' @return NNMDOOLSE of variance parameters \eqn{\sigma_1^2, ..., \sigma_l^2, \sigma^2}.
#'
#' @note Ver.: 06-Feb-2019 09:21:03.
#'
#' @import CVXR
#'
#' @example R/Examples/example_NNMDOOLSE.R
#'
#' @export
#'
NNMDOOLSE <- function(X, F, V){

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
        p_obj <- Minimize(sum_squares(SX - (s[1] %*% MF) - (MFV %*% diag(s[2:(l+1)]) %*% t(MFV))))
        constr <- list(s >= 0)
        prob <- Problem(p_obj, constr)

        sol <- solve(prob)

        return(c(sol$getValue(s)))

}
