#' @title Maximum likelihood estimates
#'
#' @description
#' \code{MLE_orth(X, F, V)} calculates maximum likelihood estimates (MLE) of variance parameters in orthogonal FDSLRM.
#'
#' @param X time series realization.
#' @param F design matrix for fixed effects.
#' @param V design matrix for random effects.
#'
#' @return MLE of variance parameters \eqn{\sigma_1^2, ..., \sigma_l^2, \sigma^2}.
#'
#' @note Ver.: 06-Feb-2019 09:31:12.
#'
#' @import CVXR
#'
#' @example R/Examples/example_MLE_orth.R
#'
#' @export
#'
MLE_orth <- function(X, F, V){

        n <- length(X)
        # k <- ncol(F)
        l <- ncol(V)

        # GF <- t(F) %*% F
        # InvGF <- solve(GF)
        # I <- diag(n)
        # PF <- F %*% InvGF %*% t(F)
        #MF <- I - PF
        MF <- makeM_F(F)
        GV <- t(V) %*% V
        p <- n - l

        e <- as.vector(MF %*% X)
        # ec <- t(X) %*% MF
        ee <- as.numeric(t(e) %*% e)
        eV <- t(e) %*% V
        Ve <- t(V) %*% e
        d <- Variable(l+1)
        logdetS <- p * log(d[1]) + sum(log(d[1] - GV %*% d[2:(l+1)]))
        obj <- Maximize(logdetS - ((d[1] * ee) - (eV %*% diag(d[2:(l+1)]) %*% Ve)))
        constr <- list(d[2:(l+1)] >= 0, d[1] >= max_entries(GV %*% d[2:(l+1)]))
        p_MLE <- Problem(obj, constr)

        sol <- solve(p_MLE)

        s <- 1 / sol$getValue(d)[1]
        sj <- vector()
        for(j in 2:(l+1)) {
                sj <- c(sj, sol$getValue(d)[j]/(sol$getValue(d)[1] * (sol$getValue(d)[1] - sol$getValue(d)[j] * GV[j-1,j-1])))
        }
        result <- c(s, sj)

        return(as.vector(result))

}
