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
        logdetS <- p * log(d[l+1]) + sum(log(d[l+1] - GV %*% d[1:l]))
        obj <- Maximize(logdetS - ((d[l+1] * ee) - (eV %*% diag(d[1:l]) %*% Ve)))
        constr <- list(d[1:l] >= 0, d[l+1] >= max_entries(GV %*% d[1:l]))
        p_MLE <- Problem(obj, constr)

        sol <- solve(p_MLE)

        s <- 1 / sol$getValue(d)[l+1]
        sj <- vector()
        for(j in 1:l) {
                sj <- c(sj, sol$getValue(d)[j]/(sol$getValue(d)[l+1] * (sol$getValue(d)[l+1] - sol$getValue(d)[j] * GV[j,j])))
        }
        result <- c(sj, s)

        return(result)

}
