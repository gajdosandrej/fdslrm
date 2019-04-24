#' @title Modified natural estimates
#'
#' @description
#' \code{NE(X, F, V)} calculates modified natural estimates (MNE) of variance parameters in general FDSLRM.
#' MNE are an unbiased version of NE.
#'
#' @param X time series realization.
#' @param F design matrix for fixed effects.
#' @param V design matrix for random effects.
#'
#' @return Modified natural estimates of variance parameters \eqn{\sigma^2, \sigma_1^2, ..., \sigma_l^2}.
#'
#' @note Ver.: 21-Jan-2019 20:48:29.
#'
#' @example R/Examples/example_MNE.R
#'
#' @export
#'
MNE <- function(X, F, V){

        M_F <- makeM_F(F)
        W_inv <- solve(t(V) %*% M_F %*% V)

        Ne <- NE(X, F, V)
        Mne <- Ne[1]

        for(i in 2:NROW(Ne)){
                Mne <- append(Mne, Ne[i] - W_inv[i-1,i-1] * Ne[1])
        }

        return(as.vector(Mne))
}


