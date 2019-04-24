#' @title Modified double ordinary least squares estimates
#'
#' @description
#' \code{MDOOLSE_orth(X, F, V)} calculates modified double ordinary least squares estimates (MDOOLSE) of variance parameters in orthogonal FDSLRM.
#' MDOOLSE are an unbiased version of DOOLSE.
#'
#' @param X time series realization.
#' @param F design matrix for fixed effects.
#' @param V design matrix for random effects.
#'
#' @return MDOOLSE of variance parameters \eqn{\sigma^2, \sigma_1^2, ..., \sigma_l^2}.
#'
#' @note Ver.: 21-Jan-2019 21:56:04.
#'
#' @example R/Examples/example_MDOOLSE_orth.R
#'
#' @export
#'
MDOOLSE_orth <- function(X, F, V){

        V2 <- t(V) %*% V

        Ne <- NE(X, F, V)
        Mdoolse <- Ne[1]

        for(i in 2:length(Ne)){
                Mdoolse <- append(Mdoolse, Ne[i] - Ne[1] / V2[i-1,i-1])
        }

        return(as.vector(Mdoolse))

}

