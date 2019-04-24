#' @title Double ordinary least squares estimates
#'
#' @description
#' \code{DOOLSE_orth(X, F, V)} calculates double ordinary least squares estimates (DOOLSE) of variance parameters in orthogonal FDSLRM.
#'
#' @param X time series realization.
#' @param F design matrix for fixed effects.
#' @param V design matrix for random effects.
#'
#' @return DOOLSE of variance parameters \eqn{\sigma^2, \sigma_1^2, ..., \sigma_l^2}.
#'
#' @note Ver.: 21-Jan-2019 21:50:44.
#'
#' @example R/Examples/example_DOOLSE_orth.R
#'
#' @export
#'
DOOLSE_orth <- function(X, F, V){

        V2 <- t(V) %*% V

        num_of_random <- ncol(V)
        num_of_effects <- ncol(F) + ncol(V)
        n <- NROW(X)
        const <- (n - num_of_effects) / (n - num_of_random)

        Ne <- NE(X, F, V)
        Doolse <- Ne[1] * const

        for(i in 2:length(Ne)){
                Doolse <- append(Doolse, Ne[i] - const * Ne[1] / V2[i-1,i-1])
        }

        return(as.vector(Doolse))
}
