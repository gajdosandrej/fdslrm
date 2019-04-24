#' @title Natural estimates
#'
#' @description
#' \code{NE(X, F, V)} calculates natural estimates (NE) of variance parameters in general FDSLRM.
#'
#' @param X time series realization.
#' @param F design matrix for fixed effects.
#' @param V design matrix for random effects.
#'
#' @return Natural estimates of variance parameters \eqn{\sigma^2, \sigma_1^2, ..., \sigma_l^2}.
#'
#' @note Ver.: 21-Jan-2019 20:40:07.
#'
#' @example R/Examples/example_NE.R
#'
#' @export
#'
NE <- function(X, F, V){

        num_of_fixed_ef <- ncol(F)

        start_of_random_ef <- num_of_fixed_ef + 1
        num_of_effects <- ncol(F) + ncol(V)
        n <- NROW(X)

        M <- data.frame(F, V)
        LinModel <- lm(X~0+., data=M)
        coef <- coefficients(LinModel)

        Ne <- coef[start_of_random_ef:num_of_effects]
        Ne <- Ne^2

        resid <- residuals(LinModel)
        var <- t(resid) %*% resid
        var <- var / (n - num_of_effects)

        Ne <- append(var, Ne)

        return(as.vector(Ne))
}

