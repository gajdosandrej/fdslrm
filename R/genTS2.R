#' @title Generate FDSLRM time series realization
#'
#' @description
#' \code{genTS2(model_formula, beta, nu)} generates a realization of a time series (FDSLRM) with mean value \code{mean},
#' variance parameters of white noise and random effects \code{var} and design matrix for random effects \code{V}.
#' Normal distribution of white noise and random effects is assumed.
#'
#' @param model_formula object of class \code{formula} which contains a symbolic model formula
#' including functions \eqn{f_i(t)} and \eqn{v_j(t)} neccessary to create design matrices \eqn{F}, \eqn{V}
#' for fixed and random effects in FDSLRM.
#' @param beta vector of mean value parameters.
#' @param nu vector of variance parameters,
#' the first one is for the white noise and the next are for the random effects.
#' @param n length of realization.
#'
#' @details explain the form of FDSLRM time series realization as linear mixed model?
#'
#' @return Realization of particular FDSLRM time series.
#'
#' @note Ver.: 06-Jan-2019 19:15:30.
#'
#' @example R/Examples/example_genTS2.R
#'
#' @export
#'
genTS2 <- function(model_formula, beta, nu, n){

        tsmf <- model.frame(model_formula, data.frame(t = 1:n, x = rep(1,n)))
        mat_FV <- model.matrix(model_formula, tsmf)

        k <- length(beta)

        if(length(nu) > 1) {

                # l <- length(nu)-1
                mat_F <- mat_FV[,1:k]
                meanv <- mat_F %*% beta
                mat_V <- mat_FV[,(k+1):ncol(mat_FV)]
                X <- genTS(meanv, mat_V, nu)

                return(X)

        } else {

                mat_F <- mat_FV
                meanv <- mat_F %*% beta
                X <- genTS(meanv, var = nu)

                return(X)

        }

}
