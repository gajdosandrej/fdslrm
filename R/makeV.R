#' @title Design matrix for random effects in FDSLRM
#'
#' @description
#' \code{makeV(times, freqs)} creates a design matrix \eqn{V} for random effects (including an intercept)
#' in FDSLRM (for each frequency \eqn{f} produces two columns of \eqn{V}:
#' \eqn{sin(2*\pi*f*times}), \eqn{cos(2*\pi*f*times}).
#' \code{makeV} is a special case of \code{makeF} except for the first column of \eqn{1} for the intercept.
#'
#' @param times vector of observation times.
#' @param freqs vector of (Fourier) frequencies.
#'
#' @return Design matrix \eqn{V} for random effects in FDSLRM.
#'
#' @note Ver.: 02-Jan-2019 13:41:24.
#'
#' @example R/Examples/example_makeV.R
#'
#' @export
#'
makeV <- function(times, freqs) {

        V <- makeF(times, freqs)
        V <- V[, -1]

        return(V)
}
