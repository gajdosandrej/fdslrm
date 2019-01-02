#' @title Design matrix for fixed effects in FDSLRM
#'
#' @description
#' \code{makeF(times, freqs)} creates a design matrix \eqn{F} for fixed effects (including an intercept)
#' in FDSLRM (for each frequency \eqn{f} produces two columns of \eqn{F}:
#' \eqn{sin(2*\pi*f*times}), \eqn{cos(2*\pi*f*times}).
#'
#' @param times vector of observation times.
#' @param freqs vector of (Fourier) frequencies.
#'
#' @return Design matrix \eqn{F} for fixed effects in FDSLRM.
#'
#' @note Ver.: 02-Jan-2019 13:28:18.
#'
#' @example R/Examples/example_makeF.R
#'
#' @export
#'
makeF <- function(times, freqs) {

        n <- length(times)
        f <- length(freqs)
        c1 <- matrix(1, n)
        F <- matrix(c1, n, 1)

        for (i in 1:f) {
                F <- cbind(F, cos(2 * pi * freqs[i] * times))
                F <- cbind(F, sin(2 * pi * freqs[i] * times))
        }

        return(F)
}
