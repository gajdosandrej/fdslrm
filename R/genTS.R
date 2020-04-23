#' @title Generate FDSLRM time series realization
#'
#' @description
#' \code{genTS(mean, V, var)} generates a realization of a time series (FDSLRM) with mean value \code{mean},
#' variance parameters of white noise and random effects \code{var} and design matrix for random effects \code{V}.
#' Normal distribution of white noise and random effects is assumed.
#'
#' @param mean vector of mean values, \eqn{mean = F\beta}.
#' @param V design matrix for random effects.
#' @param var vector of variance parameters,
#' the first one is for the white noise and the next are for the random effects.
#'
#' @details explain the form of FDSLRM time series realization as linear mixed model?
#'
#' @return Realization of particular FDSLRM time series.
#'
#' @note Ver.: 03-Jan-2019 18:56:51.
#'
#' @example R/Examples/example_genTS.R
#'
#' @export
#'
genTS <- function(mean, V, var) {

        if(missing(mean)) {

                stop("Enter the mean values vector.")

        }

        if(missing(V) && length(var) > 1) {

                stop("Design matrix for random effects is missing.")

        } else if(missing(V) && length(var) == 1) {

                sd <- sqrt(var)
                n <- NROW(mean)
                wn <- rnorm(n, 0, sd[1])
                sim <- c()

                for(i in 1:n){
                        p <- mean[i] + wn[i]
                        sim <- append(sim, p)
                }

                return(sim)

        } else {

                sd <- sqrt(var)
                n <- NROW(mean)
                l <- ncol(V)
                wn <- rnorm(n,0,sd[1])
                randEffects <- c()

                for (i in 1:l){
                        randEffects <- append(randEffects, rnorm(1, 0, sd[i+1]))
                }

                sim <- c()

                for(i in 1:n){
                        p <- mean[i] + V[i,] %*% randEffects + wn[i]
                        sim <- append(sim, p)
                }

                return(sim)
        }


}
