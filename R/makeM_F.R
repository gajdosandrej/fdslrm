#' @title Projection matrix \eqn{M_F}
#'
#' @description
#' \code{makeM_F(F)} creates a projection matrix \eqn{M_F}. It is definied by the formula
#' \eqn{M_F = F((F'F)^(-1))F'} and projects into the space orthogonal to the column space
#' spanned by the columns of matrix \eqn{F}.
#'
#' @param F matrix.
#'
#' @return Projection matrix \eqn{M_F}.
#'
#' @note Ver.: 02-Jan-2019 13:49:17.
#'
#' @example R/Examples/example_makeM_F.R
#'
#' @export
#'
makeM_F <- function(F) {

        n <- nrow(F)
        c1 <- rep(1, n)
        I <- diag(c1)
        M_F <- I - F %*% solve((t(F) %*% F)) %*% t(F)

        return(M_F)
}
