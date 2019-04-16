#' @title Initialize R packages and functions
#'
#' @description
#' \code{initialFDSLRM(supress_messages)} loads, installs all necessary R packages and functions
#' for FDSLRM fitting and diagnostics.
#'
#' @param supress_messages input parameter of type logical indicating
#' if messages during packages installation and loading should be hidden or not.
#' Default is \code{supress_messages = TRUE}.
#'
#' @note Ver.: 11-Jan-2019 20:33:50.
#'
#' @export
#'
initialFDSLRM <- function(supress_messages = TRUE) {

        # load utility functions
        # source("functions_fdslrm.R")
        # source("fit_diag_fdslrm.R")
        # source("drawTable.R")
        # source("drawDiagPlots.R")
        # source("SingerEtAl_LMM_diag.R")
        # source("SingerEtAl_justPlots256.R")
        # source("R/utilities.R")

        if(supress_messages == TRUE) {

                # install and/or load packages
                if("kableExtra" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("kableExtra")))
                }
                if("IRdisplay" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("IRdisplay")))
                }
                if("MASS" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("MASS")))
                }
                if("Matrix" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("Matrix")))
                }
                if("car" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("car")))
                }
                if("nlme" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("nlme")))
                }
                if("stats" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("stats")))
                }
                if("stats" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("forecast")))
                }
                if("stats" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("fpp2")))
                }
                if("matrixcalc" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("matrixcalc")))
                }
                if("sommer" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("sommer")))
                }
                if("gnm" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("gnm")))
                }
                if("pracma" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("pracma")))
                }
                if("CVXR" %in% rownames(installed.packages()) == FALSE) {
                        suppressWarnings(suppressMessages(install.packages("CVXR")))
                }

                suppressWarnings(suppressMessages(library(kableExtra)))
                suppressWarnings(suppressMessages(library(IRdisplay)))
                suppressWarnings(suppressMessages(library(MASS)))
                suppressWarnings(suppressMessages(library(Matrix)))
                suppressWarnings(suppressMessages(library(car)))
                suppressWarnings(suppressMessages(library(nlme)))
                suppressWarnings(suppressMessages(library(stats)))
                suppressWarnings(suppressMessages(library(forecast)))
                suppressWarnings(suppressMessages(library(fpp2)))
                suppressWarnings(suppressMessages(library(matrixcalc)))
                suppressWarnings(suppressMessages(library(sommer)))
                suppressWarnings(suppressMessages(library(gnm)))
                suppressWarnings(suppressMessages(library(pracma)))
                suppressWarnings(suppressMessages(library(CVXR)))


        } else {

                # install and/or load packages
                if("kableExtra" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("kableExtra")
                }
                if("IRdisplay" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("IRdisplay")
                }
                if("MASS" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("MASS")
                }
                if("Matrix" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("Matrix")
                }
                if("car" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("car")
                }
                if("nlme" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("nlme")
                }
                if("stats" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("stats")
                }
                if("stats" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("forecast")
                }
                if("stats" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("fpp2")
                }
                if("matrixcalc" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("matrixcalc")
                }
                if("sommer" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("sommer")
                }
                if("gnm" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("gnm")
                }
                if("pracma" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("pracma")
                }
                if("CVXR" %in% rownames(installed.packages()) == FALSE) {
                        install.packages("CVXR")
                }

                library(kableExtra)
                library(IRdisplay)
                library(MASS)
                library(Matrix)
                library(car)
                library(nlme)
                library(stats)
                library(forecast)
                library(fpp2)
                library(matrixcalc)
                library(sommer)
                library(gnm)
                library(pracma)
                library(CVXR)

        }


}
