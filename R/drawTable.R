#' @title Draw table of FDSLRM parameters estimates and frequencies
#'
#' @description
#' \code{drawTable(type, periodogram, frequencies, fixed_eff, random_eff, variances)}
#' creates a good looking table for estimates of FDSLRM parameters and frequencies.
#' This function is designed for Jupyter notebooks.
#'
#' @param type one of the following values if acceptable: \code{"periodogram", "fixed", "random", "variance"}.
#' @param periodogram output of function \code{spec.pgram}.
#' @param frequencies vector of (Fourier) frequencies in the form of fractions.
#' @param fixed_eff vector of fixed effects.
#' @param random_eff vector of random effects.
#' @param variances vector of variances. First element should be a variance of white noise
#' and the next elements should be variances of random effects.
#'
#' @return Table with corresponding values.
#'
#' @note Ver.: 12-Jan-2019 17:20:30.
#'
#' @example R/Examples/example_drawTable.R
#'
#' @export
#'
drawTable <- function(type, periodogram, frequencies = NULL, fixed_eff, random_eff, variances) {

        suppressWarnings(suppressMessages(require(kableExtra)))
        suppressWarnings(suppressMessages(require(IRdisplay)))

        if(!(type %in% c("periodogram", "fixed", "random", "variance"))) {
                stop("type should be specified as frequency/fixed/random/variance.")
        }

        if(type == "periodogram") {

                if(is.null(frequencies)) {

                        dtf_spec_freq <- t(data.frame(sort(periodogram$spec, decreasing = TRUE)[1:6], periodogram$freq[order(periodo$spec, decreasing = TRUE)][1:6]))
                        row.names(dtf_spec_freq) <- c("spectrum", "frequency (raw)")
                        kable(dtf_spec_freq, format = "html", caption = 'Frequencies by spectrum') %>%
                                kable_styling("striped") %>%
                                as.character() %>%
                                display_html()

                } else {

                        frequencies_names <- vector()
                        for(i in 1:length(frequencies)) {
                                frequencies_names <- c(frequencies_names, sprintf(gsub("pi", "\\\\pi", paste("$",frequencies[i],"$", sep = ""))))
                        }
                        dtf_spec_freq <- t(data.frame(sort(periodogram$spec, decreasing = TRUE)[1:length(frequencies)], periodogram$freq[order(periodo$spec, decreasing = TRUE)][1:length(frequencies)], frequencies_names))
                        row.names(dtf_spec_freq) <- c("spectrum", "frequency (raw)", "frequency")
                        kable(dtf_spec_freq, format = "html", caption = 'Frequencies by spectrum') %>%
                                kable_styling("striped") %>%
                                as.character() %>%
                                display_html()

                }

        } else if(type == "fixed") {

                dtf_fix_eff <- t(data.frame(fixed_eff))
                fixed_eff_names <- vector()
                for(i in 1:length(fixed_eff)) {
                        fixed_eff_names <- c(fixed_eff_names, sprintf(paste("$\\beta_{",i,"}$", sep = "")))
                }
                colnames(dtf_fix_eff) <- fixed_eff_names
                row.names(dtf_fix_eff) <- c("")
                kable(dtf_fix_eff, format = "html") %>%
                        kable_styling("striped") %>%
                        as.character() %>%
                        display_html()

        } else if(type == "random") {

                random_eff_vec <- as.vector(unlist(random_eff))
                dtf_rnd_eff <- t(data.frame(random_eff_vec))
                Y_names <- vector()
                for(i in 1:length(random_eff_vec)) {
                        Y_names <- c(Y_names, sprintf(paste("$Y_{",i,"}$", sep = "")))
                }
                colnames(dtf_rnd_eff) <- Y_names
                row.names(dtf_rnd_eff) <- c("")
                kable(dtf_rnd_eff, format = "html") %>%
                        kable_styling("striped") %>%
                        as.character() %>%
                        display_html()

        } else {
                variances_vec <- c(variances)
                # variances_vec <- c(variances_vec[length(variances_vec)], variances_vec[-length(variances_vec)])
                dtf_var <- t(data.frame(variances_vec))
                col_names <- c(sprintf("$\\sigma_{0}^2$"))
                for(i in 1:(length(variances_vec)-1)) {
                        col_names <- c(col_names, sprintf(paste("$\\sigma_{",i,"}^2$", sep = "")))
                }
                colnames(dtf_var) <- col_names
                row.names(dtf_var) <- NULL
                kable(dtf_var, format = "html") %>%
                        kable_styling("striped") %>%
                        as.character() %>%
                        display_html()

        }

}
