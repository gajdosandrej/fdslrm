#' @title Draw diagnostic plots for FDSLRM
#'
#' @description
#' \code{drawDiagPlots(plots_names, fit_diag_fdslrm_output)} draws diagnostics plots of fitted FDSLRM.
#' Plot with original time series, fitted values and fitted mean values is included.
#'
#' @param plots_names vector of type character containing names of available diagnostic plots.
#' Default is \code{plots_names = "all"}.
#' See details section of documentation and example for more information.
#' @param fit_diag_fdslrm_output list - an output of function \code{fitDiagFDSLRM}.
#'
#' @details
#' You can find names of all available diagnostic plots in the list, which is an ouput of function \code{fitDiagFDSLRM}.
#' These names are under the item "diagnostic_plots_names".
#'
#' @return Diagnostic plots of particular fitted FDSLRM.
#'
#' @note Ver.: 12-Jan-2019 17:18:53.
#'
#' @example R/Examples/example_fitDiagFDSLRM.R
#'
#' @export
#'
drawDiagPlots <- function(plots_names = "all", fit_diag_fdslrm_output) {

        output <- fit_diag_fdslrm_output

        if(!is.null(output$lm_fit)) {

                        #lm
                        if(plots_names == "all") {

                                # windows()
                                par(mfrow = c(3, 4))

                                acf(output$raw_residuals, main = "ACF of raw residuals")
                                pacf(output$raw_residuals, main = "PACF of  raw residuals")
                                hist(output$raw_residuals, main = "Distribution of raw residuals", xlab = "Raw residuals")

                                acf(output$stud_residuals, main = "ACF of studentized residuals")
                                pacf(output$stud_residuals, main = "PACF of  studentized residuals")
                                hist(output$stud_residuals, main = "Distribution of studentized residuals", xlab = "Studentized residuals")

                                plot(output$stud_residuals, type = "o", xlab = "Time", ylab = "Studentized residuals", main = "")
                                grid()
                                plot(fitted(output$lm_fit), resid(output$lm_fit), xlab = "Fitted values", ylab = "Raw residuals", main = "")
                                grid()
                                plot(fitted(output$lm_fit), output$stud_residuals, xlab = "Fitted values", ylab = "Studentized residuals", main = "")
                                grid()

                                plot(output$raw_residuals, type = "o", xlab = "Time", ylab = "Raw residuals", main = "")
                                grid()
                                qqnorm(resid(output$lm_fit), ylab = "Raw residuals")
                                qqline(resid(output$lm_fit))
                                grid()
                                qqnorm(output$stud_residuals, ylab = "Studentized residuals")
                                qqline(output$stud_residuals)
                                grid()

                                par(mfrow = c(1, 1))
                                # car::spreadLevelPlot(output$lm_fit)
                                # grid()

                        } else {

                                if("FittedTimeSeries" %in% plots_names) {
                                        x_ts <- ts(output$time_series)
                                        fitted_values_ts <- ts(output$fitted_values)
                                        ts.plot(x_ts, fitted_values_ts, gpars = list(col = c("black", "blue"), lwd = 2,
                                                                                     xlab="Time", ylab="",main="Original time series vs fitted values"))
                                        legend("topleft", legend = c("time series", "fitted values"), col = c("black", "blue"),
                                               lwd = 2, cex = 0.5, bty = "n", bg="transparent")
                                        grid()
                                        par(mfrow = c(1, 1))
                                }

                                if("RawResid" %in% plots_names) {
                                        plot(output$raw_residuals, type = "o", xlab = "Time", ylab = "Raw residuals", main = "")
                                        grid()
                                        par(mfrow = c(1, 1))
                                }

                                if("StdResid" %in% plots_names) {
                                        plot(output$stud_residuals, type = "o", xlab = "Time", ylab = "Studentized residuals", main = "")
                                        grid()
                                        par(mfrow = c(1, 1))
                                }

                                if("FittedValuesVsRawResid" %in% plots_names) {
                                        plot(fitted(output$lm_fit), resid(output$lm_fit), xlab = "Fitted values", ylab = "Raw residuals", main = "")
                                        grid()
                                        par(mfrow = c(1, 1))
                                }

                                if("FittedValuesVsStdResid" %in% plots_names) {
                                        plot(fitted(output$lm_fit), output$stud_residuals, xlab = "Fitted values", ylab = "Studentized residuals", main = "")
                                        grid()
                                        par(mfrow = c(1, 1))
                                }

                                if("SpreadLevelPlot" %in% plots_names) {
                                        car::spreadLevelPlot(output$lm_fit)
                                        grid()
                                }

                                if("QQPlotRawResid" %in% plots_names) {
                                        qqnorm(resid(output$lm_fit), ylab = "Raw residuals")
                                        qqline(resid(output$lm_fit))
                                        grid()
                                        par(mfrow = c(1, 1))
                                }

                                if("QQPlotStdResid" %in% plots_names) {
                                        qqnorm(output$stud_residuals, ylab = "Studentized residuals")
                                        qqline(output$stud_residuals)
                                        grid()
                                        par(mfrow = c(1, 1))
                                }

                                if("ACFRawResid" %in% plots_names) {
                                        acf(output$raw_residuals, main = "ACF of raw residuals")
                                        par(mfrow = c(1, 1))
                                }

                                if("PACFRawResid" %in% plots_names) {
                                        pacf(output$raw_residuals, main = "PACF of  raw residuals")
                                        par(mfrow = c(1, 1))
                                }

                                if("HistogramRawResid" %in% plots_names) {
                                        hist(output$raw_residuals, main = "Distribution of raw residuals", xlab = "Raw residuals")
                                        par(mfrow = c(1, 1))
                                }

                                if("ACFStdResid" %in% plots_names) {
                                        acf(output$stud_residuals, main = "ACF of studentized residuals")
                                        par(mfrow = c(1, 1))
                                }

                                if("PACFStdResid" %in% plots_names) {
                                        pacf(output$stud_residuals, main = "PACF of  studentized residuals")
                                        par(mfrow = c(1, 1))
                                }

                                if("HistogramStdResid" %in% plots_names) {
                                        hist(output$stud_residuals, main = "Distribution of studentized residuals", xlab = "Studentized residuals")
                                        par(mfrow = c(1, 1))
                                }

                                if("CumulatPeriodogStudResid" %in% plots_names) {
                                        cpgram(output$stud_residuals, main = "Cumulative periodogram of studentized residuals", ci.col = "red")
                                        par(mfrow = c(1, 1))
                                }

                        }

                } else if (!is.null(output$lme_fit)) {
                        #lme
                        if(plots_names == "all") {

                                # windows()
                                par(mfrow = c(3, 3))

                                SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(2), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V))
                                plot(output$diag_resid$std.marginal.residuals[,2], type = "o", ylab = "Standardized marginal residuals")
                                abline(0,0)
                                grid()
                                acf(output$cond_residuals, main = "ACF of conditional residuals")

                                SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(5), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V))
                                plot(output$diag_resid$std.conditional.residuals[,2], type = "o", ylab = "Standardized conditional residuals")
                                abline(0,0)
                                grid()
                                pacf(output$cond_residuals, main = "PACF of conditional residuals")

                                hist(output$norm_cond_resid, main = "Distribution of norm. cond. residuals", xlab = "Norm. cond. residuals")
                                hist(output$diag_resid$least.confounded.residuals, main = "Distribution of stand. least conf. residuals", xlab = "Stand. least conf. residuals")
                                SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(6), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V))

                                par(mfrow = c(1, 1))

                        } else {

                                if("StdMarginalResidVsFittedValues" %in% plots_names) {
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(2), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                }

                                if("StdConditionalResidVsFittedValues" %in% plots_names) {
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(5), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                }

                                if("NormalQQPlotStdLeastConfCondResid" %in% plots_names) {
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(6), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                }

                                if("FittedTimeSeries" %in% plots_names) {
                                        x_ts <- ts(output$time_series)
                                        fitted_values_ts <- ts(output$cond_fitted_values)
                                        trend_ts <- ts(as.vector(output$matrix_F%*%output$fixed_effects))
                                        ts.plot(x_ts, fitted_values_ts, trend_ts, gpars = list(col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = 2,
                                                                                               xlab="Time", ylab="",main="Original time series vs fitted values"))
                                        legend("topleft", legend = c("time series", "fitted values", "estimated trend"), col = c("black", "blue", "red"),
                                               lty = c(1, 1, 2), lwd = 2, cex = 0.5, bty = "n", bg="transparent")
                                        grid()
                                        par(mfrow = c(1, 1))
                                }

                                if("StdMarginalResid" %in% plots_names) {
                                        plot(output$diag_resid$std.marginal.residuals[,2], type = "o", ylab = "Standardized marginal residuals")
                                        abline(0,0)
                                        grid()
                                        par(mfrow = c(1, 1))
                                }

                                if("StdConditionalResid" %in% plots_names) {
                                        plot(output$diag_resid$std.conditional.residuals[,2], type = "o", ylab = "Standardized conditional residuals")
                                        abline(0,0)
                                        grid()
                                        par(mfrow = c(1, 1))
                                }

                                if("ACFCondtResid" %in% plots_names) {
                                        acf(output$cond_residuals, main = "ACF of conditional residuals")
                                        par(mfrow = c(1, 1))
                                }

                                if("PACFCondResid" %in% plots_names) {
                                        pacf(output$cond_residuals, main = "PACF of conditional residuals")
                                        par(mfrow = c(1, 1))
                                }

                                if("HistNormCondResid" %in% plots_names) {
                                        hist(output$norm_cond_resid, main = "Distribution of norm. cond. residuals", xlab = "Norm. cond. residuals")
                                        par(mfrow = c(1, 1))
                                }

                                if("HistLeastConfCondResid" %in% plots_names) {
                                        hist(output$diag_resid$least.confounded.residuals, main = "Distribution of stand. least conf. residuals", xlab = "Stand. least conf. residuals")
                                        par(mfrow = c(1, 1))
                                }

                                if("CumulatPeriodogCondResid" %in% plots_names) {
                                        cpgram(output$cond_residuals, main = "Cumulative periodogram of conditional residuals", ci.col = "red")
                                        par(mfrow = c(1, 1))
                                }

                        }

        } else if(!is.null(output$mixed_fit)){
                # mixed

                if(plots_names == "all") {

                        # windows()
                        par(mfrow = c(3, 3))

                        x_ts <- ts(output$time_series)
                        fitted_values_ts <- ts(output$cond_fitted_values)
                        trend_ts <- ts(as.vector(output$matrix_F%*%output$fixed_effects))
                        ts.plot(x_ts, fitted_values_ts, trend_ts, gpars = list(col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = 2,
                                                                               xlab="Time", ylab="",main="Original time series vs fitted values"))
                        legend("topleft", legend = c("time series", "fitted values", "estimated trend"), col = c("black", "blue", "red"),
                               lty = c(1, 1, 2), lwd = 2, cex = 0.5, bty = "n", bg="transparent")
                        grid()

                        plot(output$marg_residuals, type = "o", ylab = "Marginal residuals")
                        abline(0,0)
                        grid()

                        plot(output$marg_fitted_values, output$marg_residuals, ylab = "Marginal residuals", xlab = "Marginal fitted values")
                        abline(0,0)
                        grid()

                        plot(output$cond_residuals, type = "o", ylab = "Conditional residuals")
                        abline(0,0)
                        grid()

                        plot(output$cond_fitted_values, output$cond_residuals, ylab = "Conditional residuals", xlab = "Conditional fitted values")
                        abline(0,0)
                        grid()

                        acf(output$cond_residuals, main = "ACF of conditional residuals")
                        pacf(output$cond_residuals, main = "PACF of conditional residuals")
                        hist(output$marg_residuals, main = "Distribution of marginal residuals", xlab = "Marginal residuals")
                        hist(output$cond_residuals, main = "Distribution of conditional residuals", xlab = "Conditional residuals")

                        par(mfrow = c(1, 1))

                } else {

                        if("MarginalResidVsFittedValues" %in% plots_names) {
                                plot(output$marg_fitted_values, output$marg_residuals, ylab = "Marginal residuals", xlab = "Marginal fitted values")
                                abline(0,0)
                                grid()
                                par(mfrow = c(1, 1))
                        }

                        if("ConditionalResidVsFittedValues" %in% plots_names) {
                                plot(output$cond_fitted_values, output$cond_residuals, ylab = "Conditional residuals", xlab = "Conditional fitted values")
                                abline(0,0)
                                grid()
                                par(mfrow = c(1, 1))
                        }

                        if("FittedTimeSeries" %in% plots_names) {
                                x_ts <- ts(output$time_series)
                                fitted_values_ts <- ts(output$cond_fitted_values)
                                trend_ts <- ts(as.vector(output$matrix_F%*%output$fixed_effects))
                                ts.plot(x_ts, fitted_values_ts, trend_ts, gpars = list(col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = 2,
                                                                                       xlab="Time", ylab="",main="Original time series vs fitted values"))
                                legend("topleft", legend = c("time series", "fitted values", "estimated trend"), col = c("black", "blue", "red"),
                                       lty = c(1, 1, 2), lwd = 2, cex = 0.5, bty = "n", bg="transparent")
                                grid()
                                par(mfrow = c(1, 1))
                        }

                        if("MarginalResid" %in% plots_names) {
                                plot(output$marg_residuals, type = "o", ylab = "Marginal residuals")
                                abline(0,0)
                                grid()
                                par(mfrow = c(1, 1))
                        }

                        if("ConditionalResid" %in% plots_names) {
                                plot(output$cond_residuals, type = "o", ylab = "Conditional residuals")
                                abline(0,0)
                                grid()
                                par(mfrow = c(1, 1))
                        }

                        if("ACFCondtResid" %in% plots_names) {
                                acf(output$cond_residuals, main = "ACF of conditional residuals")
                                par(mfrow = c(1, 1))
                        }

                        if("PACFCondResid" %in% plots_names) {
                                pacf(output$cond_residuals, main = "PACF of conditional residuals")
                                par(mfrow = c(1, 1))
                        }

                        if("HistCondResid" %in% plots_names) {
                                hist(output$cond_residuals, main = "Distribution of conditional residuals", xlab = "Conditional residuals")
                                par(mfrow = c(1, 1))
                        }

                        if("HistMargResid" %in% plots_names) {
                                hist(output$marg_residuals, main = "Distribution of marginal residuals", xlab = "Marginal residuals")
                                par(mfrow = c(1, 1))
                        }

                        if("CumulatPeriodogCondResid" %in% plots_names) {
                                cpgram(output$cond_residuals, main = "Cumulative periodogram of conditional residuals", ci.col = "red")
                                par(mfrow = c(1, 1))
                        }
                }

        } else {
                # mmer
                if(plots_names == "all") {

                        # windows()
                        par(mfrow = c(3, 3))

                        x_ts <- ts(output$time_series)
                        fitted_values_ts <- ts(output$cond_fitted_values)
                        trend_ts <- ts(as.vector(output$matrix_F%*%output$fixed_effects))
                        ts.plot(x_ts, fitted_values_ts, trend_ts, gpars = list(col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = 2,
                                                                               xlab="Time", ylab="",main="Original time series vs fitted values"))
                        legend("topleft", legend = c("time series", "fitted values", "estimated trend"), col = c("black", "blue", "red"),
                               lty = c(1, 1, 2), lwd = 2, cex = 0.5, bty = "n", bg="transparent")
                        grid()

                        plot(output$marg_residuals, type = "o", ylab = "Marginal residuals")
                        abline(0,0)
                        grid()

                        plot(output$marg_fitted_values, output$marg_residuals, ylab = "Marginal residuals", xlab = "Marginal fitted values")
                        abline(0,0)
                        grid()

                        plot(output$cond_residuals, type = "o", ylab = "Conditional residuals")
                        abline(0,0)
                        grid()

                        plot(output$cond_fitted_values, output$cond_residuals, ylab = "Conditional residuals", xlab = "Conditional fitted values")
                        abline(0,0)
                        grid()

                        acf(output$cond_residuals, main = "ACF of conditional residuals")
                        pacf(output$cond_residuals, main = "PACF of conditional residuals")
                        hist(output$marg_residuals, main = "Distribution of marginal residuals", xlab = "Marginal residuals")
                        hist(output$cond_residuals, main = "Distribution of conditional residuals", xlab = "Conditional residuals")

                        par(mfrow = c(1, 1))

                } else {

                        if("MarginalResidVsFittedValues" %in% plots_names) {
                                plot(output$marg_fitted_values, output$marg_residuals, ylab = "Marginal residuals", xlab = "Marginal fitted values")
                                abline(0,0)
                                grid()
                                par(mfrow = c(1, 1))
                        }

                        if("ConditionalResidVsFittedValues" %in% plots_names) {
                                plot(output$cond_fitted_values, output$cond_residuals, ylab = "Conditional residuals", xlab = "Conditional fitted values")
                                abline(0,0)
                                grid()
                                par(mfrow = c(1, 1))
                        }

                        if("FittedTimeSeries" %in% plots_names) {
                                x_ts <- ts(output$time_series)
                                fitted_values_ts <- ts(output$cond_fitted_values)
                                trend_ts <- ts(as.vector(output$matrix_F%*%output$fixed_effects))
                                ts.plot(x_ts, fitted_values_ts, trend_ts, gpars = list(col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = 2,
                                                                                       xlab="Time", ylab="",main="Original time series vs fitted values"))
                                legend("topleft", legend = c("time series", "fitted values", "estimated trend"), col = c("black", "blue", "red"),
                                       lty = c(1, 1, 2), lwd = 2, cex = 0.5, bty = "n", bg="transparent")
                                grid()
                                par(mfrow = c(1, 1))
                        }

                        if("MarginalResid" %in% plots_names) {
                                plot(output$marg_residuals, type = "o", ylab = "Marginal residuals")
                                abline(0,0)
                                grid()
                                par(mfrow = c(1, 1))
                        }

                        if("ConditionalResid" %in% plots_names) {
                                plot(output$cond_residuals, type = "o", ylab = "Conditional residuals")
                                abline(0,0)
                                grid()
                                par(mfrow = c(1, 1))
                        }

                        if("ACFCondtResid" %in% plots_names) {
                                acf(output$cond_residuals, main = "ACF of conditional residuals")
                                par(mfrow = c(1, 1))
                        }

                        if("PACFCondResid" %in% plots_names) {
                                pacf(output$cond_residuals, main = "PACF of conditional residuals")
                                par(mfrow = c(1, 1))
                        }

                        if("HistCondResid" %in% plots_names) {
                                hist(output$cond_residuals, main = "Distribution of conditional residuals", xlab = "Conditional residuals")
                                par(mfrow = c(1, 1))
                        }

                        if("HistMargResid" %in% plots_names) {
                                hist(output$marg_residuals, main = "Distribution of marginal residuals", xlab = "Marginal residuals")
                                par(mfrow = c(1, 1))
                        }

                        if("CumulatPeriodogCondResid" %in% plots_names) {
                                cpgram(output$cond_residuals, main = "Cumulative periodogram of conditional residuals", ci.col = "red")
                                par(mfrow = c(1, 1))
                        }
                }
        }
}
