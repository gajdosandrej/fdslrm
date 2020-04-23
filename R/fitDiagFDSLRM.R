#' @title Fitting and diagnostic of FDSLRM
#'
#' @description
#' \code{fitDiagFSDLRM(x, times, freq_mean, poly_trend_degree = 0, include_fixed_eff, freq_random, include_random_eff, var_estim_method = "REML", fit_function = "lme", season_period = 0, display_plots = FALSE)}
#' fits particular FDSLRM (in the form of linear mixed model - LMM or classical linear regression model - CLRM)
#' to time series data, returns estimates of all model's parameters
#' and it also returns the model diagnostic (numerical and graphical).
#'
#' @param x time series observation.
#' @param times vector of observation times.
#' @param freq_mean vector of (Fourier) frequencies to create the design matrix for fixed effects.
#' @param poly_trend_degree degree of polynomial in linear regression representing the mean value.
#' Default is \code{poly_trend_degree = 0}. Available just for \code{fit_function = lme}.
#' @param include_fixed_eff vector of zeros / ones to indicate which cosine and sine component
#' corresponding to particular frequency in \code{freq_mean} should be excluded / included.
#' @param freq_random vector of (Fourier) frequencies to create the design matrix for random effects.
#' @param include_random_eff vector of zeros / ones to indicate which cosine and sine component
#' corresponding to particular frequency in \code{freq_random} should be excluded / included.
#' @param fit_function type of function to fit the model. One of the following types is acceptable: \code{"lme", "mixed", "mmer"}.
#' Default is \code{fit_function = "lme"}.
#' @param var_estim_method type of method for estimating variance parameters. Default is \code{"REML"}.
#' For \code{fit_function = "lme"} the \code{"REML"} and \code{"ML"} estimates are available.
#' For \code{fit_function = "mixed"} the following estimates are available: \code{"ML", "REML", "MINQEI", "MINQEUI"}.
#' @param season_period input parameter of type numeric expressing the period of seasonality.
#' Default is \code{season_period = 0}.
#' @param display_plots input parameter of type logical indicating if all diagnostic plots should be drawn.
#' Default value is \code{display_plots = FALSE}.
#'
#' @details
#' Fitting functions \code{lme} (\code{lm} in case of CLRM), \code{mixed} (Witkovsky, 2002; Gajdos 2018),
#' \code{mmer} (Covarrubias-Pazaran, 2016)
#' consequently offer a little bit different diagnostic tools and outputs.
#' For example when \code{lme} is used the residual diagnostic is created with the help
#' of open source R code from Singer et al. (2017).
#'
#' If there are no frequencies for random effects as input, just the frequencies for fixed effects
#' then the CLRM is fitted automatically by the \code{lm} function.
#'
#' @references
#' Covarrubias-Pazaran G. 2016: Genome assisted prediction of quantitative traits using the R package sommer.
#' PLoS ONE 11(6):1-15.
#'
#' Gajdos, A. (2018): MMEinR: R version of Witkovsky's Matlab implementation
#' of iterative solving of Henderson's mixed model equations.
#' https://github.com/gajdosandrej/MMEinR.
#'
#' Pinheiro J., Bates D., DebRoy S., Sarkar D. and R Core Team (2019): _nlme:
#' Linear and Nonlinear Mixed Effects Models_. R package version 3.1-131,
#' <URL: https://CRAN.R-project.org/package=nlme>.
#'
#' Singer, J. M., Rocha, F. M. M., and Nobre, J. S. (2017): Graphical Tools for Detecting Departures
#' from Linear Mixed Model Assumptions and Some Remedial Measures.
#' \emph{International Statistical Review}, 85: 290-324.
#'
#' Witkovsky, V.: MATLAB Algorithm for solving Henderson's Mixed Model Equations. Technical Report No. 3/2000,
#' Institute of Measurement Science, Slovak Academy of Sciences, Bratislava, 2000.
#'
#' @return This function returns an output - \code{result}, which is a list with the following elements (for LMM):
#' \itemize{
#' \item \code{result$fixed_effects} estimates of fixed effects (regression parameter \eqn{\beta}),
#' \item \code{result$random_effects} BLUPs of random effects \eqn{Y},
#' \item \code{result$error_variance} estimate of variance of white noise component,
#' \item \code{result$rand_eff_variance} estimates of variances of random effects,
#' \item \code{result$marg_fitted_values} marginal fitted values (i.e. trend i.e. mean value estimate),
#' \item \code{result$cond_fitted_values} conditional fitted values (estimate of trend plus BLUP of random effects),
#' \item \code{result$marg_residuals} marginal residuals,
#' \item \code{result$cond_residuals} conditional residuals,
#' \item \code{result$norm_cond_resid} normalized conditional residuals,
#' \item \code{result$pearson_resid} Pearson's residuals,
#' \item \code{result$Box_test} Box test of conditional residuals independence (not adjusted for time series (non)-seasonality),
#' \item \code{result$BoxLjung_test} Box-Ljung test of conditional residuals independence (not adjusted for time series (non)-seasonality),
#' \item \code{result$ShapiroWilk_test_norm_cond_resid} Shapiro-Wilk normality test of normalized residuals,
#' \item \code{result$ShapiroWilk_test_stand_least_conf_resid}Shapiro-Wilk normality test of least confounded conditional residuals,
#' \item \code{result$lme_fit} object inheriting from class \code{lme}, an output of function \code{lme} containing details about model fitting,
#' \item \code{result$fit_summary} summary info about model fitting,
#' \item \code{result$diag_resid} an ouput (list) of diagnostic function created by Singer et al. (2017),
#' \item \code{result$time_series} original time series data (input for model fitting),
#' \item \code{result$matrix_F} design matrix for fixed effects,
#' \item \code{result$matrix_V} design matrix for random effects,
#' \item \code{result$data_frame} data frame containing original time series, design matrices, grouping factors (for lme fitting),
#' \item \code{result$output$Box_test_season_resid} Box test of conditional residuals independence (accounts for time series seasonality),
#' available if \code{season_period} input argument is greater than zero,
#' \item \code{result$output$BoxLjung_test_season_resid} Box-Ljung test of conditional residuals independence (accounts for time series seasonality),
#' available if \code{season_period} input argument is greater than zero,
#' \item \code{result$Box_test_lag10_resid} Box test of conditional residuals independence (accounts for time series non-seasonality),
#' available if \code{season_period = 0},
#' \item \code{result$BoxLjung_test_lag10_resid} Box-Ljung test of conditional residuals independence (accounts for time series non-seasonality),
#' available if \code{season_period = 0},
#' \item \code{result$diagnostic_plots_names} vector of available diagnostic plots names, these can be displayed by function \code{drawDiagPlots}.
#' }
#' For CLRM (using \code{lm}) and also for other fitting functions than \code{lme} the output looks slightly different.
#' If set the diagnostic plots could be drawn.
#'
#' @note Ver.: 13-Jan-2019 11:15:48.
#'
#' @example R/Examples/example_fitDiagFDSLRM.R
#'
#' @export
#'
fitDiagFDSLRM <- function(x, times, freq_mean, poly_trend_degree = 0, include_fixed_eff, freq_random, include_random_eff, fit_function = "lme", var_estim_method = "REML", season_period = 0, display_plots = FALSE) {

        # checking input parameters
        if(!(fit_function %in% c("lme", "mixed", "mmer"))) {
                stop("fit_function input must be specified as lme/mixed/mmer.")
        }
        if((fit_function == "mixed" || fit_function == "mmer") && missing(freq_random)) {
                stop("Frequencies for random effects are missing!
                     Only mixed effects models can be fitted by mixed() and mmer().")
        }
        if(missing(x)) {
                stop("Data missing!")
        }
        if(missing(freq_mean) && missing(freq_random)) {
                stop("Frequencies for mean value and for random component missing!")
        }
        # '%!in%' <- function(x, y)!('%in%'(x, y))
        # if(var_estim_method %!in% c("ML", "REML")) {
        #         stop("Set var_estim_method = ML or var_estim_method = REML.")
        # }
        if(poly_trend_degree %% 1 != 0 || poly_trend_degree < 0) {
                stop("poly_trend_degree must be a non-negative integer!")
        }
        if(season_period %% 1 != 0 || season_period < 0) {
                stop("season_period must be a non-negative integer!")
        }
        if(!missing(include_fixed_eff) && (length(include_fixed_eff) != 2*length(freq_mean) || !all(include_fixed_eff %in% c(0, 1)))) {
                stop("include_fixed_eff must be a vector of lenght = 2*length(freq_mean), consisting of 0 and 1.")
        }
        if(!missing(include_random_eff) && (length(include_random_eff) != 2*length(freq_random) || !all(include_random_eff %in% c(0, 1)))) {
                stop("include_random_eff must be a vector of lenght = 2*length(freq_random), consisting of 0 and 1.")
        }
        if(!missing(freq_mean) && missing(include_fixed_eff)) {
                include_fixed_eff <- rep(1, 2 * length(freq_mean))
        }
        if(!missing(freq_random) && missing(include_random_eff)) {
                include_random_eff <- rep(1, 2 * length(freq_random))
        }

        if(fit_function == "lme") {

                if(missing(freq_random) && !missing(freq_mean)) {

                        if(poly_trend_degree > 0) {
                                auxDTF <- data.frame(times)
                                if(poly_trend_degree > 1) {
                                        for(j in 2:poly_trend_degree) {
                                                auxDTF <- cbind(auxDTF, times^j)
                                        }
                                }
                                auxF <- makeF(times, freq_mean)
                                if((!missing(include_fixed_eff) || length(include_fixed_eff) > 0) && sum(include_fixed_eff) < length(include_fixed_eff)) {
                                        auxF <- as.matrix(auxF[,-(which(0 == include_fixed_eff)+1)])
                                }
                                F <- cbind(cbind(auxF[,1], auxDTF), auxF[,2:ncol(auxF)])
                                colnames(F) <- NULL
                                kk <- ncol(F)
                                nn <- length(times)
                                d <- data.frame(rep(0,nn))
                                for(i in 2:kk) {
                                        d <- cbind(d, F[,i])
                                }
                                d <- d[,-1]
                                d <- cbind(d, x)

                                names(d) <- c(paste(rep("f", kk-1), as.character(2:kk), sep = ""), "x")

                                fit <- lm(as.formula(paste("x~1+", paste(names(d)[1:(kk-1)], collapse = "+"))), data = d, model = TRUE, x = TRUE)

                                output <- list(fixed_effects = coefficients(fit),
                                               error_variance = (summary(fit)$sigma)^2,
                                               fitted_values = fitted(fit),
                                               raw_residuals = resid(fit),
                                               stand_residuals = MASS::stdres(fit),
                                               stud_residuals = MASS::studres(fit),
                                               non_const_err_var_test = car::ncvTest(fit),
                                               test_autocor_errors = durbinWatsonTest(fit),
                                               Box_test_raw_resid = Box.test(resid(fit)),
                                               BoxLjung_test_raw_resid = Box.test(resid(fit), type = "Ljung-Box"),
                                               ShapiroWilk_test_raw_resid = shapiro.test(resid(fit)),
                                               Box_tests_stud_resid = Box.test(resid(fit)),
                                               BoxLjung_test_stud_resid = Box.test(MASS::studres(fit), type = "Ljung-Box"),
                                               ShapiroWilk_test_stud_resid = shapiro.test(MASS::studres(fit)),
                                               fit_summary = summary(fit),
                                               lm_fit = fit,
                                               time_series = x,
                                               matrix_F = as.matrix(F),
                                               data_frame = d
                                )

                                if(season_period > 0) {
                                        if(2* season_period > nn / 5) {
                                                output$Box_test_season_raw_resid <- Box.test(resid(fit), lag = nn / 5)
                                                output$BoxLjung_test_season_raw_resid <- Box.test(resid(fit), lag = nn / 5, type = "Ljung-Box")
                                                output$Box_test_season_stud_resid <- Box.test(MASS::studres(fit), lag = nn / 5)
                                                output$BoxLjung_test_season_stud_resid <- Box.test(MASS::studres(fit), lag = nn / 5, type = "Ljung-Box")
                                        } else {
                                                output$Box_test_season_raw_resid <- Box.test(resid(fit), lag = 2 * season_period)
                                                output$BoxLjung_test_season_raw_resid <- Box.test(resid(fit), lag = 2 * season_period, type = "Ljung-Box")
                                                output$Box_test_season_stud_resid <- Box.test(MASS::studres(fit), lag = 2 * season_period)
                                                output$BoxLjung_test_season_stud_resid <- Box.test(MASS::studres(fit), lag = 2 * season_period, type = "Ljung-Box")
                                        }
                                } else {
                                        output$Box_test_lag10_raw_resid <- Box.test(resid(fit), lag = 10)
                                        output$BoxLjung_test_lag10_raw_resid <- Box.test(resid(fit), lag = 10, type = "Ljung-Box")
                                        output$Box_test_lag10_stud_resid <- Box.test(MASS::studres(fit), lag = 10)
                                        output$BoxLjung_test_lag10_stud_resid <- Box.test(MASS::studres(fit), lag = 10, type = "Ljung-Box")
                                }

                                if(display_plots == TRUE) {

                                        x_ts <- ts(x)
                                        fitted_values_ts <- ts(output$fitted_values)
                                        ts.plot(x_ts, fitted_values_ts, gpars = list(col = c("black", "blue"), lwd = 2,
                                                                                     xlab="Time", ylab="",main="Original time series vs fitted values"))
                                        legend("topleft", legend = c("time series", "fitted values"), col = c("black", "blue"),
                                               lwd = 2, cex = 0.5, bty = "n", bg="transparent")
                                        grid()
                                        plot(output$raw_residuals, type = "o", xlab = "Time", ylab = "Raw residuals", main = "")
                                        grid()
                                        plot(output$stud_residuals, type = "o", xlab = "Time", ylab = "Studentized residuals", main = "")
                                        grid()
                                        plot(fitted(fit), resid(fit), xlab = "Fitted values", ylab = "Raw residuals", main = "")
                                        grid()
                                        plot(fitted(fit), output$stud_residuals, xlab = "Fitted values", ylab = "Studentized residuals", main = "")
                                        grid()
                                        car::spreadLevelPlot(fit)
                                        grid()
                                        qqnorm(resid(fit), ylab = "Raw residuals")
                                        qqline(resid(fit))
                                        grid()
                                        qqnorm(output$stud_residuals, ylab = "Studentized residuals")
                                        qqline(output$stud_residuals)
                                        grid()
                                        acf(output$raw_residuals, main = "ACF of raw residuals")
                                        pacf(output$raw_residuals, main = "PACF of  raw residuals")
                                        hist(output$raw_residuals, main = "Distribution of raw residuals", xlab = "Raw residuals")
                                        acf(output$stud_residuals, main = "ACF of studentized residuals")
                                        pacf(output$stud_residuals, main = "PACF of  studentized residuals")
                                        hist(output$stud_residuals, main = "Distribution of studentized residuals", xlab = "Studentized residuals")
                                        cpgram(output$stud_residuals, main = "Cumulative periodogram of studentized residuals", ci.col = "red")

                                }

                                diagnostic_plots_names <- list(FittedTimeSeries = "FittedTimeSeries",
                                                               RawResid = "RawResid",
                                                               StdResid = "StdResid",
                                                               FittedValuesVsRawResid = "FittedValuesVsRawResid",
                                                               FittedValuesVsStdResid = "FittedValuesVsStdResid",
                                                               SpreadLevelPlot = "SpreadLevelPlot",
                                                               QQPlotRawResid = "QQPlotRawResid",
                                                               QQPlotStdResid = "QQPlotStdResid",
                                                               ACFRawResid = "ACFRawResid",
                                                               PACFRawResid = "PACFRawResid",
                                                               HistogramRawResid = "HistogramRawResid",
                                                               ACFStdResid = "ACFStdResid",
                                                               PACFStdResid = "PACFStdResid",
                                                               HistogramStdResid = "HistogramStdResid",
                                                               CumulatPeriodogCondResid = "CumulatPeriodogStdResid"
                                )

                                output$diagnostic_plots_names <- diagnostic_plots_names

                                return(output)
                        } else {
                                F <- makeF(times, freq_mean)
                                if((!missing(include_fixed_eff) || length(include_fixed_eff) > 0) && sum(include_fixed_eff) < length(include_fixed_eff)) {
                                        F <- as.matrix(F[,-(which(0 == include_fixed_eff)+1)])
                                }
                                kk <- ncol(F)
                                nn <- length(times)
                                d <- data.frame(rep(0,nn))
                                for(i in 2:kk) {
                                        d <- cbind(d, F[,i])
                                }
                                d <- d[,-1]
                                d <- cbind(d, x)
                                names(d) <- c(paste(rep("f", kk-1), as.character(2:kk), sep = ""), "x")

                                fit <- lm(as.formula(paste("x~1+", paste(names(d)[1:(kk-1)], collapse = "+"))), data = d, model = TRUE, x = TRUE)

                                output <- list(fixed_effects = coefficients(fit),
                                               error_variance = (summary(fit)$sigma)^2,
                                               fitted_values = fitted(fit),
                                               raw_residuals = resid(fit),
                                               stand_residuals = MASS::stdres(fit),
                                               stud_residuals = MASS::studres(fit),
                                               non_const_err_var_test = car::ncvTest(fit),
                                               test_autocor_errors = durbinWatsonTest(fit),
                                               Box_test_raw_resid = Box.test(resid(fit)),
                                               BoxLjung_test_raw_resid = Box.test(resid(fit), type = "Ljung-Box"),
                                               ShapiroWilk_test_raw_resid = shapiro.test(resid(fit)),
                                               Box_tests_stud_resid = Box.test(resid(fit)),
                                               BoxLjung_test_stud_resid = Box.test(MASS::studres(fit), type = "Ljung-Box"),
                                               ShapiroWilk_test_stud_resid = shapiro.test(MASS::studres(fit)),
                                               fit_summary = summary(fit),
                                               lm_fit = fit,
                                               time_series = x,
                                               matrix_F = as.matrix(F),
                                               data_frame = d
                                )

                                if(season_period > 0) {
                                        if(2* season_period > nn / 5) {
                                                output$Box_test_season_raw_resid <- Box.test(resid(fit), lag = nn / 5)
                                                output$BoxLjung_test_season_raw_resid <- Box.test(resid(fit), lag = nn / 5, type = "Ljung-Box")
                                                output$Box_test_season_stud_resid <- Box.test(MASS::studres(fit), lag = nn / 5)
                                                output$BoxLjung_test_season_stud_resid <- Box.test(MASS::studres(fit), lag = nn / 5, type = "Ljung-Box")
                                        } else {
                                                output$Box_test_season_raw_resid <- Box.test(resid(fit), lag = 2 * season_period)
                                                output$BoxLjung_test_season_raw_resid <- Box.test(resid(fit), lag = 2 * season_period, type = "Ljung-Box")
                                                output$Box_test_season_stud_resid <- Box.test(MASS::studres(fit), lag = 2 * season_period)
                                                output$BoxLjung_test_season_stud_resid <- Box.test(MASS::studres(fit), lag = 2 * season_period, type = "Ljung-Box")
                                        }
                                } else {
                                        output$Box_test_lag10_raw_resid <- Box.test(resid(fit), lag = 10)
                                        output$BoxLjung_test_lag10_raw_resid <- Box.test(resid(fit), lag = 10, type = "Ljung-Box")
                                        output$Box_test_lag10_stud_resid <- Box.test(MASS::studres(fit), lag = 10)
                                        output$BoxLjung_test_lag10_stud_resid <- Box.test(MASS::studres(fit), lag = 10, type = "Ljung-Box")
                                }

                                if(display_plots == TRUE) {

                                        x_ts <- ts(x)
                                        fitted_values_ts <- ts(output$fitted_values)
                                        ts.plot(x_ts, fitted_values_ts, gpars = list(col = c("black", "blue"), lwd = 2,
                                                                                     xlab="Time", ylab="",main="Original time series vs fitted values"))
                                        legend("topleft", legend = c("time series", "fitted values"), col = c("black", "blue"),
                                               lwd = 2, cex = 0.5, bty = "n", bg="transparent")
                                        grid()
                                        plot(output$raw_residuals, type = "o", xlab = "Time", ylab = "Raw residuals", main = "")
                                        grid()
                                        plot(output$stud_residuals, type = "o", xlab = "Time", ylab = "Studentized residuals", main = "")
                                        grid()
                                        plot(fitted(fit), resid(fit), xlab = "Fitted values", ylab = "Raw residuals", main = "")
                                        grid()
                                        plot(fitted(fit), output$stud_residuals, xlab = "Fitted values", ylab = "Studentized residuals", main = "")
                                        grid()
                                        car::spreadLevelPlot(fit)
                                        grid()
                                        qqnorm(resid(fit), ylab = "Raw residuals")
                                        qqline(resid(fit))
                                        grid()
                                        qqnorm(output$stud_residuals, ylab = "Studentized residuals")
                                        qqline(output$stud_residuals)
                                        grid()
                                        acf(output$raw_residuals, main = "ACF of raw residuals")
                                        pacf(output$raw_residuals, main = "PACF of  raw residuals")
                                        hist(output$raw_residuals, main = "Distribution of raw residuals", xlab = "Raw residuals")
                                        acf(output$stud_residuals, main = "ACF of studentized residuals")
                                        pacf(output$stud_residuals, main = "PACF of  studentized residuals")
                                        hist(output$stud_residuals, main = "Distribution of studentized residuals", xlab = "Studentized residuals")
                                        cpgram(output$stud_residuals, main = "Cumulative periodogram of studentized residuals", ci.col = "red")

                                }

                                diagnostic_plots_names <- list(FittedTimeSeries = "FittedTimeSeries",
                                                               RawResid = "RawResid",
                                                               StdResid = "StdResid",
                                                               FittedValuesVsRawResid = "FittedValuesVsRawResid",
                                                               FittedValuesVsStdResid = "FittedValuesVsStdResid",
                                                               SpreadLevelPlot = "SpreadLevelPlot",
                                                               QQPlotRawResid = "QQPlotRawResid",
                                                               QQPlotStdResid = "QQPlotStdResid",
                                                               ACFRawResid = "ACFRawResid",
                                                               PACFRawResid = "PACFRawResid",
                                                               HistogramRawResid = "HistogramRawResid",
                                                               ACFStdResid = "ACFStdResid",
                                                               PACFStdResid = "PACFStdResid",
                                                               HistogramStdResid = "HistogramStdResid",
                                                               CumulatPeriodogStudResid = "CumulatPeriodogStudResid"
                                )

                                output$diagnostic_plots_names <- diagnostic_plots_names

                                return(output)
                        }

                } else if (missing(freq_mean) && !missing(freq_random)) {
                        if(poly_trend_degree > 0) {
                                auxDTF <- data.frame(times)
                                if(poly_trend_degree > 1) {
                                        for(j in 2:poly_trend_degree) {
                                                auxDTF <- cbind(auxDTF, times^j)
                                        }
                                }
                                auxF <- as.matrix(rep(1, length(x)))
                                F <- cbind(auxF, auxDTF)
                                colnames(F) <- NULL
                                V <- makeV(times, freq_random)
                                if((!missing(include_random_eff) || length(include_random_eff) > 0) && sum(include_random_eff) < length(include_random_eff)) {
                                        V <- as.matrix(V[,-which(0 == include_random_eff)])
                                }
                                kk <- ncol(F)
                                ll <- ncol(V)
                                nn <- length(times)
                                g <- rep(1, nn)
                                d <- data.frame(g)
                                for(i in 2:kk) {
                                        d <- cbind(d, F[,i])
                                }
                                for(j in 1:ll) {
                                        d <- cbind(d, V[,j])
                                }
                                d <- cbind(d, x)
                                names(d) <- c("g", paste(rep("f", kk-1), as.character(2:kk), sep = ""),
                                              paste(rep("v", ll), as.character(1:ll), sep = ""), "x")

                                fit <- nlme::lme(fixed = as.formula(paste("x~1+", paste(names(d)[2:kk], collapse = "+"))),
                                                 random = list(g = pdDiag(as.formula(paste("~-1+", paste(names(d)[(kk+1):(kk+ll)], collapse = "+"))))), data = d, method = var_estim_method)

                                invisible(capture.output(SingerEtAl_resid_diag <- residdiag.nlme(fit, limit = 2, plotid = 1:6, d = d, kk = kk, ll = ll)))

                                output <- list(fixed_effects = fixed.effects(fit),
                                               random_effects = random.effects(fit),
                                               error_variance = fit$sigma ^ 2,
                                               rand_eff_variance = getVarCov(fit),
                                               marg_fitted_values = fitted(fit, level=0),
                                               cond_fitted_values = fitted(fit),
                                               marg_residuals = resid(fit, level=0),
                                               cond_residuals = resid(fit),
                                               norm_cond_resid = resid(fit, type="normalized"),
                                               pearson_resid = resid(fit, type="pearson"),
                                               Box_test = Box.test(resid(fit)),
                                               BoxLjung_test = Box.test(resid(fit), type = "Ljung-Box"),
                                               ShapiroWilk_test_norm_cond_resid = shapiro.test(resid(fit, type="normalized")),
                                               ShapiroWilk_test_stand_least_conf_resid = shapiro.test(SingerEtAl_resid_diag$least.confounded.residuals),
                                               lme_fit = fit,
                                               fit_summary = summary(fit),
                                               diag_resid = SingerEtAl_resid_diag,
                                               time_series = x,
                                               matrix_F = as.matrix(F),
                                               matrix_V = as.matrix(V),
                                               data_frame = d
                                )

                                if(season_period > 0) {
                                        if(2* season_period > nn / 5) {
                                                output$Box_test_season_resid <- Box.test(resid(fit), lag = nn / 5)
                                                output$BoxLjung_test_season_resid <- Box.test(resid(fit), lag = nn / 5, type = "Ljung-Box")
                                        } else {
                                                output$Box_test_season_resid <- Box.test(resid(fit), lag = 2 * season_period)
                                                output$BoxLjung_test_season_resid <- Box.test(resid(fit), lag = 2 * season_period, type = "Ljung-Box")
                                        }
                                } else {
                                        output$Box_test_lag10_resid <- Box.test(resid(fit), lag = 10)
                                        output$BoxLjung_test_lag10_resid <- Box.test(resid(fit), lag = 10, type = "Ljung-Box")
                                }

                                if(display_plots == TRUE) {

                                        par(mfrow=c(1,1))
                                        x_ts <- ts(x)
                                        fitted_values_ts <- ts(output$cond_fitted_values)
                                        trend_ts <- ts(as.vector(output$matrix_F%*%output$fixed_effects))
                                        ts.plot(x_ts, fitted_values_ts, trend_ts, gpars = list(col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = 2,
                                                                                               xlab="Time", ylab="",main="Original time series vs fitted values"))
                                        legend("topleft", legend = c("time series", "fitted values", "estimated trend"), col = c("black", "blue", "red"),
                                               lty = c(1, 1, 2), lwd = 2, cex = 0.5, bty = "n", bg="transparent")
                                        grid()
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(2), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(5), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(6), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                        plot(output$diag_resid$std.marginal.residuals[,2], type = "o", ylab = "Standardized marginal residuals")
                                        abline(0,0)
                                        grid()
                                        plot(output$diag_resid$std.conditional.residuals[,2], type = "o", ylab = "Standardized conditional residuals")
                                        abline(0,0)
                                        grid()
                                        acf(output$cond_residuals, main = "ACF of conditional residuals")
                                        pacf(output$cond_residuals, main = "PACF of conditional residuals")
                                        hist(output$norm_cond_resid, main = "Distribution of norm. cond. residuals", xlab = "Norm. cond. residuals")
                                        hist(output$diag_resid$least.confounded.residuals, main = "Distribution of stand. least conf. residuals", xlab = "Stand. least conf. residuals")
                                        cpgram(output$cond_residuals, main = "Cumulative periodogram of conditional residuals", ci.col = "red")
                                }

                                diagnostic_plots_names <- list(StdMarginalResidVsFittedValues = "StdMarginalResidVsFittedValues",
                                                               StdConditionalResidVsFittedValues = "StdConditionalResidVsFittedValues",
                                                               NormalQQPlotStdLeastConfResid = "NormalQQPlotStdLeastConfResid",
                                                               FittedTimeSeries = "FittedTimeSeries",
                                                               StdMarginalResid = "StdMarginalResid",
                                                               StdConditionalResid = "StdConditionalResid",
                                                               ACFCondtResid = "ACFCondtResid",
                                                               PACFCondResid = "PACFCondResid",
                                                               HistNormCondResid = "HistNormCondResid",
                                                               HistLeastConfResid = "HistLeastConfResid",
                                                               CumulatPeriodogStudResid = "CumulatPeriodogStudResid")

                                output$diagnostic_plots_names <- diagnostic_plots_names

                                return(output)

                        } else {
                                F <- as.matrix(rep(1, length(x)))
                                V <- makeV(times, freq_random)
                                if((!missing(include_random_eff) || length(include_random_eff) > 0) && sum(include_random_eff) < length(include_random_eff)) {
                                        V <- as.matrix(V[,-which(0 == include_random_eff)])
                                }
                                kk <- ncol(F)
                                ll <- ncol(V)
                                nn <- length(times)
                                g <- rep(1, nn)
                                d <- data.frame(g)
                                for(j in 1:ll) {
                                        d <- cbind(d, V[,j])
                                }
                                d <- cbind(d, x)
                                names(d) <- c("g", paste(rep("v", ll), as.character(1:ll), sep = ""), "x")

                                fit <- nlme::lme(fixed = x~1, random = list(g = pdDiag(as.formula(paste("~-1+", paste(names(d)[2:(ll)], collapse = "+"))))), data = d, method = var_estim_method)

                                invisible(capture.output(SingerEtAl_resid_diag <- residdiag.nlme(fit, limit = 2, plotid = 1:6, d = d, kk = kk, ll = ll)))

                                output <- list(fixed_effects = fixed.effects(fit),
                                               random_effects = random.effects(fit),
                                               error_variance = fit$sigma ^ 2,
                                               rand_eff_variance = getVarCov(fit),
                                               marg_fitted_values = fitted(fit, level=0),
                                               cond_fitted_values = fitted(fit),
                                               marg_residuals = resid(fit, level=0),
                                               cond_residuals = resid(fit),
                                               norm_cond_resid = resid(fit, type="normalized"),
                                               pearson_resid = resid(fit, type="pearson"),
                                               Box_test = Box.test(resid(fit)),
                                               BoxLjung_test = Box.test(resid(fit), type = "Ljung-Box"),
                                               ShapiroWilk_test_norm_cond_resid = shapiro.test(resid(fit, type="normalized")),
                                               ShapiroWilk_test_stand_least_conf_resid = shapiro.test(SingerEtAl_resid_diag$least.confounded.residuals),
                                               lme_fit = fit,
                                               fit_summary = summary(fit),
                                               diag_resid = SingerEtAl_resid_diag,
                                               time_series = x,
                                               matrix_F = as.matrix(F),
                                               matrix_V = as.matrix(V),
                                               data_frame = d
                                )

                                if(season_period > 0) {
                                        if(2* season_period > nn / 5) {
                                                output$Box_test_season_resid <- Box.test(resid(fit), lag = nn / 5)
                                                output$BoxLjung_test_season_resid <- Box.test(resid(fit), lag = nn / 5, type = "Ljung-Box")
                                        } else {
                                                output$Box_test_season_resid <- Box.test(resid(fit), lag = 2 * season_period)
                                                output$BoxLjung_test_season_resid <- Box.test(resid(fit), lag = 2 * season_period, type = "Ljung-Box")
                                        }
                                } else {
                                        output$Box_test_lag10_resid <- Box.test(resid(fit), lag = 10)
                                        output$BoxLjung_test_lag10_resid <- Box.test(resid(fit), lag = 10, type = "Ljung-Box")
                                }

                                if(display_plots == TRUE) {

                                        par(mfrow=c(1,1))
                                        x_ts <- ts(x)
                                        fitted_values_ts <- ts(output$cond_fitted_values)
                                        trend_ts <- ts(as.vector(output$matrix_F%*%output$fixed_effects))
                                        ts.plot(x_ts, fitted_values_ts, trend_ts, gpars = list(col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = 2,
                                                                                               xlab="Time", ylab="",main="Original time series vs fitted values"))
                                        legend("topleft", legend = c("time series", "fitted values", "estimated trend"), col = c("black", "blue", "red"),
                                               lty = c(1, 1, 2), lwd = 2, cex = 0.5, bty = "n", bg="transparent")
                                        grid()
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(2), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(5), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(6), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                        plot(output$diag_resid$std.marginal.residuals[,2], type = "o", ylab = "Standardized marginal residuals")
                                        abline(0,0)
                                        grid()
                                        plot(output$diag_resid$std.conditional.residuals[,2], type = "o", ylab = "Standardized conditional residuals")
                                        abline(0,0)
                                        grid()
                                        acf(output$cond_residuals, main = "ACF of conditional residuals")
                                        pacf(output$cond_residuals, main = "PACF of conditional residuals")
                                        hist(output$norm_cond_resid, main = "Distribution of norm. cond. residuals", xlab = "Norm. cond. residuals")
                                        hist(output$diag_resid$least.confounded.residuals, main = "Distribution of stand. least conf. residuals", xlab = "Stand. least conf. residuals")
                                        cpgram(output$cond_residuals, main = "Cumulative periodogram of conditional residuals", ci.col = "red")
                                }

                                diagnostic_plots_names <- list(StdMarginalResidVsFittedValues = "StdMarginalResidVsFittedValues",
                                                               StdConditionalResidVsFittedValues = "StdConditionalResidVsFittedValues",
                                                               NormalQQPlotStdLeastConfResid = "NormalQQPlotStdLeastConfResid",
                                                               FittedTimeSeries = "FittedTimeSeries",
                                                               StdMarginalResid = "StdMarginalResid",
                                                               StdConditionalResid = "StdConditionalResid",
                                                               ACFCondtResid = "ACFCondtResid",
                                                               PACFCondResid = "PACFCondResid",
                                                               HistNormCondResid = "HistNormCondResid",
                                                               HistLeastConfResid = "HistLeastConfResid",
                                                               CumulatPeriodogCondResid = "CumulatPeriodogCondResid")

                                output$diagnostic_plots_names <- diagnostic_plots_names

                                return(output)

                        }

                } else {
                        if(poly_trend_degree > 0) {
                                auxDTF <- data.frame(times)
                                if(poly_trend_degree > 1) {
                                        for(j in 2:poly_trend_degree) {
                                                auxDTF <- cbind(auxDTF, times^j)
                                        }
                                }
                                auxF <- makeF(times, freq_mean)
                                if((!missing(include_fixed_eff) || length(include_fixed_eff) > 0) && sum(include_fixed_eff) < length(include_fixed_eff)) {
                                        auxF <- as.matrix(auxF[,-(which(0 == include_fixed_eff)+1)])
                                }
                                F <- cbind(cbind(auxF[,1], auxDTF), auxF[,2:ncol(auxF)])
                                colnames(F) <- NULL
                                V <- makeV(times, freq_random)
                                if((!missing(include_random_eff) || length(include_random_eff) > 0) && sum(include_random_eff) < length(include_random_eff)) {
                                        V <- as.matrix(V[,-which(0 == include_random_eff)])
                                }
                                kk <- ncol(F)
                                ll <- ncol(V)
                                nn <- length(times)
                                g <- rep(1, nn)
                                d <- data.frame(g)
                                for(i in 2:kk) {
                                        d <- cbind(d, F[,i])
                                }
                                for(j in 1:ll) {
                                        d <- cbind(d, V[,j])
                                }
                                d <- cbind(d, x)
                                names(d) <- c("g", paste(rep("f", kk-1), as.character(2:kk), sep = ""), paste(rep("v", ll), as.character(1:ll), sep = ""), "x")

                                fit <- nlme::lme(fixed = as.formula(paste("x~", paste(names(d)[2:kk], collapse = "+"))),
                                                 random = list(g = pdDiag(as.formula(paste("~-1+", paste(names(d)[(kk+1):(kk+ll)], collapse = "+"))))), data = d, method = var_estim_method)

                                invisible(capture.output(SingerEtAl_resid_diag <- residdiag.nlme(fit, limit = 2, plotid = c(2,5,6), d = d, kk = kk, ll = ll)))

                                output <- list(fixed_effects = fixed.effects(fit),
                                               random_effects = random.effects(fit),
                                               error_variance = fit$sigma ^ 2,
                                               rand_eff_variance = getVarCov(fit),
                                               marg_fitted_values = fitted(fit, level=0),
                                               cond_fitted_values = fitted(fit),
                                               marg_residuals = resid(fit, level=0),
                                               cond_residuals = resid(fit),
                                               norm_cond_resid = resid(fit, type="normalized"),
                                               pearson_resid = resid(fit, type="pearson"),
                                               Box_test = Box.test(resid(fit)),
                                               BoxLjung_test = Box.test(resid(fit), type = "Ljung-Box"),
                                               ShapiroWilk_test_norm_cond_resid = shapiro.test(resid(fit, type="normalized")),
                                               ShapiroWilk_test_stand_least_conf_resid = shapiro.test(SingerEtAl_resid_diag$least.confounded.residuals),
                                               lme_fit = fit,
                                               fit_summary = summary(fit),
                                               diag_resid = SingerEtAl_resid_diag,
                                               time_series = x,
                                               matrix_F = as.matrix(F),
                                               matrix_V = as.matrix(V),
                                               data_frame = d
                                )

                                if(season_period > 0) {
                                        if(2* season_period > nn / 5) {
                                                output$Box_test_season_resid <- Box.test(resid(fit), lag = nn / 5)
                                                output$BoxLjung_test_season_resid <- Box.test(resid(fit), lag = nn / 5, type = "Ljung-Box")
                                        } else {
                                                output$Box_test_season_resid <- Box.test(resid(fit), lag = 2 * season_period)
                                                output$BoxLjung_test_season_resid <- Box.test(resid(fit), lag = 2 * season_period, type = "Ljung-Box")
                                        }
                                } else {
                                        output$Box_test_lag10_resid <- Box.test(resid(fit), lag = 10)
                                        output$BoxLjung_test_lag10_resid <- Box.test(resid(fit), lag = 10, type = "Ljung-Box")
                                }

                                if(display_plots == TRUE) {

                                        par(mfrow=c(1,1))
                                        x_ts <- ts(x)
                                        fitted_values_ts <- ts(output$cond_fitted_values)
                                        trend_ts <- ts(as.vector(output$matrix_F%*%output$fixed_effects))
                                        ts.plot(x_ts, fitted_values_ts, trend_ts, gpars = list(col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = 2,
                                                                                               xlab="Time", ylab="",main="Original time series vs fitted values"))
                                        legend("topleft", legend = c("time series", "fitted values", "estimated trend"), col = c("black", "blue", "red"),
                                               lty = c(1, 1, 2), lwd = 2, cex = 0.5, bty = "n", bg="transparent")
                                        grid()
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(2), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(5), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(6), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                        plot(output$diag_resid$std.marginal.residuals[,2], type = "o", ylab = "Standardized marginal residuals")
                                        abline(0,0)
                                        grid()
                                        plot(output$diag_resid$std.conditional.residuals[,2], type = "o", ylab = "Standardized conditional residuals")
                                        abline(0,0)
                                        grid()
                                        acf(output$cond_residuals, main = "ACF of conditional residuals")
                                        pacf(output$cond_residuals, main = "PACF of conditional residuals")
                                        hist(output$norm_cond_resid, main = "Distribution of norm. cond. residuals", xlab = "Norm. cond. residuals")
                                        hist(output$diag_resid$least.confounded.residuals, main = "Distribution of stand. least conf. residuals", xlab = "Stand. least conf. residuals")
                                        cpgram(output$cond_residuals, main = "Cumulative periodogram of conditional residuals", ci.col = "red")

                                }

                                diagnostic_plots_names <- list(StdMarginalResidVsFittedValues = "StdMarginalResidVsFittedValues",
                                                               StdConditionalResidVsFittedValues = "StdConditionalResidVsFittedValues",
                                                               NormalQQPlotStdLeastConfResid = "NormalQQPlotStdLeastConfResid",
                                                               FittedTimeSeries = "FittedTimeSeries",
                                                               StdMarginalResid = "StdMarginalResid",
                                                               StdConditionalResid = "StdConditionalResid",
                                                               ACFCondtResid = "ACFCondtResid",
                                                               PACFCondResid = "PACFCondResid",
                                                               HistNormCondResid = "HistNormCondResid",
                                                               HistLeastConfResid = "HistLeastConfResid",
                                                               CumulatPeriodogCondResid = "CumulatPeriodogCondResid")

                                output$diagnostic_plots_names <- diagnostic_plots_names

                                return(output)

                        } else {
                                F <- makeF(times, freq_mean)
                                colnames(F) <- NULL
                                if((!missing(include_fixed_eff) || length(include_fixed_eff) > 0) && sum(include_fixed_eff) < length(include_fixed_eff)) {
                                        F <- as.matrix(F[,-(which(0 == include_fixed_eff)+1)])
                                }
                                V <- makeV(times, freq_random)
                                if((!missing(include_random_eff) || length(include_random_eff) > 0) && sum(include_random_eff) < length(include_random_eff)) {
                                        V <- as.matrix(V[,-which(0 == include_random_eff)])
                                }
                                kk <- ncol(F)
                                ll <- ncol(V)
                                nn <- length(times)
                                g <- rep(1, nn)
                                d <- data.frame(g)
                                for(i in 2:kk) {
                                        d <- cbind(d, F[,i])
                                }
                                for(j in 1:ll) {
                                        d <- cbind(d, V[,j])
                                }
                                d <- cbind(d, x)
                                names(d) <- c("g", paste(rep("f", kk-1), as.character(2:kk), sep = ""), paste(rep("v", ll), as.character(1:ll), sep = ""), "x")

                                fit <- nlme::lme(fixed = as.formula(paste("x~", paste(names(d)[2:kk], collapse = "+"))), random = list(g = pdDiag(as.formula(paste("~-1+", paste(names(d)[(kk+1):(kk+ll)], collapse = "+"))))), data = d, method = var_estim_method)

                                invisible(capture.output(SingerEtAl_resid_diag <- residdiag.nlme(fit, limit = 2, plotid = c(2,5,6), d = d, kk = kk, ll = ll)))

                                output <- list(fixed_effects = fixed.effects(fit),
                                               random_effects = random.effects(fit),
                                               error_variance = fit$sigma ^ 2,
                                               rand_eff_variance = getVarCov(fit),
                                               marg_fitted_values = fitted(fit, level=0),
                                               cond_fitted_values = fitted(fit),
                                               marg_residuals = resid(fit, level=0),
                                               cond_residuals = resid(fit),
                                               norm_cond_resid = resid(fit, type="normalized"),
                                               pearson_resid = resid(fit, type="pearson"),
                                               Box_test = Box.test(resid(fit)),
                                               BoxLjung_test = Box.test(resid(fit), type = "Ljung-Box"),
                                               ShapiroWilk_test_norm_cond_resid = shapiro.test(resid(fit, type="normalized")),
                                               ShapiroWilk_test_stand_least_conf_resid = shapiro.test(SingerEtAl_resid_diag$least.confounded.residuals),
                                               lme_fit = fit,
                                               fit_summary = summary(fit),
                                               diag_resid = SingerEtAl_resid_diag,
                                               time_series = x,
                                               matrix_F = as.matrix(F),
                                               matrix_V = as.matrix(V),
                                               data_frame = d
                                )

                                if(season_period > 0) {
                                        if(2* season_period > nn / 5) {
                                                output$Box_test_season_resid <- Box.test(resid(fit), lag = nn / 5)
                                                output$BoxLjung_test_season_resid <- Box.test(resid(fit), lag = nn / 5, type = "Ljung-Box")
                                        } else {
                                                output$Box_test_season_resid <- Box.test(resid(fit), lag = 2 * season_period)
                                                output$BoxLjung_test_season_resid <- Box.test(resid(fit), lag = 2 * season_period, type = "Ljung-Box")
                                        }
                                } else {
                                        output$Box_test_lag10_resid <- Box.test(resid(fit), lag = 10)
                                        output$BoxLjung_test_lag10_resid <- Box.test(resid(fit), lag = 10, type = "Ljung-Box")
                                }

                                if(display_plots == TRUE) {

                                        par(mfrow=c(1,1))
                                        x_ts <- ts(x)
                                        fitted_values_ts <- ts(output$cond_fitted_values)
                                        trend_ts <- ts(as.vector(output$matrix_F%*%output$fixed_effects))
                                        ts.plot(x_ts, fitted_values_ts, trend_ts, gpars = list(col = c("black", "blue", "red"), lty = c(1, 1, 2), lwd = 2,
                                                                                               xlab="Time", ylab="",main="Original time series vs fitted values"))
                                        legend("topleft", legend = c("time series", "fitted values", "estimated trend"), col = c("black", "blue", "red"),
                                               lty = c(1, 1, 2), lwd = 2, cex = 0.5, bty = "n", bg="transparent")
                                        grid()
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(2), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(5), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                        SingerEtAl_justPlots256(output$lme_fit, limit = 2, plotid = c(6), d = output$data_frame, kk = ncol(output$matrix_F), ll = ncol(output$matrix_V), histogram = TRUE)
                                        par(mfrow=c(1,1))
                                        plot(output$diag_resid$std.marginal.residuals[,2], type = "o", ylab = "Standardized marginal residuals")
                                        abline(0,0)
                                        grid()
                                        plot(output$diag_resid$std.conditional.residuals[,2], type = "o", ylab = "Standardized conditional residuals")
                                        abline(0,0)
                                        grid()
                                        acf(output$cond_residuals, main = "ACF of conditional residuals")
                                        pacf(output$cond_residuals, main = "PACF of conditional residuals")
                                        hist(output$norm_cond_resid, main = "Distribution of norm. cond. residuals", xlab = "Norm. cond. residuals")
                                        hist(output$diag_resid$least.confounded.residuals, main = "Distribution of stand. least conf. residuals", xlab = "Stand. least conf. residuals")
                                        cpgram(output$cond_residuals, main = "Cumulative periodogram of conditional residuals", ci.col = "red")

                                }

                                diagnostic_plots_names <- list(StdMarginalResidVsFittedValues = "StdMarginalResidVsFittedValues",
                                                               StdConditionalResidVsFittedValues = "StdConditionalResidVsFittedValues",
                                                               NormalQQPlotStdLeastConfResid = "NormalQQPlotStdLeastConfResid",
                                                               FittedTimeSeries = "FittedTimeSeries",
                                                               StdMarginalResid = "StdMarginalResid",
                                                               StdConditionalResid = "StdConditionalResid",
                                                               ACFCondtResid = "ACFCondtResid",
                                                               PACFCondResid = "PACFCondResid",
                                                               HistNormCondResid = "HistNormCondResid",
                                                               HistLeastConfResid = "HistLeastConfResid",
                                                               CumulatPeriodogCondResid = "CumulatPeriodogCondResid")

                                output$diagnostic_plots_names <- diagnostic_plots_names

                                return(output)

                        }
                }

        } else if(fit_function == "mixed") {

                if(!(var_estim_method %in% c("ML", "REML", "MINQEI", "MINQEUI"))) {
                        stop("var_estim_method for mixed() must be ML / REML / MINQEI / MINQEUI.")
                }

                n <- length(x)
                F <- matrix()

                if(missing(freq_mean)) {
                        F <- rep(1, n)
                } else {
                        F <- makeF(times, freq_mean)
                        if((!missing(include_fixed_eff) || length(include_fixed_eff) > 0) && sum(include_fixed_eff) < length(include_fixed_eff)) {
                                F <- as.matrix(F[,-(which(0 == include_fixed_eff)+1)])
                        }
                }

                V <- makeV(times, freq_random)
                if((!missing(include_random_eff) || length(include_random_eff) > 0) && sum(include_random_eff) < length(include_random_eff)) {
                        V <- as.matrix(V[,-which(0 == include_random_eff)])
                }
                l <- ncol(V)
                variances_start <- rep(1, l+1)
                dimensions <- rep(1, l)
                var_est_met <- NULL
                if(var_estim_method == "ML") {
                        var_est_met <- 1
                } else if (var_estim_method == "REML") {
                        var_est_met <- 2
                } else if (var_estim_method == "MINQEI") {
                        var_est_met <- 3
                } else {
                        var_est_met <- 4
                }

                fit <- mixed(x, F, V, dimensions, variances_start, var_est_met)

                output <- list(fixed_effects = fit$b,
                               random_effects = fit$u,
                               error_variance = fit$s2[l+1],
                               rand_eff_variance = fit$s2[1:l],
                               marg_fitted_values = F %*% fit$b,
                               cond_fitted_values = F %*% fit$b + V %*% fit$u,
                               marg_residuals = x - F %*% fit$b,
                               cond_residuals = x - F %*% fit$b - V %*% fit$u,
                               Box_test = Box.test(x - F %*% fit$b - V %*% fit$u),
                               BoxLjung_test = Box.test(x - F %*% fit$b - V %*% fit$u, type = "Ljung-Box"),
                               ShapiroWilk_test_cond_resid = shapiro.test(x - F %*% fit$b - V %*% fit$u),
                               time_series = x,
                               matrix_F = F,
                               matrix_V = V,
                               mixed_fit = fit)

                if(season_period > 0) {
                        if(2* season_period > n / 5) {
                                output$Box_test_season_resid <- Box.test(output$cond_residuals, lag = n / 5)
                                output$BoxLjung_test_season_resid <- Box.test(output$cond_residuals, lag = n / 5, type = "Ljung-Box")
                        } else {
                                output$Box_test_season_resid <- Box.test(output$cond_residuals, lag = 2 * season_period)
                                output$BoxLjung_test_season_resid <- Box.test(output$cond_residuals, lag = 2 * season_period, type = "Ljung-Box")
                        }
                } else {
                        output$Box_test_lag10_resid <- Box.test(output$cond_residuals, lag = 10)
                        output$BoxLjung_test_lag10_resid <- Box.test(output$cond_residuals, lag = 10, type = "Ljung-Box")
                }

                if(display_plots == TRUE) {

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
                        cpgram(output$cond_residuals, main = "Cumulative periodogram of conditional residuals", ci.col = "red")

                }

                diagnostic_plots_names <- list(MarginalResidVsFittedValues = "MarginalResidVsFittedValues",
                                               ConditionalResidVsFittedValues = "ConditionalResidVsFittedValues",
                                               FittedTimeSeries = "FittedTimeSeries",
                                               MarginalResid = "MarginalResid",
                                               ConditionalResid = "ConditionalResid",
                                               ACFCondtResid = "ACFCondtResid",
                                               PACFCondResid = "PACFCondResid",
                                               HistCondResid = "HistCondResid",
                                               HistMargResid = "HistMargResid",
                                               CumulatPeriodogCondResid = "CumulatPeriodogCondResid")

                output$diagnostic_plots_names <- diagnostic_plots_names

                return(output)

        } else {

                # if(!(var_estim_method %in% c("ML", "REML"))) {
                #         stop("var_estim_method for mmer() must be ML or REML.")
                # }
                #
                n <- length(times)
                F <- matrix()

                if(missing(freq_mean)) {
                        F <- rep(1, n)
                } else {
                        F <- makeF(times, freq_mean)
                        if((!missing(include_fixed_eff) || length(include_fixed_eff) > 0) && sum(include_fixed_eff) < length(include_fixed_eff)) {
                                F <- as.matrix(F[,-(which(0 == include_fixed_eff)+1)])
                        }
                }

                k <- ncol(F)

                V <- makeV(times, freq_random)
                if((!missing(include_random_eff) || length(include_random_eff) > 0) && sum(include_random_eff) < length(include_random_eff)) {
                        V <- as.matrix(V[,-which(0 == include_random_eff)])
                }
                l <- ncol(V)
                #
                # ETA <- list()
                # for(i in 1:l) {
                #         Zi <- as.matrix(V[,i])
                #         ETA[[i]] <- list(Z=Zi, K=diag(1))
                # }
                #
                # fit <- list()
                #
                # if(var_estim_method == "ML") {
                #
                #         fit <- sommer::mmer(Y = x, X = F, Z = ETA, REML = FALSE)
                #
                # } else {
                #
                #         fit <- sommer::mmer(Y = x, X = F, Z = ETA)
                #
                # }

                d <- data.frame(F[,1])
                for(i in 2:k) {
                        d <- cbind(d, F[,i])
                }
                for(j in 1:l) {
                        d <- cbind(d, V[,j])
                }
                d <- cbind(d, x)
                names(d) <- c(paste(rep("f", k), as.character(1:k), sep = ""), paste(rep("v", l), as.character(1:l), sep = ""), "x")
                names_aux <- paste(paste(rep("vs(ds(v", l), as.character(1:l), sep = ""), rep("),1)", l), sep = "")

                fit <- mmer(fixed = as.formula(paste("x~", paste(names(d)[1:k], collapse = "+"))), random = as.formula(paste("~", paste(names_aux, collapse = "+"))), data = d, verbose = FALSE)

                output <- list(fixed_effects = as.vector(fit$Beta$Estimate),
                               random_effects = as.vector(unlist(fit$U)),
                               error_variance = as.vector(unlist(fit$sigma))[l+1],
                               rand_eff_variance = as.vector(unlist(fit$sigma))[1:l],
                               marg_fitted_values = F %*% as.vector(fit$Beta$Estimate),
                               cond_fitted_values = F %*% as.vector(fit$Beta$Estimate) + V %*% as.vector(unlist(fit$U)),
                               marg_residuals = as.vector(fit$residuals),
                               cond_residuals = x - (F %*% as.vector(fit$Beta$Estimate) + V %*% as.vector(unlist(fit$U))),
                               Box_test = Box.test(x - (F %*% as.vector(fit$Beta$Estimate) + V %*% as.vector(unlist(fit$U)))),
                               BoxLjung_test = Box.test(x - (F %*% as.vector(fit$Beta$Estimate) + V %*% as.vector(unlist(fit$U))), type = "Ljung-Box"),
                               ShapiroWilk_test_cond_resid = shapiro.test(x - (F %*% as.vector(fit$Beta$Estimate) + V %*% as.vector(unlist(fit$U)))),
                               AIC = fit$AIC,
                               BIC = fit$BIC,
                               time_series = x,
                               matrix_F = F,
                               matrix_V = V,
                               fit_summary = summary(fit),
                               mmer_fit = fit)

                if(season_period > 0) {
                        if(2* season_period > n / 5) {
                                output$Box_test_season_resid <- Box.test(output$cond_residuals, lag = n / 5)
                                output$BoxLjung_test_season_resid <- Box.test(output$cond_residuals, lag = n / 5, type = "Ljung-Box")
                        } else {
                                output$Box_test_season_resid <- Box.test(output$cond_residuals, lag = 2 * season_period)
                                output$BoxLjung_test_season_resid <- Box.test(output$cond_residuals, lag = 2 * season_period, type = "Ljung-Box")
                        }
                } else {
                        output$Box_test_lag10_resid <- Box.test(output$cond_residuals, lag = 10)
                        output$BoxLjung_test_lag10_resid <- Box.test(output$cond_residuals, lag = 10, type = "Ljung-Box")
                }

                if(display_plots == TRUE) {

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
                        cpgram(output$cond_residuals, main = "Cumulative periodogram of conditional residuals", ci.col = "red")

                }

                diagnostic_plots_names <- list(MarginalResidVsFittedValues = "MarginalResidVsFittedValues",
                                               ConditionalResidVsFittedValues = "ConditionalResidVsFittedValues",
                                               FittedTimeSeries = "FittedTimeSeries",
                                               MarginalResid = "MarginalResid",
                                               ConditionalResid = "ConditionalResid",
                                               ACFCondtResid = "ACFCondtResid",
                                               PACFCondResid = "PACFCondResid",
                                               HistCondResid = "HistCondResid",
                                               HistMargResid = "HistMargResid",
                                               CumulatPeriodogCondResid = "CumulatPeriodogCondResid")

                output$diagnostic_plots_names <- diagnostic_plots_names

                return(output)
        }

}
