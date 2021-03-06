## EXAMPLE 1 (LMM)
# load all necessary R packages and functions to fit and make diagnostic of FDSLRM
initialFDSLRM()
# the following data set was adapted from Hyndman's package fpp2 (Hyndman 2018)
# data represent total quarterly visitor nights (in millions) from 1998-2016 in one of the regions
# of Australia - inner zone of Victoria state
# the number of time series observations is  n = 76
# data visualization
autoplot(window(visnights)[,15], ylab = "Visnights")
# time series realization, observation times
dt <- as.numeric((window(visnights)[,15]))
t <- 1:length(dt)
# periodogram to identify the most significant frequencies to build the model (design matrices)
periodo <- spec.pgram(dt, log="no")
# according to periodogram we choose the frequencies 1/80, 2/80 for trend
# and the frequencies 20/80, 40/80 for random component
# fit the model with lme function
output_lme <- fitDiagFDSLRM(dt, t, c(1/80, 2/80), include_fixed_eff = c(1,0,0,1),
                          freq_random = c(20/80, 40/80), include_random_eff = c(1,1,1,0),
                          poly_trend_degree = 0, season_period = 4)
# draw diagnostic plots
drawDiagPlots("all", output_lme)
drawDiagPlots(output_lme$diagnostic_plots_names$FittedTimeSeries, output_lme)
# show summary of model fit
str(output_lme$fit_summary)
# fit the model with mixed function
output_mixed <- fitDiagFDSLRM(dt, t, c(1/80, 2/80), include_fixed_eff = c(1,0,0,1),
                            freq_random = c(20/80, 40/80), include_random_eff = c(1,1,1,0),
                            var_estim_method = "MINQEI",
                            fit_function = "mixed", season_period = 4)
# draw diagnostic plots
drawDiagPlots("all", output_mixed)
drawDiagPlots(output_mixed$diagnostic_plots_names$CumulatPeriodogCondResid, output_lme)
# show summary of model fit
str(output_mixed$mixed_fit)

# fit the model with mmer function
output_mmer <- fitDiagFDSLRM(dt, t, c(1/80, 2/80), include_fixed_eff = c(1,0,0,1),
                            freq_random = c(20/80, 40/80), include_random_eff = c(1,1,1,0),
                            fit_function = "mmer", season_period = 4)
# draw diagnostic plots
drawDiagPlots("all", output_mmer)
drawDiagPlots(output_mmer$diagnostic_plots_names$CumulatPeriodogCondResid, output_mmer)
# show summary of model fit
str(output_mmer$mmer_fit)

## EXAMPLE 2 (CLRM)
# load all necessary R packages and functions to fit and make diagnostic of FDSLRM
initialFDSLRM()
# the following data set was adapted from Hyndman's package fpp2 (Hyndman 2018)
# data represent total quarterly visitor nights (in millions) from 1998-2016 in one of the regions
# of Australia - inner zone of Victoria state
# the number of time series observations is  n = 76
# data visualization
autoplot(window(visnights)[,15], ylab = "Visnights")
# time series realization, observation times
dt <- as.numeric((window(visnights)[,15]))
t <- 1:length(dt)
# periodogram to identify the most significant frequencies to build the model (design matrices)
periodo <- spec.pgram(dt, log="no")
# according to periodogram we choose the frequencies 1/80, 2/80 for trend
# and the frequencies 20/80, 40/80 for random component
# fit the model with lme function
output_lm <- fitDiagFDSLRM(dt, t, c(1/80, 2/80, 20/80, 40/80), include_fixed_eff = c(1,0,0,1,0,1,1,0), season_period = 4)
# draw diagnostic plots
drawDiagPlots("all", output_lm)
drawDiagPlots(output_lm$diagnostic_plots_names$FittedTimeSeries, output_lme)
# show summary of model fit
str(output_lm$fit_summary)
