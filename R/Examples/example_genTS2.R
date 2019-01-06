## EXAMPLE 1
# realization from general FDSLRM
model_formula <- formula(x ~ cos(2 * pi / 24 * t) + sin(2 * pi / 24 * t) +
                                 cos(2 * pi / 8 * t) + sin(2 * pi / 8 * t) +
                                 cos(2 * pi / 6 * t) + sin(2 * pi / 6 * t))
beta <- c (-3, 8, 7)
nu <- c(3, 2.7, 1.5, 1, 2)
times <- 1:24
n <- length(times)
realization <- genTS2(model_formula, beta, nu, n)
print(realization)
plot(times, realization, xlab = "time", type = "o")

## EXAMPLE 2
# realization from classical linear regression time series model (special case of FDSLRM)
model_formula <- formula(x ~ cos(2 * pi / 24 * t) + sin(2 * pi / 24 * t))
beta <- c (-3, 8, 7)
nu <- 3
times <- 1:24
n <- length(times)
realization <- genTS2(model_formula, beta, nu, n)
print(realization)
plot(times, realization, xlab = "time", type = "o")
