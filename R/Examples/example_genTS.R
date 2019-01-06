## EXAMPLE 1
# realization from general FDSLRM
beta <- c (-3, 8, 7)
times <- 1:24
matF <- makeF(times, c(1/24))
matV <- makeV(times, c(1/8, 1/6))
meanv <- matF %*% beta
variances <- c(3, 2.7, 1.5, 1, 2)
realization <- genTS(meanv, matV, variances)
print(realization)
plot(times, realization, xlab = "time", type = "o")

## EXAMPLE 2
# realization from classical linear regression time series model (special case of FDSLRM)
beta <- c (-3, 8, 7)
times <- 1:24
matF <- makeF(times, c(1/24))
meanv <- matF %*% beta
variances <- 3
realization <- genTS(meanv, var = variances)
print(realization)
plot(times, realization, xlab = "time", type = "o")
