## EXAMPLE 1
times <- c(1:24)
freq_fixed <- c(1/24)
freq_rand <- c(1/12, 1/8)
matF <- makeF(times, freq_fixed)
matV <- makeV(times, freq_rand)
X <- c(40.3, 40.7, 38.5, 37.9, 38.6, 41.1, 45.2, 45.7, 46.7, 46.5, 45.2, 45.1, 45.8, 46.3, 47.5, 48.5, 49.1, 51.7, 50.6, 48, 44.7, 41.2, 40, 40.3)
doolse <- DOOLSE_orth(X, matF, matV)
print(doolse)
