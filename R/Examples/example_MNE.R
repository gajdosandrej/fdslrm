## EXAMPLE 1
times <- c(1:40)
freq_fixed <- c(1/24)
freq_rand <- c(1/8, 1/6)
matF <- makeF(times, freq_fixed)
matV <- makeV(times, freq_rand)
X <- rnorm(40, 2, 1.25)
mne <- MNE(X, matF, matV)
print(mne)
