## EXAMPLE 1
times <- c(1:24)
frequencies <- c(1/24, 1/8, 1/6)
matF <- makeF(times, frequencies)
matMF <- makeM_F(matF)
print(matMF)
