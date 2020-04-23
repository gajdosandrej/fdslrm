## EXAMPLE 1:
IncN <- t(matrix(c(4, 5, 8, 9, 5, 10, 15, 20), 4, 2))
matrices <- Design2(IncN)
n <- nrow(matrices$A)
n1 <- ncol(matrices$B)
n2 <- ncol(matrices$C)
X <- cbind(matrix(1, n, 1), matrices$A)
Z <- cbind(matrices$B, matrices$C)
btrue <- Conj(c(1, 2, 3))
s2true <- c(0.5, 3, 1)
u1 <- sqrt(s2true[1]) * rnorm(n1)
u2 <- sqrt(s2true[2]) * rnorm(n2)
u <- c(u1, u2)
e <- sqrt(s2true[3]) * rnorm(n)
y <- as.vector(X %*% btrue + Z %*% u + e)
dim <- c(n1, n2)
s20 <- c(1, 1, 1)
method <- 2     # 0:NONE, 1:ML, 2:REML, 3:MINQE(I), 4:MINQE(U,I)
result1 <- mixed(y, X, Z, dim, s20, method)
result1$s2
result1$b
result1$u
result1$Is2
result1$C
result1$H
result1$q
result1$loglik
result1$loops
