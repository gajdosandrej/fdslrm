#' @title Estimating parameters of a linear mixed model (LMM) with a simple variance components structure
#'
#' @description
#' \code{mixed(y, X, Z, dim, s20, method, lambda, adaptRW)} computes ML, REML, MINQE(I), MINQE(U,I), BLUE(b), BLUP(u)
#' by Henderson's Mixed Model Equations Algorithm. This is R version of the original Matlab program
#' created by Viktor Witkovsky (Witkovsky, 2000).
#'
#' Model: \eqn{Y=X*b+Z*u+e,}
#'
#' \eqn{b=(b_1',...,b_f')'} and \eqn{u=(u_1',...,u_r')'},\eqn{ E(u)=0, Var(u)=diag(sigma^2_i*I_{m_i}), i=1,...,r
#' E(e)=0, Var(e)=sigma^2_{r+1}*I_n, Var(y)=Sig=sum_{i=1}^{r+1} sigma^2_i*Sig_i}.
#' We assume normality and independence of \eqn{u} and \eqn{e}.
#'
#' @param y n-dimensional vector of observations.
#' @param X (n * k)-design matrix for fixed effects b=(b_1,...,b_f), typically X=[X_1,...,X_f] for some X_i.
#' @param Z (n * m)-design matrix for random efects u=(u_1,...,u_r), typically Z=[Z_1,...,Z_r] for some Z_i.
#' @param dim vector of dimensions of u_i, i=1,...,r, dim=(m_1,...,m_r), m=sum(dim).
#' @param s20 a prior choice of the variance components, s20=(s20_1,...,s20_r,s20_{r+1}). SHOULD BE POSITIVE for \code{method>0}.
#' @param method method of estimation of variance components; 0:NO estimation, 1:ML, 2:REML, 3:MINQE(I), 4:MINQE(U,I)
#' @param lambda regularization parameter used for ridge regression weights,
#' default value is \code{lambda = numeric()} (use standatd estimation procedure, i.e.
#' no regularized ridge estimation procedure).
#' @param adaptRW flag for using adaptive method for the ridge matrix weights.
#' \code{adaptRW = FALSE}, the used ridge matrix is \code{Rw = lambda * diag(rep(1,k))}.
#' If \code{adaptRW = TRUE}, \code{Rw = lambda * diag(weights)}, where \code{weights = rep(1,k)/abs(b)}
#' and b is fixed effect estimate in current iteration. Default value is \code{adaptRW = FALSE}.
#'
#' @references
#' Witkovsky, V.: MATLAB Algorithm for solving Henderson's Mixed Model Equations. Technical Report No. 3/2000,
#' Institute of Measurement Science, Slovak Academy of Sciences, Bratislava, 2000.
#' https://www.mathworks.com/matlabcentral/fileexchange/200-mixed?s_tid=prof_contriblnk
#'
#' Searle, S.R., Cassela, G., McCulloch, C.E.: Variance Components.
#' John Wiley & Sons, INC., New York, 1992. (pp. 275-286).
#'
#' @return List with the following components:
#' \itemize{
#' \item \code{s2} - estimated vector of variance components \eqn{(sigma^2_1,..., sigma^2_{r+1})}'.
#' A warning message appears if some of the estimated variance components is negative or equal to zero.
#' In such cases the calculated Fisher information matrices are inadequate.
#' \item \code{b} - k-dimensional vector of estimated fixed effects beta,
#' \eqn{b=(b_1,...,b_f)=(X'Sig^{-1}X)^{+}X'Sig^{-1}y}.
#' \item \code{u} - m-dimensional vector of EBLUP's of random effects U, \eqn{u=(u_1,...,u_r)}.
#' \item \code{Is2} - Fisher information matrix for variance components; if \code{method=0} the output is empty matrix Is2;
#' if \code{metod=3} or \code{method=4} the output is inversion of the covariance matrix of MINQE calculated at estimated s2.
#' \item \code{C} - g-inverse of Henderson's MME matrix, where \eqn{C=MPinv([XX XZ; XZ' ZZ+inv(D)*s0]/s0)}, if inv(D) exists
#' or C=s0*[I 0; 0 D]*pinv([XX XZ*D; XZ' V]) otherwise
#' \item \code{H} - Criterial matrix for MINQE calculated at priors s20;
#' if method=3: \eqn{H_ij=trace(Sig_0^{-1}*Sig_i*Sig_0^{-1}*Sig_j)},
#' if method=4: \eqn{H_ij=trace((M*Sig_0*M)^{+}*Sig_i*(M*Sig_0*M)^{+}*Sig_j)}
#' \item \code{q} - (r+1)-dimensional vector of MINQE(U,I) quadratic forms calculated at prior values s20;
#' if \code{method=0,1,2} the output is empty vector q, otherwise \eqn{q_i=y'*(M*Sig_0*M)^{+}*Sig_i*(M*Sig_0*M)^{+}*y}.
#' \item \code{loglik} - Log-likelihood evaluated at the estimated parameters;
#' if \code{method=1}: loglik=log-likelihood(ML),
#' if \code{method=2}: loglik=log-likelihood(REML),
#' if \code{method=3}: or \code{method=4} loglik=numeric(),
#' if \code{method=0}: loglik=informative value of log of the joint pdf of (y,u).
#' \item \code{loops} - Number of loops.
#' }
#'
#' @note Ver.: 23-Apr-2020 19:44:40.
#'
#' @example R/Examples/example_mixed.R
#'
#' @export
#'
mixed <- function(y, X, Z, dim, s20, method, lambda, adaptRW) {

        # Adaptive ridge-estimation weights
        if(missing(adaptRW)) {
                adaptRW <- logical()
        }

        if(length(adaptRW) == 0) {
                adaptRW <- FALSE
        }

        # Parameter lambda for computing the regularized (ridge-estimation) ML/REML
        if(missing(lambda)) {
                lambda <- numeric()
        }

        # The required input parameters
        yy <- t(y)%*%y
        Xy <- t(X)%*%y
        Zy <- t(Z)%*%y
        XX <- t(X)%*%X
        XZ <- t(X)%*%Z
        ZZ <- t(Z)%*%Z
        a <- rbind(Xy,Zy)
        # Other input parameters
        n <- length(y)
        k <- ncol(X)
        m <- ncol(Z)
        rx <- Matrix::rankMatrix(XX)
        r <- length(s20)-1
        Im <- diag(rep(1,m))
        loops <- 0

        ## Method 0 (No estimation of variance components)
        # ======================================================================
        # METHOD=0:
        # No estimation of variance components
        # Output is BLUE(b), BLUP(u), and C,
        # calculated at chosen values s20
        # ======================================================================
        if(method == "0") {
                s0 <- s20[r+1]
                d <- rep(1,dim[1])
                id0 <- 0
                if(r >= 1) {
                        for(i in 1:r) {
                                id <- 1:dim[i]
                                d[id0 + id] <- s20[i] * d[id0 + id]
                                id0 <- id0 + dim[i]
                        }
                }
                D <- diag(d)
                V <- s0 * Im + ZZ %*% D
                A <- rbind(cbind(XX, XZ %*% D), cbind(t(XZ), V))
                A <- gnm::MPinv(A) # A = A \ Ikm
                if(k == 1) {
                        C <- s0*rbind(cbind(A[1:k,1:k],t(A[1:k,(k+1):(k+m)])),cbind(D%*%A[(k+1):(k+m),1:k],D%*%A[(k+1):(k+m),(k+1):(k+m)]))
                } else {
                        C <- s0*rbind(cbind(A[1:k,1:k],A[1:k,(k+1):(k+m)]),cbind(D%*%A[(k+1):(k+m),1:k],D%*%A[(k+1):(k+m),(k+1):(k+m)]))
                }
                bb <- A %*% a
                b <- bb[1:k]
                v <- bb[(k+1):(k+m)]
                u <- D %*% v
                Aux <- yy-t(b) %*% Xy - t(u) %*% Zy
                if(all(s20 != 0)) {
                        loglik <- (-1)*((n + m) * log(2 * pi) + n * log(s0) + log(prod(d)) + Aux / s0) / 2
                } else {
                        loglik <- numeric()
                }

                s2 <- s20
                Is2 <- matrix()
                H <- matrix()
                q <- vector()
                return(list("s2" = s2, "b" = b, "u" = u, "Is2" = Is2, "C" = C, "H" = H, "q" = q, "loglik" = loglik, "loops" = loops))
        }

        # ======================================================================
        # METHOD=1,2,3,4: ESTIMATION OF VARIANCE COMPONENTS
        # ======================================================================
        fk <- which(s20 <= 0)
        if(length(fk) > 0) {
                s20[fk] <- 100 * .Machine$double.eps * rep(1, length(fk))
                warning("Priors in s20 are negative or zeros !CHANGED!")
        }
        s21 <- s20
        ZMZ <- ZZ - t(XZ) %*% gnm::MPinv(XX) %*% XZ # ZMZ = ZZ - XZ' * (XX \ XZ)
        q <- rep(0, r+1)
        H <- matrix(nrow = r+1,ncol = r+1)
        # ======================================================================
        # START OF THE MAIN LOOP
        # ======================================================================
        epss <- 1e-8
        crit <- 1
        weight <- rep(1, k)
        small <- 0.1
        # H <- matrix()
        # q <- vector()

        while (crit > epss) {
                loops <- loops + 1
                sigaux <- s20
                s0 <- s20[r+1]
                d <- rep(1, sum(dim))
                id0 <- 0
                if(r >= 1) {
                        for(i in 1:r) {
                                id <- 1:dim[i]
                                d[id0 + id] <- s20[i] * d[id0 + id]
                                id0 <- id0 + dim[i]
                        }
                }

                D <- diag(d)
                V <- s0 * Im + ZZ %*% D
                # W <- s0 * solve(V)
                W <- s0 * pracma::mldivide(V, Im)
                # T <- solve(Im + ZMZ %*% D / s0)
                T <- pracma::mldivide((Im + ZMZ %*% D / s0), Im)


                if(length(lambda) == 0) {
                        A <- rbind(cbind(XX, XZ %*% D),cbind(t(XZ), V))
                        bb <- gnm::MPinv(A) %*% a # bb <- pracma::mldivide(A, a)
                } else {
                        # The regularized version
                        A <- rbind(cbind(XX + diag(weight)*lambda, XZ %*% D),cbind(t(XZ), V))
                        bb <- gnm::MPinv(A) %*% a # bb <- pracma::mldivide(A, a)
                }
                b <- bb[1:k]

                if(length(lambda) > 0 && adaptRW == TRUE) {
                        small0 <- small / 1.02
                        small <- small / 1.3
                        weight[abs(b) <= small0] <- 1 / small
                        weight[abs(b) > small0] <- 1 / abs(b[(abs(b) > small0)])
                }
                v <- bb[(k+1):(k+m)]
                u <- D %*% v

                # ======================================================================
                # ESTIMATION OF ML AND REML OF VARIANCE COMPONENTS
                # ======================================================================
                iupp <- 0
                q <- rep(1, r+1)
                for(i in 1:r){
                        ilow <- iupp+1
                        iupp <- iupp+dim[i]
                        Wii <- W[ilow:iupp,ilow:iupp]
                        Tii <- T[ilow:iupp,ilow:iupp]
                        w <- u[ilow:iupp]
                        ww <- t(w)%*%w
                        q[i] <- ww/(s20[i]*s20[i])
                        s20[i] <- ww/(dim[i]-matrixcalc::matrix.trace(as.matrix(Wii)))
                        if(!is.finite(s20[i])) {
                                s20[i] <- .Machine$double.eps
                        }
                        s21[i] <- ww/(dim[i]-matrixcalc::matrix.trace(as.matrix(Tii)))
                        if(!is.finite(s21[i])) {
                                s21[i] <- .Machine$double.eps
                        }
                }

                Aux <- yy-t(b)%*%Xy-t(u)%*%Zy
                Aux1 <- Aux-(t(u)%*%v)*s20[r+1]
                q[r+1] <- Aux1/(s20[r+1]*s20[r+1])
                s20[r+1] <- Aux/n
                s21[r+1] <- Aux/(n-rx)

                if(method == "1") {
                        crit <- sqrt(sum((sigaux-s20)^2))
                        H <- matrix(nrow = r+1,ncol = r+1)
                } else if(method == "2") {
                        s20 <- s21
                        crit <- sqrt(sum((sigaux-s20)^2))
                        H <- matrix(nrow = r+1,ncol = r+1)
                } else {
                        crit <- 0
                }
                # stop()
                # if(loops==max_iter) {
                #         warning("Maximum number of iterations reached.")
                #         break
                # }


        }
        # ======================================================================
        # END OF THE MAIN LOOP
        # ======================================================================
        # COMPUTING OF THE MINQE CRITERIAL MATRIX H
        # ======================================================================
        if (method=="3" || method=="4") {
                H <- diag(r+1)
                if(method=="4") {
                        W <- T
                        H[r+1,r+1] <- (n-rx-m+matrixcalc::matrix.trace(W%*%W))/(sigaux[r+1]*sigaux[r+1]) #%VW
                } else {
                        H[r+1,r+1] <- (n-m+matrixcalc::matrix.trace(W%*%W))/(sigaux[r+1]*sigaux[r+1])
                }

                iupp <- 0
                for(i in 1:r) {
                        ilow <- iupp+1
                        iupp <- iupp+dim[i]
                        trii <- matrixcalc::matrix.trace(as.matrix(W[ilow:iupp,ilow:iupp]))
                        trsum <- 0
                        jupp <- 0
                        for(j in 1:r) {
                                jlow <- jupp+1
                                jupp <- jupp+dim[j]
                                tr <- matrixcalc::matrix.trace(as.matrix(W[ilow:iupp,jlow:jupp]%*%W[jlow:jupp,ilow:iupp]))
                                trsum <- trsum+tr
                                H[i,j] <- (1*(i==j)*(dim[i]-2*trii)+tr)/(sigaux[i]*sigaux[j])
                        }
                        H[r+1,i] <- (trii-trsum)/(sigaux[r+1]*sigaux[i])
                        H[i,r+1] <- H[r+1,i]
                }
        }
        s2 <- numeric()

        # ======================================================================
        # MINQE(I), MINQE(U,I), ML, AND REML
        # ======================================================================
        if(method=="3" || method=="4") {
                s2 <- gnm::MPinv(H)%*%q
                loglik <- numeric()
        } else {
                s2<- s20
        }

        fk <- which(s2 < 10*epss)
        if(length(fk) > 0) {
                warning("Estimated variance components are negative or zeros!")
        }
        # ======================================================================
        # BLUE, BLUP, THE MME'S C MATRIX AND THE LOG-LIKELIHOOD
        # ======================================================================
        s0 <- s2[r+1]
        d <- rep(1, sum(dim))
        id0 <- 0
        if(r >= 1) {
                for(i in 1:r) {
                        id <- 1:dim[i]
                        d[id0+id] <- s20[i] * d[id0 + id]
                        id0 <- id0 + dim[i]
                }
        }

        D <- diag(d)
        V <- s0*Im+ZZ%*%D
        # W<-s0*solve(V)
        W <- s0 * pracma::mldivide(V, Im)
        # T<-solve(Im+ZMZ%*%D/s0)
        T <- pracma::mldivide(Im + ZMZ %*% D / s0, Im)
        if(length(lambda) == 0) {
                A <- rbind(cbind(XX,XZ%*%D),cbind(t(XZ),V))
                A <- gnm::MPinv(A)
        } else {
                # Regularized version
                A <- rbind(cbind(XX + diag(weight) * lambda, XZ%*%D),cbind(t(XZ),V))
                A <- gnm::MPinv(A)
        }

        if(k==1) {
                C <- s0*rbind(cbind(A[1:k,1:k],t(A[1:k,(k+1):(k+m)])),cbind(D%*%A[(k+1):(k+m),1:k],D%*%A[(k+1):(k+m),(k+1):(k+m)]))
        } else {
                C <- s0*rbind(cbind(A[1:k,1:k],A[1:k,(k+1):(k+m)]),cbind(D%*%A[(k+1):(k+m),1:k],D%*%A[(k+1):(k+m),(k+1):(k+m)]))
        }

        bb <- A%*%a
        b <- bb[1:k]
        v <- bb[(k+1):(k+m)]
        u <- D%*%v
        if(method == "1") {
                loglik <- -(n*log(2*pi*s0)-log(det(W))+n)/2
        } else if(method == "2") {
                # loglik=-((n-rx)*log(2*pi*s0)-log(det(T))+(n-rx))/2
                # Here is the adjusted loglik version - like in SAS
                loglik <- -((n-rx) * log(2*pi*s0) - log(det(T)) + (n-rx))/2 - log(det(t(X)%*%X))/2
                loglik <- -2 * loglik
        }
        # ======================================================================
        # FISHER INFORMATION MATRIX FOR VARIANCE COMPONENTS
        # ======================================================================
        Is2 <- diag(r+1)
        if (method == "2" || method == "4") {
                W<-T
                Is2[r+1,r+1] <- (n-rx-m+matrixcalc::matrix.trace(W%*%W))/(s2[r+1]*s2[r+1]) #%VW
        } else {
                Is2[r+1,r+1] <- (n-m+matrixcalc::matrix.trace(W%*%W))/(s2[r+1]*s2[r+1])
        }

        iupp <- 0
        for(i in 1:r) {
                ilow <- iupp+1
                iupp <- iupp+dim[i]
                trii <- matrixcalc::matrix.trace(as.matrix(W[ilow:iupp,ilow:iupp]))
                trsum <- 0
                jupp <- 0
                for(j in 1:r) {
                        jlow <- jupp+1
                        jupp <- jupp+dim[j]
                        tr <- matrixcalc::matrix.trace(as.matrix(W[ilow:iupp,jlow:jupp]%*%W[jlow:jupp,ilow:iupp]))
                        trsum <- trsum+tr
                        Is2[i,j] <- (1*(i == j)*(dim[i]-2*trii)+tr)/(s2[i]*s2[j])
                }
                Is2[r+1,i] <- (trii-trsum)/(s2[r+1]*s2[i])
                Is2[i,r+1] <- Is2[r+1,i]
        }
        Is2 <- Is2/2
        # ======================================================================
        # EOF MIXED.M
        # ======================================================================
        result <- list("s2" = as.vector(s2), "b" = b, "u" = u, "Is2" = Is2, "C" = C, "H" = H, "q" = q, "loglik" = loglik, "loops" = loops)
        return(result)
}



Design2 <- function(N) {
        # DESIGN2 Creates the design matrices for two-way classification
        # model y_ijk = a_i + b_j + c_ij + e_ijk,
        # given by its incidence matrix N (ixj).
        # The incidence matrix N could have some cells equal to zero.

        sizeN <- vector()
        if(is.vector(N)) {
                sizeN <- c(1, length(N))
        } else if (is.matrix(N)) {
                sizeN <- dim(N)
        }

        M <- t(Conj(N))
        v <- c(M)
        v <- v[v > 0]
        oj <- matrix(1, sizeN[2], 1)
        n <- N %*% oj

        A <- Diagm(n)
        B <- matrix()
        if(is.vector(N)) {
                B <- Diagm(N)
        } else if (is.matrix(N)) {
                B <- Diagm(N[1,])
        }

        if(is.matrix(N) && sizeN[1] > 1) {
                for(i in 2:sizeN[1]) {
                        B <- rbind(B, Diagm(N[i,]))
                }
        }

        C <- Diagm(v)

        return(list("A" = A, "B" = B, "C" = C))

}



Diagm <- function(n) {
        # DIAGM Creates block-diagonal matrix from 1-vectors
        # with dimensions given in n=c(n_1,...,n_k).

        ones_list <- list()
        for (i in 1:length(n)) {
                ones_list[[i]] <- rep(1,n[i])
        }

        D <- as.matrix(Matrix::bdiag(ones_list))

        return(D)
}
