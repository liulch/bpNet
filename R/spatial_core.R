## tscs spatial: core function

## Y: outcome a N*T matrix, each column represents observations for all units at a given period 
## W: W matrix: a N*N*T cube
## X: covariates with constant effect, a N*p1*T cube
## Z: covariates with unit-level random effect, a T*p2*N cube
## A: covariates with time-level random effect, a N*p3*T cube

## niter: number of mcmc iterations

## constantRho: a logical flag indicating whether to incorporate a systematic rho
## randomRho: a logical flag indicating whether to incorporate time-varying rho_ts

## spec: specifications for rho_t: multilevel, random walk or ar(1)

bspanel.core <- function(Y,  
	                     W,
		                 lagY = NULL,
		                 X = NULL,
		                 Z = NULL,
		                 A = NULL,
		                 rhoZ = NULL, ## covariates used to interprete rho_t
		                 r = 0,
		                 flasso = 1, 
		                 burn = 1000,
		                 niter = 10000,
		                 constantRho = TRUE,
		                 randomRho = FALSE, 
		                 spec = "rw" ## rw, multilevel, ar1
		                 ){

	if (sum(constantRho + randomRho) > 1) {
		stop("Either fixed or random.\n")
	}

	TT <- dim(Y)[2]
	N <- dim(Y)[1]

	WY <- matrix(NA, N, TT)
	for (i in 1:TT) {
		WY[, i] <- W[,,i] %*% as.matrix(Y[, i])
	}

	p1 <- ifelse(!is.null(X), dim(X)[2], 0)
	p2 <- ifelse(!is.null(Z), dim(Z)[2], 0)
	p3 <- ifelse(!is.null(A), dim(A)[2], 0)
	p4 <- ifelse(!is.null(rhoZ), dim(rhoZ)[2], 0)
	p5 <- ifelse(!is.null(lagY), 1, 0)

	lagy <- NULL
	if (p5 == 1) {
		lagy <- array(NA, dim = c(N, 1, TT))
		for (i in 1:TT) {
			lagy[,,i] <- as.matrix(lagY[, i])
		}
	}

	res <- matrix(0, N, TT)

	xfit <- matrix(0, N, TT)
	zfit <- matrix(0, N, TT)
	afit <- matrix(0, N, TT)
	fefit <- matrix(0, N, TT)
	yfit <- matrix(0, N, TT)

	r0fit <- matrix(0, N, TT)  ## constant rho fit 
	rfit <- matrix(0, N, TT)   ## random rho fit

	Lambda <- matrix(NA, N, TT)
	for (i in 1:TT) {
		Lambda[, i] <- eigen(W[,,i])[[1]]
	}

	## spec <- "ar1"  ## ar1, rw, multilevel
	sigma_n2 <- sigma_n20 <- 1

	## type : 1: rw, 2: multi-level, 3: ar1
	type <- NULL
	
	kappa <- kappa0 <- 0
	if (spec == "rw") {
		kappa <- 1
		type <- 1
	} else if (spec == "multilevel") {
		type <- 2
	} else {
		type <- 3
	}

	## initialize 
	b0.all <- b0 <- B0 <- B0.all <- NULL
	if (p1 > 0) {
		beta0 <- olsBeta(Y, X)
		beta <- beta0
		xfit <- Xfit(X, beta)

		b0 <- as.matrix(rep(0, p1))
	    B0 <- matrix(0, p1, p1)

		for (i in 1:p1) {
			B0[i, i] <- 10000
		}
	}

	k0 <- K0 <- NULL
	if (p5 > 0) {
		k0 <- as.matrix(0)
		K0 <- matrix(1000, 1, 1)
	}

	ra0 <- RA0 <- NULL
	if (p4 > 0) {
		ra0 <- as.matrix(rep(0, p4))
		RA0 <- matrix(0, p4, p4)

		for (i in 1:p4) {
			RA0[i, i] <- 10000
		}
	}

	## 3. Xi
	#Xi <- NULL
	#VX <- c0 <- C0 <- NULL
	#if (p3 > 0) {

	#	Xi <- Xi0 <- matrix(0, p3, TT)
	#	afit <- Afit(A, Xi)
	#	VX <- diag(10000, p3, p3)

	#	c0 <- 0.001
	#	C0 <- diag(0.001, p3, p3)

	#}

	## 4. multi-level alpha
	#Alpha <- NULL
	#VA <- a0 <- A0 <- NULL
	#if (p2 > 0) {
		
	#	Alpha <- Alpha0 <- matrix(0, p2, N)
	#	zfit <- Zfit(Z, Alpha)
	#	VA <- diag(10000, p2, p2)

	#	a0 <- 0.001
	#	A0 <- diag(0.001, p2, p2)

	#}

	## 5. Factor
	#L <- NULL
	#ng <- max(c(index))
	#if (r > 0) {
	#	L <- matrix(rnorm(N * r), N, r)
	#}

	#Factor <- NULL
	#if (r > 0) {
	#	Factor <- matrix(rnorm(TT * r), TT, r)
	#}

	omega <- NULL
	lg2 <- 10
	if (flasso == 1 && r > 0) {
		omega <- matrix(1, r, 1)

		b0.all <- as.matrix(rep(0, p1 + r))
	    B0.all <- matrix(0, p1 + r, p1 + r)

		for (i in 1:(p1 + r)) {
			B0.all[i, i] <- 10000
		}

	}

	## factor structure
	#if (r > 0) {
	#	fefit <- L %*% t(Factor)
	#}

	L <- VA <- a0 <- A0 <- NULL
	if ((p2 + r) > 0) {
		
		if (p2 > 0) {
			Alpha <- Alpha0 <- matrix(0, p2, N)
			zfit <- Zfit(Z, Alpha)
		}
		
		L <- matrix(0, N, r) ## initialize loadings matrix
		if (r > 0) {
			for (i in 1:r) {
				L[, i] <- c(rep(0, i-1), rnorm(N - i + 1))
			}
		}
		
		VA <- matrix(0, p2 + r, p2 + r)


		a0 <- 0.001
		A0 <- NULL
		
        if (flasso == 0 || r == 0) {
        	A0 <- matrix(0, p2 + r, p2 + r)

			for (i in 1:(p2 + r)) {
				A0[i, i] <- 0.001
				VA[i, i] <- 1000
			}
        } else {

        	if (p2 > 0) {
        		A0 <- matrix(0, p2, p2)
        		for (i in 1:p2) {
	        		A0[i, i] <- 0.001
					VA[i, i] <- 1000
	        	}
        	}
        	
        	for (i in (p2+1):(p2 + r)) {
        		VA[i, i] <- 1
        	}
        }
		
	}

	Factor <- VX <- c0 <- C0 <- NULL
	if ((p3 + r) > 0) {
	    if (p3 > 0) {
	    	Xi <- Xi0 <- matrix(0, p3, TT)
			afit <- Afit(A, Xi)
	    }

	    Factor <- matrix(0, TT, r) ## initialize loadings matrix
		if (r > 0) {
			Factor <- matrix(rnorm(TT*r), TT, r)
			
		}

		VX <- matrix(0, p3 + r, p3 + r)

		c0 <- 0.001
		C0 <- NULL 

		if (flasso == 0 || r == 0) {
			C0 <- matrix(0, p3 + r, p3 + r)

			for (i in 1:(p3 + r)) {
				C0[i, i] <- 0.001
				VX[i, i] <- 1000
			}

		} else {

			if (p3 > 0) {
				C0 <- matrix(0, p3, p3)
				for (i in 1:p3) {
					C0[i, i] <- 0.001
					VX[i, i] <- 1000
				}
			}

			for (i in (p3+1):(p3 + r)) {
        		VX[i, i] <- 1
        	}

		}
		
	}

	#O0 <- NULL
	#if (r > 0) {
	#	O0 <- diag(10000, r, r)
	#}


	#L <- VA <- a0 <- A0 <- NULL
	#if ((p2 + r) > 0) {
		
	#	if (p2 > 0) {
	#		Alpha <- Alpha0 <- matrix(0, p2, N)
	#		zfit <- Zfit(Z, Alpha)
	#	}
		
	#	L <- matrix(0, N, r) ## initialize loadings matrix
	#	if (r > 0) {
	#		for (i in 1:r) {
	#			L[, i] <- c(rep(0, i-1), 1, rnorm(N - i))
	#		}
	#	}
		
	#	VA <- matrix(0, p2 + r, p2 + r)

	#	a0 <- 0.001
	#	A0 <- matrix(0, p2 + r, p2 + r)

	#	for (i in 1:(p2 + r)) {
	#		A0[i, i] <- 0.001
	#		VA[i, i] <- 1000
	#	}
	#}

	#Factor <- VX <- c0 <- C0 <- NULL
	#if ((p3 + r) > 0) {
	#    if (p3 > 0) {
	#    	Xi <- Xi0 <- matrix(0, p3, TT)
	#		afit <- Afit(A, Xi)
	#    }

	#   Factor <- matrix(0, TT, r) ## initialize loadings matrix
	#	if (r > 0) {
	#		Factor <- matrix(rnorm(TT*r), TT, r)
			
	#	}

	#	VX <- matrix(0, p3 + r, p3 + r)

	#	c0 <- 0.001
	#	C0 <- matrix(0, p3 + r, p3 + r)

	#	for (i in 1:(p3 + r)) {
	#		C0[i, i] <- 0.001
	#		VX[i, i] <- 1000
	#	}
	#}

	## factor structure
	#if (r > 0) {
	#	fefit <- L %*% t(Factor)
	#}

	rho0_old <- rho0 <- 0
	#if (constantRho == 1 && rho0 != 0) {
	#	r0fit <- rhofit(WY, matrix(rho0, TT, 1))
	#}

	Rho_old <- Rho <- Rho0 <- matrix(0, TT, 1)
	#if (randomRho == 1 && sum(abs(c(Rho0))) != 0) {
	#	rfit <- rhofit(WY, Rho)
	#}

	## covariates for rho_t
	rhoA <- matrix(0, p4, 1)

	## initialize error term variance
	res <- Y - xfit - zfit - afit - r0fit - rfit - fefit - yfit
	sigma2 <- sampleSigma2(res)
	## sigma2 <- 1 


	## store results
	beta_i <- NULL
	if (p1 > 0) {
		beta_i <- matrix(NA, p1, niter - burn)
	}

	rhoy_i <- NULL
	if (p5 > 0) {
		rhoy_i <- matrix(NA, 1, niter - burn)
	}

	rhoA_i <- NULL
	if (p4 > 0) {
		rhoA_i <- matrix(NA, p4, niter - burn)
	}

	Alpha_i <- NULL
	if (p2 > 0) {
		Alpha_i <- array(NA, dim = c(p2, N, niter - burn))
	}

	Xi_i <- NULL
	if (p3 > 0) {
		Xi_i <- array(NA, dim = c(p3, TT, niter - burn))
	}

	L_i <- NULL
	if (r > 0) {
		L_i <- array(NA, dim = c(N, r, niter - burn))
	}

	Factor_i <- NULL
	if (r > 0) {
		Factor_i <- array(NA, dim = c(TT, r, niter - burn))
	}

	omega_i <- NULL
	#o_i <- NULL
	if (flasso == 1 && r > 0) {
		omega_i <- matrix(NA, r, niter - burn)
		#o_i <- matrix(NA, r, niter - burn)
	} 

	rho0_i <- NULL
	if (constantRho == TRUE) {
		rho0_i <- matrix(NA, 1, niter - burn)
	}

	Rho_i <- NULL
	kappa_i <- NULL
	sigma_n2_i <- NULL
	if (randomRho == TRUE) {
		Rho_i <- matrix(NA, TT, niter - burn)
		sigma_n2_i <- matrix(NA, 1, niter - burn)
		if (spec == "ar1") {
			kappa_i <- matrix(NA, 1, niter - burn)
		}
	}


	sigma2_i <- matrix(NA, 1, niter - burn)
	## yfit_i <- array(NA, dim = c(N, TT, niter - burn))


	## iterations start!!!
	for (i in 1:niter) {

		## ----------- 1. update beta ---------- ##	
	    if (flasso == 0 || r == 0) {
	    	## no shrinkage on factor loadings
	    	if (p1 > 0) {
		    	res <- Y - zfit - afit - r0fit - rfit - fefit - yfit
			    beta <- sampleBeta(X, res, B0, b0, sigma2)
			    xfit <- Xfit(X, beta)

			    if (i > burn) {
			    	beta_i[, i - burn] <- beta
			    }	    
		    }
	    } else {
	    	## jointly update beta and omega
	    	res <- Y - zfit - afit - r0fit - rfit - yfit 

	    	X.all <- array(NA, dim = c(N, p1 + r, TT))

	    	if (p1 > 0) {

	    		for (j in 1:TT) {
	    			X.all[,, j] <- cbind(as.matrix(X[,,j]), L * matrix(rep(Factor[j, ], each = N), N, r))
	    		}

	    		#B0.all <- matrix(0, p1 + r, p1 + r) 
	    		#b0.all <- matrix(0, p1 + r, 1)
	    		#b0.all[1:p1, ] <- c(b0)

	    		#B0.all[1:p1, 1:p1] <- B0
	    		#for (i in (p1 + 1):(p1 + r)) {
	    		#	B0.all[i, i] <- 
	    		#}

	    		beta.all <- sampleBeta(X.all, res, B0.all, b0.all, sigma2)

	    		beta <- as.matrix(beta.all[1:p1, ])
			    xfit <- Xfit(X, beta)

			    if (i > burn) {
			    	beta_i[, i - burn] <- beta
			    }

			    omega <- as.matrix(beta.all[(p1 + 1):(p1 + r), ])

			    for (j in 1:r) {
                    pa1 <- sqrt(lg2/omega[j,1]^2)
                    pa2 <- lg2
                    B0.all[p1 + j, p1 + j] <- 1 / rrinvgauss(pa1, pa2)  
                }
                ## sample lambda_gamma2 
                pa1 <- 0.001 + r
                pa2 <- 1/(0.001 + sum(diag(B0.all)[(p1+1):(p1+r)])/2)
                lg2 <- sampleG(pa1, pa2)

	    	} else {

	    		for (j in 1:TT) {
	    			X.all[,, j] <- L * matrix(rep(Factor[j, ], each = N), N, r)
	    		}

	    		#b0.all <- matrix(0, r, 1)
	    		#B0.all <- matrix(0, r, r)

	    		#for (i in (p1 + 1):(p1 + r)) {
	    		#	B0.all[i, i] <- 
	    		#}

	    		omega <- sampleBeta(X.all, res, B0.all, b0.all, sigma2)

			    for (j in 1:r) {
                    pa1 <- sqrt(lg2/omega[j,1]^2)
                    pa2 <- lg2
                    B0.all[j, j] <- 1 / rrinvgauss(pa1, pa2)  
                }
                ## sample lambda_gamma2 
                pa1 <- 0.001 + r
                pa2 <- 1/(0.001 + sum(diag(B0.all))/2)
                lg2 <- sampleG(pa1, pa2)

	    	}

	    	fefit <- (L * matrix(rep(c(omega), each = N), N, r)) %*% t(Factor)
	    }

	    ## ar1 coef
	    if (p5 == 1) {
	    	res <- Y - zfit - afit - r0fit - rfit - fefit - xfit
		    rhoy <- sampleBeta(lagy, res, K0, k0, sigma2)
		    yfit <- Xfit(lagy, rhoy)

		    if (i > burn) {
		    	rhoy_i[, i - burn] <- rhoy[1, 1]
		    }

	    }
		

		## ------------ 2. update unit-level random coefficients --------- ##
	    #if (p2 > 0) {
	    #	res <- Y - xfit - afit - r0fit - rfit - fefit - yfit
	 
			### 1. update coef
		#	Alpha <- sampleAlpha(Z, res, VA, 0, sigma2) ## (r + p2) * N 
			#if (p2 > 0) {
			#	if (p2 > 1) {
			#		Alpha <- as.matrix(Alpha_all[1:p2,])
			#	} else {
			#		Alpha <- matrix(Alpha_all[1,], 1, N)
			#	}
			#	zfit <- Zfit(Z, Alpha)
			#	Alpha_i[,,i] <- Alpha
			#}

		#	zfit <- Zfit(Z, Alpha)
		#	if (i > burn) {
		#		Alpha_i[,,i - burn] <- Alpha
		#	}
			

			### 2. update coef
		#	VA <- sampleD(A0, Alpha, a0)
			
	    #}

	    ## ------------ 2. update unit-level random coefficients --------- ##
	    if ((p2 + r) > 0) {
	    	res <- Y - xfit - afit - r0fit - rfit - yfit

	    	Z_all <- NULL
	    	if (r == 0) {
	    		Z_all <- Z
	    	} else {
	    		nFactor <- Factor
    			if (flasso == 1) {
    				nFactor <- nFactor * matrix(rep(c(omega), each = TT), TT, r)
    			} 
	    		if (p2 == 0) {
	    			Z_all <- array(NA, dim = c(TT, r, N))
	    			for (ii in 1:N) {
	    				Z_all[,, ii] <- nFactor
	    			}
    			} else {
    				Z_all <- array(NA, dim = c(TT, (r + p2), N)) 
    				for (ii in 1:N) {
    					Z_all[,, ii] <- cbind(Z[,,ii], nFactor)
    				}
    			}
	    	}
	 
			### 1. update coef
			Alpha_all <- sampleAlpha(Z_all, res, VA, 1, r, sigma2) ## (r + p2) * N 
			if (p2 > 0) {
				if (p2 > 1) {
					Alpha <- as.matrix(Alpha_all[1:p2,])
				} else {
					Alpha <- matrix(Alpha_all[1,], 1, N)
				}
				zfit <- Zfit(Z, Alpha)
				if (i > burn) {
					Alpha_i[,,i - burn] <- Alpha
				}	
			}

			if (r > 0) {
				if (r == 1) {
					L <- matrix(Alpha_all[(p2+1),], 1, N)
				} else {
					L <- as.matrix(Alpha_all[(p2+1):(p2+r),])
				}
				#if (r > burn) {
				#	L_i[,,i - burn] <- L
				#}
				L <- t(L)
			}

			### 2. update coef
			if (flasso == 0 || r ==0) {
				VA <- sampleD(A0, Alpha_all, a0)
			} else {
				if (p2 > 0) {
					VA.sub <- sampleD(A0, Alpha, a0)
				    VA[1:p2, 1:p2] <- VA.sub
				}				
			}			
	    }

		

		## ------------ 3. update time-level random coefficients ------------ ##
		#if (p3 > 0) {
		#	res <- Y - xfit - zfit - r0fit - rfit - fefit - yfit


	    	### 1. update coef
		#	Xi <- sampleAlpha(A, res, VX, 1, sigma2) ## (r + p2) * N 
			#if (p3 > 0) {
			#	if (p3 > 1) {
			#		Xi <- as.matrix(Xi_all[1:p3,])
			#	} else {
			#		Xi <- matrix(Xi_all[1,], 1, TT)
			#	}

			#	afit <- Afit(A, Xi)
			#	Xi_i[,,i] <- Xi
			#}

		#	afit <- Afit(A, Xi)
		#	if (i > burn) {
		#		Xi_i[,,i - burn] <- Xi
		#	}
			

			### 2. update coef
		#	VX <- sampleD(C0, Xi, c0)

		#}

		## ------------ 3. update time-level random coefficients ------------ ##
		if ((p3 + r) > 0) {
			res <- Y - xfit - zfit - r0fit - rfit - yfit

			A_all <- NULL
	    	if (r == 0) {
	    		A_all <- A
	    	} else {
	    		nL <- L
    			if (flasso == 1) {
    				nL <- nL * matrix(rep(c(omega), each = N), N, r)
    			} 
	    		if (p3 == 0) {
	    			A_all <- array(NA, dim = c(N, r, TT))
	    			for (ii in 1:TT) {
	    				A_all[,, ii] <- nL
	    			}
    			} else {
    				A_all <- array(NA, dim = c(N, (r + p3), TT)) 
    				for (ii in 1:TT) {
    					A_all[,, ii] <- cbind(A[,,ii], nL)
    				}
    			}
	    	}

	    	### 1. update coef
			Xi_all <- sampleAlpha(A_all, res, VX, 0, r, sigma2) ## (r + p2) * N 
			if (p3 > 0) {
				if (p3 > 1) {
					Xi <- as.matrix(Xi_all[1:p3,])
				} else {
					Xi <- matrix(Xi_all[1,], 1, TT)
				}

				afit <- Afit(A, Xi)
				if (i > burn) {
					Xi_i[,,i - burn] <- Xi
				}
			}

			if (r > 0) {
				if (r == 1) {
					Factor <- matrix(Xi_all[(p3+1),], 1, TT)
				} else {
					Factor <- as.matrix(Xi_all[(p3+1):(p3+r),])
				}
				#if (i > burn) {
				#	Factor_i[,,i - burn] <- Factor
				#}
				
				Factor <- t(Factor)
			}

			### 2. update coef
			if (flasso == 0 || r == 0) {
				VX <- sampleD(C0, Xi_all, c0)
			} else {
				if (p3 > 0) {
					VX.sub <- sampleD(C0, Xi, c0)
				    VX[1:p3, 1:p3] <- VX.sub
				}	
			}
			

		}

		## ----------- update factor structure fitting -------------- ##
		if (r > 0) {
			if (flasso == 0) {
				fefit <- L %*% t(Factor)
			} else {
				fefit <- nL %*% t(Factor)
                ## permutation 
				
				#permu <- permute(omega, Factor)
                #omega <- permu$omega
                #Factor <- permu$Xi 

                for (j in 1:r) {
                    con <- runif(1, 0, 1)
                    if (con > 1/4 && con <= 1/2) {
                        L[, j] <- -L[, j]
                        Factor[, j] <- -Factor[, j]
                    } else if (con > 1/2 && con <= 3/4) {
                        L[, j] <- -L[, j]
                        omega[j, 1] <- -omega[j, 1] 

                    } else if (con > 3/4) {
                         Factor[, j] <- -Factor[, j]
                         omega[j, 1] <- -omega[j, 1] 
                    }
                }


                if (i > burn) {
			    	omega_i[, i - burn] <- omega
			    }
			}

			if (i > burn) {
		    	L_i[,,i - burn] <- t(L)
		    	Factor_i[,,i - burn] <- t(Factor)
		    }
		}

		## ----------- update factor structure fitting -------------- ##
		#if (r > 0) {
		#	res <- Y - xfit - zfit - afit - r0fit - rfit - yfit

			## 1. update loadings
		#	L <- sampleL(r, sigma2, res, Factor, omega, index)

			## 2. update factor 
		#	Factor <- sampleF(r, sigma2, res, L, omega, index)

			## 3. update omega
		#	omega <- sampleOmega(res, Factor, L, index, O0, sigma2) 

			## 4. permutation
		#	 permu <- permute(omega, Factor)
        #    omega <- permu$omega
        #    Factor <- permu$Xi


        #    if (flasso == TRUE) {
        #        for (j in 1:r) {
        #            pa1 <- sqrt(lg2/omega[j,1]^2)
        #            pa2 <- lg2
        #            O0[j, j] <- 1 / rrinvgauss(pa1, pa2)  
        #        }
                ## sample lambda_gamma2 
        #        pa1 <- 0.001 + r
        #        pa2 <- 1/(0.001 + sum(diag(O0))/2)
        #        lg2 <- sampleG(pa1, pa2)   
        #    } 

			## 5. update omega variance

		#	fefit <- Fefit(Factor, L, omega, index)
			
		#	if (i > burn) {
		#		Factor_i[,,i - burn] <- Factor
		#		L_i[,,i - burn] <- L
		#		omega_i[, i - burn] <- omega
		#		o_i[, i - burn] <- diag(O0)
		#	}

		#}		


		## ------------- 4. update rho0: constant part --------- ## 
		if (constantRho == TRUE) {
			res <- Y - xfit - zfit - afit - fefit - yfit - rfit

			### 1. sample from proposed distribution
			srho <- samplesRho0(WY, res, Lambda, sigma2)
			## 2. whether to accept 
			#rho0 <- srho
			rho0 <- sampleRho0(WY, res, Lambda, rho0_old, srho, sigma2)

			#a0 <- abvalue0(W, Y, res, Lambda, rho0_old, sigma2)
			#a1 <- abvalue0(W, Y, res, Lambda, srho, sigma2)

			#ra <- a1 / a0
			#rv <- runif(1)
			#if (a0 == 0) {
			#	rho0 <- srho
			#} else {
			#	if (ra > rv) {
			#		rho0 <- srho
			#	} else {
			#		rho0 <- rho0_old
			#	}
			#}
			
			## 3. fitted 
			r0fit <- rhofit(WY, matrix(rho0, TT, 1))

			rho0_old <- rho0

			if (i > burn) {
				rho0_i[,i - burn] <- rho0
			}
			
		}

		## ------------- 4. update rho --------- ## 
		if (randomRho == TRUE) {
			res <- Y - xfit - zfit - afit - fefit - r0fit - yfit

			### 1. sample from proposed distribution
			Rho <- sampleRho(WY, res, Rho_old, Lambda, rhoZ, rhoA, type, kappa, sigma2, sigma_n2)
			## 2. calculate acceptance rate
			#arate <- genRate(res, sRho, Rho_old, Lambda, multilevel, kappa, sigma2, sigma_n2)
			## 3. gen new Rhos
			#Rho <- updateRho(sRho, Rho_old, rate)
			Rho_old <- Rho
			if (i > burn) {
				Rho_i[,i - burn] <- Rho
			}
			

			### 2. sample rhoA
			if (p4 > 0) {
				rhoA <- sampleRhoA(Rho, rhoZ, RA0, ra0, type, kappa, sigma_n2)
				if (i > burn) {
					rhoA_i[, i - burn] <- rhoA
				}
			}

			### 3. sample ar1 parameters
			if (type == 3) {
				kappa <- sampleKappa(Rho, rhoZ, rhoA, sigma_n2)
				if (i > burn) {
					kappa_i[,i - burn] <- kappa
				}
			}

			## 4. update state space error variance 
			sigma_n2 <- sampleSigman2(Rho, rhoZ, rhoA, type, kappa)
			if (i > burn) {
				sigma_n2_i[,i - burn] <- sigma_n2
			}
			
			## 6. update fitted values
			rfit <- rhofit(WY, Rho)
		}
		

		## ---------------- 5. update error term variance ------------ ##
		res <- Y - xfit - zfit - afit - r0fit - rfit - fefit - yfit
		sigma2 <- sampleSigma2(res)
		if (i > burn) {
			sigma2_i[,i - burn] <- sigma2
			## yfit_i[,, i - burn] <- Y - res
		}
		

		if (i %% 100 == 0) {
			cat(paste("Simulated rounds: ", i, sep = ""))
			cat("\n")
		}


	}

	out <- list(beta = beta_i, 
		        rhoy = rhoy_i,
		        omega = omega_i,
		        sigma2 = sigma2_i, 
		        Alpha = Alpha_i,
		        Xi = Xi_i,
		        rho0 = rho0_i,
		        Rho = Rho_i,
		        kappa = kappa_i,
		        sigma_n2 = sigma_n2_i,
		        rhoA = rhoA_i,
		        L = L_i,
		        Factor = Factor_i
		        )

	return(out)

} 




olsBeta <- function(Y, X) {

	TT <- dim(Y)[2]
	N <- dim(Y)[1]
	p <- dim(X)[2]

	xx <- matrix(0, p, p)
	xy <- matrix(0, p, 1)

	for (i in 1:TT) {
		subx <- X[,,i]
		suby <- as.matrix(Y[, i])

		xx <- xx + t(subx) %*% subx
		xy <- xy + t(subx) %*% suby
	}

	beta <- solve(xx) %*% xy
	return(beta)
}




