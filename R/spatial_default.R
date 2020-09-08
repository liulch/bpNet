
bpNet <- function(data, ## sort by id, time 
	              W,    ## N * N * T
	              index, ## id and time 
	              Yname, ## outcome
	              Xname, ## covar
	              Zname, ## unit-level random effect
	              Aname, ## time-level random effect
	              Contextual = NULL, ## contextual effect
	              Contextual.effect = "time", ## random effect of contextual variable
	              lagY, ## lagged outcome, bool value
	              force, ## two-way fe
	              r = 0, ## factor number
	              flasso = 0, ## lasso
	              rhoZ, ## time-varying indicator
	              constantRho = TRUE, ## rho
	              randomRho = FALSE, ## rho_t
	              spec = "multilevel", ## rho_t specification
	              niter = 10000, 
	              burn = 1000) {

	library(abind)
	
	id <- index[1]
	time <- index[2]

	if (sum(is.na(data)) > 0) {
		data <- na.omit(data)
	}

	data[, id] <- as.numeric(as.factor(data[, id]))
	data[, time] <- as.numeric(as.factor(data[, time]))
	## data[, group] <- as.numeric(as.factor(data[, group]))

	if (var(table(data[, id])) != 0) {
		stop("Unbalanced panels.\n")
	}

	if (var(table(data[, time])) != 0) {
		stop("Unbalanced panels.\n")
	}

	if (force == "both") {
		force <- 3
	} else if (force == "unit") {
		force <- 1
	} else if (force == "time") {
		force <- 2
	} else {
		force <- 0
	}

	if (Contextual.effect == "both") {
		Contextual.effect <- 3
	} else if (force == "unit") {
		Contextual.effect <- 1
	} else if (force == "time") {
		Contextual.effect <- 2
	} else {
		Contextual.effect <- 0
	}


	## duplicate data
	data.old <- data

	data1 <- data.old[order(data.old[, id], data.old[, time]), ]
	data2 <- data.old[order(data.old[, time], data.old[, id]), ]

	N <- length(unique(data1[, id]))
	TT <- length(unique(data1[, time]))

	## group <- matrix(data2[, group], N, TT)

	if (!is.null(rhoZ)) {
		if (class(rhoZ) != "matrix") {
			stop("\"rhoZ\" should be a matrix.\n")
		} else {
			if (dim(rhoZ)[1] != TT) {
				stop("\"rhoZ\" should have length that equals the time periods.\n")
			}
		}
	}

	## gen matrix
	Y <- matrix(data2[, Yname], N, TT)

	p1 <- length(Xname)
	p2 <- length(Zname)
	p3 <- length(Aname)
	p4 <- length(Contextual)

	X <- Z <- A <- C <- C1 <- C2 <- NULL

	if (p1 > 0) {
		X <- array(NA, dim = c(N, p1, TT))
	}

	if (p2 > 0) {
		Z <- array(NA, dim = c(TT, p2, N))
	}

	if (p3 > 0) {
		A <- array(NA, dim = c(N, p3, TT))
	}

	if (p4 > 0) {
		C <- array(NA, dim = c(N, p4, TT))
		if (Contextual.effect %in% c(1, 3)) {
			C1 <- array(NA, dim = c(TT, p4, N))
		}
		#if (Contextual.effect %in% c(2, 3)) {
		#	C2 <- array(NA, dim = c(N, p4, TT))
		#}
	}

	if (p1 > 0 || p3 > 0 || p4 > 0) {
		for (i in 1:TT) {
			subdata <- data2[which(data2[, time] == i), ]
            if (p1 > 0) {
            	X[,, i] <- as.matrix(subdata[, Xname])
            }
            if (p3 > 0) {
            	A[,, i] <- as.matrix(subdata[, Aname])
            }
            if (p4 > 0) {
            	C[,, i] <- W[,, i] %*% as.matrix(subdata[, Contextual])  
            }
		}
	}

	if (p2 > 0) {
		for (i in 1:N) {
			subdata <- data1[which(data1[, id] == i), ]
			Z[,, i] <- as.matrix(subdata[, Zname])
		}
	}

	if (p4 > 0) {
		
		if (p1 > 0) {
			X <- abind(X, C, along = 2)
		} else {
			X <- C
		}
		
		if (Contextual.effect %in% c(1, 3)) {
			for (i in 1:p4) {
				C1[,i,] <- t(C[,i,])
			}
			if (p2 > 0) {
				Z <- abind(Z, C1, along = 2)
			} else {
				Z <- C1
			}

		}

		if (Contextual.effect %in% c(2, 3)) {
			if (p3 > 0) {
				A <- abind(A, C, along = 2)
			} else {
				A <- C
			}
			
		}

	}

	## add intercept
	if (!is.null(X)) {
		X <- abind(array(1, dim = c(N, 1, TT)), X, along = 2)
	} else {
		X <- array(1, dim = c(N, 1, TT))
	}

    if (force %in% c(1, 3)) {
    	if (!is.null(Z)) {
    		Z <- abind(array(1, dim = c(TT, 1, N)), Z, along = 2)
    	} else {
    		Z <- array(1, dim = c(TT, 1, N))
    	}

    }

    if (force %in% c(2, 3)) {
    	if (!is.null(A)) {
    		A <- abind(array(1, dim = c(N, 1, TT)), A, along = 2)
    	} else {
    		A <- array(1, dim = c(N, 1, TT))
    	}
    }

    ## lagged outcome
    lagy <- NULL
    if (lagY == TRUE) {
    	
    	lagy <- Y[, -TT]
    	Y <- Y[, -1]

    	if (!is.null(X)) {
    		X <- X[,,-1, drop = FALSE]
    	}

    	if (!is.null(Z)) {
    		Z <- Z[-1,,, drop = FALSE]
    	}

    	if (!is.null(A)) {
    		A <- A[,,-1, drop = FALSE]
    	}

    	if (!is.null(rhoZ)) {
    		rhoZ <- as.matrix(rhoZ[-1, ])
    	}

    	## group <- group[, -1]

    }

    sp.fit <- bspanel.core(Y = Y,  
	                       W = W,
		                   lagY = lagy,
		                   X = X,
		                   Z = Z,
		                   A = A,
		                   rhoZ = rhoZ, ## covariates used to interprete rho_t
		                   r = r,
		                   flasso = flasso,
		                   burn = burn,
		                   niter = niter,
		                   constantRho = constantRho,
		                   randomRho = randomRho, 
		                   spec = spec
		                   )

    return(sp.fit)

}






