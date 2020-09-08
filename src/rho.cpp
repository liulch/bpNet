# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp ;

/* ------------- 1. useful sub-functions ------------- */

// generate Ayt
// [[Rcpp::export]]
arma::mat GenAY (arma::cube W, // a N*N*T cube
	             arma::mat Y, // a N*T matrix
	             arma::mat Rho) { // a T*1 coloumn vector

  int T = Y.n_cols ;
  int N = Y.n_rows ;

  arma::mat AY(N, T, arma::fill::zeros) ;
  arma::mat suby(N, 1, arma::fill::zeros) ;

  arma::mat I(N, N, arma::fill::eye) ;
  arma::mat w(N, N, arma::fill::zeros) ;
  double p = 0 ;

  for (int i = 0; i < T; i++) {
  	p = Rho(i, 0) ;
  	w = W.slice(i) ;
  	suby = Y.col(i) ;

  	AY.col(i) = (I - p * w) * suby ;
  }

  return(AY) ;

}


// gen rho fit
// [[Rcpp::export]]
arma::mat rhofit(arma::mat WY,
	             arma::mat Rho) {

  int T = WY.n_cols ;
  int N = WY.n_rows ;

  arma::mat yfit(N, T, arma::fill::zeros) ;

  for (int i = 0; i < T; i++) {
    yfit.col(i) = Rho(i, 0) * WY.col(i) ;
  }

  return(yfit) ;
}


// gen beta fit
// [[Rcpp::export]]
arma::mat Xfit(arma::cube X,
	           arma::mat beta) {

	int N = X.n_rows ;
	int p = X.n_cols ;
	int T = X.n_slices ;

    arma::mat subx(N, p, arma::fill::zeros) ;
	arma::mat yfit(N, T, arma::fill::zeros) ;

	for (int i = 0; i < T; i++) {
		subx = X.slice(i) ;
		yfit.col(i) = subx * beta ;
	}

	return(yfit) ;

}

// gen unit multilevel fit
// [[Rcpp::export]]
arma::mat Zfit(arma::cube Z,      // T * p * N
	           arma::mat Alpha) { // p * N

	int T = Z.n_rows ;
	int p = Z.n_cols ;
	int N = Z.n_slices ;

    arma::mat subz(T, p, arma::fill::zeros) ;
    arma::mat subyfit(T, 1, arma::fill::zeros) ;

	arma::mat yfit(N, T, arma::fill::zeros) ;

	for (int i = 0; i < N; i++) {
		subz = Z.slice(i) ;
		subyfit = subz * Alpha.col(i) ;
		yfit.row(i) = subyfit.t() ;
	}

	return(yfit) ;

}


// gen time multilevel fit
// [[Rcpp::export]]
arma::mat Afit(arma::cube A,      // N * p * T
	           arma::mat Alpha) { // p * T

	int N = A.n_rows ;
	int p = A.n_cols ;
	int T = A.n_slices ;

    arma::mat suba(N, p, arma::fill::zeros) ;
	arma::mat yfit(N, T, arma::fill::zeros) ;

	for (int i = 0; i < T; i++) {
		suba = A.slice(i) ;
		yfit.col(i) = suba * Alpha.col(i) ;
	}

	return(yfit) ;

}


// gen B1 
// [[Rcpp::export]]
arma::mat genB1(arma::cube X, // a N*p*T cube
	            arma::mat invB0,
	            double sigma2) {

	int p = X.n_cols ;
	int T = X.n_slices ;
	int N = X.n_rows ;

	arma::mat XX(p, p, arma::fill::zeros) ;
	arma::mat subx(N, p, arma::fill::zeros) ;

	for (int i = 0; i < T; i++) {
        subx = X.slice(i) ;
		XX = XX + subx.t() * subx ;
	}

	return(arma::inv(invB0 + XX/sigma2)) ;
}


// gen bar_b
// [[Rcpp::export]]
arma::mat genb1(arma::cube X, // a N*p*T cube
	            arma::mat Y,  // a N*T matrix
	            arma::mat invB0,
	            arma::mat B1,
	            arma::mat b0,
	            double sigma2) {

	int p = X.n_cols ;
	int N = X.n_rows ;
	int T = X.n_slices ;

	arma::mat subx(N, p, arma::fill::zeros) ;
	arma::mat suby(N, 1, arma::fill::zeros) ;

	arma::mat b1 = invB0 * b0 ;

	for (int i = 0; i < T; i++) {
		subx = X.slice(i) ;
		suby = Y.col(i) ;
		b1 = b1 + subx.t() * suby / sigma2 ;
	}

	b1 = B1 * b1 ;

	return(b1) ;
}


// sample multivariate normal
// [[Rcpp::export]]
arma::mat sampleBeta(arma::cube X,
                     arma::mat Y,
                     arma::mat B0,
                     arma::mat b0,
                     double sigma2) {  

  // arma::mat invCov0 = arma::inv(Cov0) ;
  // invCov0(0, 0) = 0 ;

  int p = X.n_cols ;

  arma::mat invB0 = arma::inv(B0) ;

  arma::mat B1 = genB1(X, invB0, sigma2) ;
  arma::mat b1 = genb1(X, Y, invB0, B1, b0, sigma2) ;

  arma::mat beta(p, 1, arma::fill::zeros) ;
  if (p > 1) {
    beta = arma::mvnrnd(b1.col(0), B1) ;
  }
  else {
    beta(0, 0) = R::rnorm(b1(0, 0), sqrt(B1(0, 0))) ; 
  }

  return(beta) ;

}


/* ----------- random effect at unit and time level -------------- */
// gen A1 
// [[Rcpp::export]]
arma::mat genA1(arma::mat X, // a T*p matrix
	            arma::mat invA0,
	            double sigma2) {

	// int p = X.n_cols ;
	// int N = X.n_rows ;

	arma::mat XX = X.t() * X ;

	return(arma::inv(invA0 + XX/sigma2)) ;
}


// gen bar_a
// [[Rcpp::export]]
arma::mat gena1(arma::mat X,  // a T*p matrix
	            arma::mat Y,  // a T*1 matrix
	            arma::mat A1,
	            double sigma2) {

	int p = X.n_cols ;
	// int N = X.n_rows ;

	arma::mat a1(p, 1, arma::fill::zeros) ;
	a1 = X.t() * Y / sigma2 ;
	a1 = A1 * a1 ;

	return(a1) ;
}


// sample conditional multivariate normal
// [[Rcpp::export]]
arma::mat sampleCN(arma::mat M,   // mean vector
                   arma::mat S,   // covariance matrix
                   arma::mat m) { // conditional on outcome

  int n1 = M.n_rows ;
  int n2 = m.n_rows ;
  int n3 = n1 - n2 ;

  arma::mat s11 = S.submat(0, 0, n3 - 1, n3 - 1) ;
  arma::mat s22 = S.submat(n3, n3, n1 - 1, n1 - 1) ;
  arma::mat s12 = S.submat(0, n3, n3 - 1, n1 - 1) ;
  arma::mat s21 = S.submat(n3, 0, n1 - 1, n3 - 1) ;

  arma::mat mu1 = M.rows(0, n3 - 1) ;
  arma::mat mu2 = M.rows(n3, n1 - 1) ;

  arma::mat invs22 = arma::inv(s22) ;
  arma::mat mu = mu1 + s12 * invs22 * (m - mu2) ;
  arma::mat s = s11 - s12 * invs22 * s21 ;

  arma::mat sv(n3, 1) ;
  sv.fill(-1) ;

  if (n3 == 1) {
    // while(arma::accu(sv - arma::abs(sv)) != 0) {
      sv(0, 0) = R::rnorm(mu(0, 0), sqrt(s(0, 0))) ;
    // }
  } 
  else {
    // while(arma::accu(sv - arma::abs(sv)) != 0) {
      sv = arma::mvnrnd(mu.col(0), s) ;
    // }
  }

  return(sv) ;

}


// sampling fixed effects
// [[Rcpp::export]]
arma::mat sampleVector(arma::mat M,  // mean vector
                       arma::mat S, // covariance matrix
                       int id, 
                       int r 
                      ) {

  int p = M.n_rows ;
  arma::mat sV(p, 1, arma::fill::zeros) ;

  if (id <= r - 1) {
  	// sV(p - r - 1 + id, 0) = 1 ;
  	sV.rows(0, p - r - 1 + id) = sampleCN(M, S, sV.rows(p - r + id, p - 1)) ;
  } else {
  	sV = arma::mvnrnd(M.col(0), S) ;
  }

  return(sV) ;

}




// sample multivariate normal
// [[Rcpp::export]]
arma::mat sampleSubAlpha(arma::mat X, // a N*p matrix
                         arma::mat Y, // a N*1 matrix
                         arma::mat A0,
                         int unit, // unit or time level random effetcs
                         int id,   // pos of the id
                         int r,    // factor numbers
                         double sigma2) {  

	// arma::mat invCov0 = arma::inv(Cov0) ;
	// invCov0(0, 0) = 0 ;

	int p = X.n_cols ;

	arma::mat invA0 = arma::inv(A0) ;

	arma::mat A1 = genA1(X, invA0, sigma2) ;
	arma::mat a1 = gena1(X, Y, A1, sigma2) ;

	arma::mat alpha(p, 1, arma::fill::zeros) ;
	  
	if (p > 1) {
		if (r > 0 && unit == 1) {
			alpha = sampleVector(a1, A1, id, r) ;
		} else {
			alpha = arma::mvnrnd(a1.col(0), A1) ;
		}
	}
	else {
	    if (r == 1 && unit == 1) {
	    	// if (id == 1) {
	    	//	alpha(0 ,0) = 1;
	    	// } else {
	    		alpha(0, 0) = R::rnorm(a1(0, 0), sqrt(A1(0, 0))) ;
	    	// }
	    } else {
	    	alpha(0, 0) = R::rnorm(a1(0, 0), sqrt(A1(0, 0))) ; 
	    }
	}

	return(alpha) ;

}


// sample random effects
// [[Rcpp::export]]
arma::mat sampleAlpha(arma::cube X,  // 
	                  arma::mat Y,   // N*T
	                  arma::mat A0,
	                  int unit, // unit or time level random effetcs
                      int r,    // factor numbers    
	                  double sigma2) {

	
	int N = Y.n_rows ;
	// int T = Y.n_cols ;
	
	int p = X.n_cols ;
	int n1 = X.n_slices ;
	int n2 = X.n_rows ;

	int stime = 1 ;
	if (n1 == N) {
		stime = 0 ; 
	}

	arma::mat subx(n2, p, arma::fill::zeros) ;
	arma::mat suby(n2, 1, arma::fill::zeros) ;

	arma::mat Alpha(p, n1, arma::fill::zeros) ;

	for (int i = 0; i < n1; i++) {
        subx = X.slice(i) ;
        if (stime == 1) {
        	suby = Y.col(i) ;
        } 
        else if (stime == 0) {
        	suby = Y.row(i).t() ;
        }

		Alpha.col(i) = sampleSubAlpha(subx, suby, A0, unit, i+1, r, sigma2) ;
	}

	return(Alpha) ;
}



/* ---------------- 2. variance in prior --------------- */
// update prior
// [[Rcpp::export]]
arma::mat genD1(arma::mat D0, arma::mat coef) {

	int p = coef.n_rows ;
	int n = coef.n_cols ;
	arma::mat D1(p, p, arma::fill::zeros) ;
	arma::mat f(1, p, arma::fill::zeros) ;

	for (int i = 0; i < n; i++) {
	    f = coef.col(i) ;   
	    D1 = D1 + f * f.t() ;
	}
	  
	D1 = arma::inv(arma::inv(D0) + D1) ;
	return(D1) ;

}

// sample covariance matrix
// [[Rcpp::export]]
arma::mat sampleD(arma::mat D0, 
                  arma::mat coef,
                  double d0) {

	int p = coef.n_rows ;
	int n = coef.n_cols ;
	arma::mat D1 = genD1(D0, coef) ;
	double d1 = d0 ;
	  
	d1 = d1 + n ;
	  
	arma::mat D(p, p, arma::fill::zeros) ;

	if (p > 1) {
	    D = arma::wishrnd(D1, d1) ;
	}
	else {
	    D(0, 0) = arma::randg<double>( arma::distr_param(d1 / 2, 1 / (D1(0, 0) / 2)) ) ;
	}
	return(arma::inv(D)) ;
}



                   /* ----------------------- */
/* ---------------- MH Algorithm1: constant rho ----------------- */
                   /* ----------------------- */

// gen Phi
// [[Rcpp::export]]         
double genPhi0(arma::mat WY,
	           double sigma2) {

	int T = WY.n_cols ;
	int N = WY.n_rows ;

	double a = 0 ;

	// arma::mat subw(N, N, arma::fill::zeros) ;
	// arma::mat suby(N, 1, arma::fill::zeros) ;
	arma::mat wy(N, 1, arma::fill::zeros) ;
	arma::mat pwy(1, 1, arma::fill::zeros) ;

	for (int i = 0; i < T; i++) {
		// subw = W.slice(i) ;
		// suby = Y.col(i) ;
		// wy = subw * suby ;
		wy = WY.col(i) ;
		pwy = wy.t() * wy ;
 
		a = a + pwy(0, 0) ;
	}

	// a = 1 / (a / sigma2) ;
	a = sigma2 / a ; 
	return(a) ;
}


// gen phi
// [[Rcpp::export]]
double genphi0(arma::mat WY,
	           arma::mat res,
	           double Phi0,
	           double sigma2) {

	int T = WY.n_cols ;
	int N = WY.n_rows ;

	double a = 0 ;

	//arma::mat subw(N, N, arma::fill::zeros) ;
	//arma::mat suby(N, 1, arma::fill::zeros) ;
	arma::mat wy(N, 1, arma::fill::zeros) ;

	arma::mat subres(N, 1, arma::fill::zeros) ;
	
	arma::mat prwy(1, 1, arma::fill::zeros) ;

	for (int i = 0; i < T; i++) {
		// subw = W.slice(i) ;
		// suby = Y.col(i) ;
		wy = WY.col(i) ;
		subres = res.col(i) ;

		prwy = subres.t() * wy ;
		
		a = a + prwy(0, 0) ;
	}

	a = a / sigma2 * Phi0 ;
	return(a) ;
	
}


// sample Rho
// [[Rcpp::export]]
double samplesRho0(arma::mat WY,
	               arma::mat res,
	               arma::mat lambda,
	               double sigma2) {


	// int T = WY.n_cols ;
	// int N = W.n_rows ;

	//  arma::mat minlambda(T, 1, arma::fill::zeros) ;
	// arma::mat maxlambda(T, 1, arma::fill::zeros) ;

	// arma::mat subl(N, 1, arma::fill::zeros) ;

	// for (int i = 0; i < T; i++) {
	//	minlambda(i, 0) = lambda.col(i).min() ;
	//	maxlambda(i, 0) = lambda.col(i).max() ;
	// }

	// arma::mat iminlambda = 1 / minlambda ;
	// arma::mat imaxlambda = 1 / maxlambda ;

	// double a1 = iminlambda.max() ;
	// double a2 = imaxlambda.min() ;

	double P0 = 0 ;
	double p0 = 0 ;

	// update rho0
	double rho0 = 0 ;

	P0 = genPhi0(WY, sigma2) ;

    p0 = genphi0(WY, res, P0, sigma2) ;

    rho0 = R::rnorm(p0, sqrt(P0)) ;

    // while (rho0 >= a2 || rho0 <= a1) {
	// 	rho0 = R::rnorm(p0, sqrt(P0)) ;
	// }

	return(rho0) ;
}


// acceptance rate: likelihood
// [[Rcpp::export]]
double abvalue0(arma::mat WY,
	            arma::mat res,
	            arma::mat lambda,
	            double rho,
	            double sigma2) {

	int m = lambda.n_rows ;
	int T = lambda.n_cols ;

	double a = 1 ;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < T; j++) {
			a = a * (1 - rho * lambda(i, j)) ;
		}		
	}

	double b = 0 ;

    int N = res.n_rows ;
	arma::mat subres(N, 1, arma::fill::zeros) ;

	arma::mat pror(1, 1, arma::fill::zeros) ;

	for (int i = 0; i < T; i++) {
		subres = res.col(i) - WY.col(i) * rho ;
		pror = subres.t() * subres ;
		b = b + pror(0, 0) ;
		// b = b + 2 * log(abs(subres(i, 0))) ;
	}

	b = -b / (2 * sigma2) ;

	// double c = a * exp(b) ;
	double c = 1 ;
	if (a > 0) {
		// c = (log(a) + b) / (N * T) ;
		c = log(a) + b ;
	}

	return(c) ;

}



// calculate acceptance rate
// [[Rcpp::export]]
double sampleRho0(arma::mat WY,
	              arma::mat res,
	              arma::mat lambda,
	              double rho0,
	              double rho1,
	              double sigma2) {

	double a0 = abvalue0(WY, res, lambda, rho0, sigma2) ;
	double a1 = abvalue0(WY, res, lambda, rho1, sigma2) ;

	// double rate = 0 ;
	double rho = 0 ;
	double rv = 2 ;

	if (a1 == 1) {
		rho = rho0 ;
	} else {
		rv = R::runif(0, 1) ;

		if (a0 == 1 || rv == 0) {
			rho = rho1 ;
		} else {
			if (a1 - a0 >= log(rv)) {
				rho = rho1 ;
			} else {
				rho = rho0 ;
			}
		}		

		/*if (a1 >= a0) {
			rho = rho1 ;
		} else {
			rv = R::runif(0, 1) ;
			if (rv == 0) {
				rho = rho1 ;
			} else {
				if (a1 - a0 >= log(rv)) {
					rho = rho1 ;
				} else {
					rho = rho0 ;
				}
			}
		}*/
	}

	/*while (a1 == 0 || a0 == 0) {
		a0 = a0 * 10000000000 ;
	    a1 = a1 * 10000000000 ;
	}

	if (a0 == 0) {
		rho = rho1 ;
	} else {
		
		rate = a1 / a0 ;

		if (rate >= 1) {
			rho = rho1 ;
			// update = 1 ;
		}
		else {
			rv = R::runif(0, 1) ;
			if (rv < rate) {
				rho = rho1 ;
			    // update = 1 ;
			}
			else {
				rho = rho0 ;
			}
		}

	}*/

	

	// double rate = a1 / a0 ;

	// double rho = 0 ;
	// double update = 0 ;

	

	



	// if (a1 < 0) {
	// 	rho = rho0 ;
	// } else {
	//	if (a0 == 0) {
	//		if (a1 == 0) {
	//		  a0 = a0 * 1000000 ;
	//		  a1 = a1 * 1000000 ;
	//		  if (a0 > a1) {
	//		  	rho = rho0 ;
	//		  } else {
	//		  	rho = rho1 ;
	//		  }	
	//		} else {
	//		  rho = rho1 ;
	//		}	
	//	}
	//	else {
	//		if (rate >= 1) {
	//			rho = rho1 ;
				// update = 1 ;
	//		}
	//		else {
	//			rv = R::runif(0, 1) ;
	//			if (rv < rate) {
	//				rho = rho1 ;
	//			    // update = 1 ;
	//			}
	//			else {
	//				rho = rho0 ;
	//			}
	//		}

	//	}
	// }
	
	// List out;

	// out["rho"] = rho ;
	// out["update"] = update ;

	return(rho) ;

}



                   /* ----------------------- */
/* ---------------- MH Algorithm2: random rho_t ----------------- */
                   /* ----------------------- */

// we consider three types: type1: random walk, type2: multilevel, type3: ar1

// gen Phi
// [[Rcpp::export]]
double genPhi(arma::mat wy,
	          int type, 
	          int time,
	          int T,
	          double kappa,
	          double sigma2,
	          double sigma_n2) {

	// arma::mat wy = W * Y ;
	arma::mat pwy = wy.t() * wy ;

	double Phi = 0 ;
	if (type == 2) {
		Phi = 1 / (pwy(0, 0) / sigma2 + 1 / sigma_n2) ;
	} else {
		if (time == 1) {
			Phi = 1 / (pwy(0, 0) / sigma2 + (1 + pow(kappa, 2)) / sigma_n2) ;
		}
		else if (time < T) {
			Phi = 1 / (pwy(0, 0) / sigma2 + (1 + pow(kappa, 2)) / sigma_n2) ;
		} else {
			Phi = 1 / (pwy(0, 0) / sigma2 + 1 / sigma_n2) ;
		}
	}
	return(Phi) ;
}


// gen phi
// [[Rcpp::export]]
double genphi(arma::mat wy,
	          arma::mat res,
	          arma::mat Rho_old,
	          arma::mat Rho_new,
	          arma::mat rhoZ, // covariates
	          arma::mat rhoA, // covar coefficients
	          int type, 
	          int time,
	          int T,
	          double kappa,
	          double Phi,
	          double sigma2,
	          double sigma_n2) {

	// arma::mat wy = W * Y ;

	arma::mat prow = wy.t() * res ;

	arma::mat m1(1, 1, arma::fill::zeros) ;
	arma::mat m2(1, 1, arma::fill::zeros) ;

	double a = prow(0, 0) / sigma2 ;

	int p = rhoZ.n_cols ;


	// Rcpp::Rcout << Rho_new.n_rows << std::endl;
	// Rcpp::Rcout << time << ": " << Phi << std::endl;
	if (type == 2) { // multi-level
		if (p > 0) {
			m1 = rhoZ.row(time - 1) * rhoA ;
		}
		a = a + m1(0, 0) / sigma_n2 ;
	} else {
		if (time == 1) {
			if (p > 0) {
				m1 = kappa * (Rho_old(time, 0) - rhoZ.row(time) * rhoA) ;
				m2 = rhoZ.row(time - 1) * rhoA ;
			}
			a = a + (m1(0, 0) + m2(0, 0)) / sigma_n2 ;
		} else if (time < T) {
			if (p > 0) {
				m1 = kappa * (Rho_old(time, 0) - rhoZ.row(time) * rhoA) ;
			    m2 = kappa * Rho_new(time - 2, 0) + rhoZ.row(time - 1) * rhoA ;
			}
			a = a + (m1(0, 0) + m2(0, 0)) / sigma_n2 ; 
		} else {
			if (p > 0) {
			    m2 = kappa * Rho_new(time - 2, 0) + rhoZ.row(time - 1) * rhoA ;
		    }
			a = a +  m2(0, 0)/ sigma_n2 ;
		}
	}

	// if (time == 1) {
	//	a = a + kappa *  Rho_old(1, 0) / sigma_n2 ;
	// }
	// else if (time < T) {
	// 	a = a + kappa * (Rho_new(time - 2, 0) + Rho_old(time, 0)) / sigma_n2 ;
	// }
	// else {
	//	Rcpp::Rcout << a << std::endl;
	// 	a = a + kappa * Rho_new(time - 2, 0) / sigma_n2 ;
	// } 

	a = a * Phi ;
	// Rcpp::Rcout << time << ": " << a << std::endl;
	return(a) ;
	
}

// sample rho_t for a given t
// [[Rcpp::export]]
double sampleSubRho(arma::mat wy,
			        arma::mat res,
			        arma::mat lambda,
			        arma::mat Rho_old,
			        arma::mat Rho_new,
			        arma::mat rhoZ, // covariates
	                arma::mat rhoA, // covar coefficients
			        int type,
			        int time,
		            int T,
			        double kappa,
			        double sigma2,
			        double sigma_n2) {

	double Phi = genPhi(wy, type, time, T, kappa, sigma2, sigma_n2) ;
	double phi = genphi(wy, res, Rho_old, Rho_new, rhoZ, rhoA, type, time, T, kappa, Phi, sigma2, sigma_n2) ;

    // double lmin = lambda.min() ;
    // double lmax = lambda.max() ;

    double srho = R::rnorm(phi, sqrt(Phi)) ;

	// while (srho >= 1/lmax || srho <= 1/lmin) {
	// 	srho = R::rnorm(phi, sqrt(Phi)) ;
	// }

	return(srho) ;
}


// acceptance rate: likelihood
// [[Rcpp::export]]
double abvalue(arma::mat WY,
	           arma::mat res,
	           arma::mat Rho_old,
	           arma::mat Rho_new,
	           arma::mat lambda,
	           arma::mat rhoZ,
	           arma::mat rhoA,
	           int type,
	           int time,
	           int T,
	           double rho,
	           double kappa,
	           double sigma2,
	           double sigma_n2) {

	int m = lambda.n_rows ;
	int p = rhoZ.n_cols ;

	double a = 1 ;
	for (int i = 0; i < m; i++) {
		a = a * (1 - rho * lambda(i, 0)) ;
	}

	double b = 0 ;

	arma::mat new_res = res - WY.col(time - 1) * rho ;

	arma::mat pror = new_res.t() * new_res ;
	b = - pror(0, 0) / (2 * sigma2) ;

	double c = 0 ;
	arma::mat m1(1, 1, arma::fill::zeros) ;
	arma::mat m2(1, 1, arma::fill::zeros) ;

	if (type == 2) {
		if (p > 0) {
			m1 = rhoZ.row(time - 1) * rhoA ;
		}
		c = - pow(rho - m1(0, 0), 2) / (2 * sigma_n2) ;
	} else {
		if (time == 1) {
			if (p > 0) {
				m1 = rhoZ.row(time) * rhoA ;
			}			
			c = - (pow(Rho_old(time, 0) - m1(0, 0) - kappa * rho, 2)) / (2 * sigma_n2) ;
		}
		else if (time < T) {
			if (p > 0) {
				m1 = rhoZ.row(time) * rhoA ;
			    m2 = rhoZ.row(time - 1) * rhoA ;
			}

			c = - (pow(Rho_old(time, 0) - m1(0, 0) - kappa * rho, 2) + pow(rho - kappa * Rho_new(time - 2, 0) - m2(0, 0), 2)) / (2 * sigma_n2) ;
		}
		else {
			if (p > 0) {
				m2 = rhoZ.row(time - 1) * rhoA ;
			}

			c = - pow(rho - kappa * Rho_new(time - 2, 0) - m2(0, 0), 2) / (2 * sigma_n2) ;
		}

	}
	
	// double d = a * exp(b + c) ;
	double d = 1 ;
	if (a > 0) {
		d = log(a) + b + c ;
	}
	
	return(d) ;

}


// sample Rho
// [[Rcpp::export]]
arma::mat sampleRho(arma::mat WY,
	                arma::mat res,
	                arma::mat Rho_old,
	                arma::mat lambda,
	                arma::mat rhoZ,
	                arma::mat rhoA,
	                int type,
	                double kappa,
	                double sigma2,
	                double sigma_n2) {


	int T = WY.n_cols ;
	int N = WY.n_rows ;

	arma::mat Rho(T, 1, arma::fill::zeros) ;

	arma::mat wy(N, 1, arma::fill::zeros) ;
	arma::mat subres(N, 1, arma::fill::zeros) ;
	arma::mat sublambda(N, 1, arma::fill::zeros) ;
	// arma::mat subw(N, N, arma::fill::zeros) ;

    double srho = 0 ;

    double l1 = 0 ;
    double l2 = 0 ;

    // double rate = 0 ;
    double rv = 2 ;

	for (int i = 0; i < T; i++) {

		// subw = W.slice(i) ;
		wy = WY.col(i) ;
		subres = res.col(i) ;
		sublambda = lambda.col(i) ;


		srho = sampleSubRho(wy, subres, sublambda, Rho_old, Rho, rhoZ, rhoA, type, i + 1, T, kappa, sigma2, sigma_n2) ;
        
        l2 = abvalue(WY, subres, Rho_old, Rho, sublambda, rhoZ, rhoA, type, i + 1, T, srho, kappa, sigma2, sigma_n2) ;
        l1 = abvalue(WY, subres, Rho_old, Rho, sublambda, rhoZ, rhoA, type, i + 1, T, Rho_old(i, 0), kappa, sigma2, sigma_n2) ;
    
        // l1 = l1 * 10000000 ;
		// l2 = l2 * 10000000 ;

		if (l2 == 1) {
			Rho(i, 0) = Rho_old(i, 0) ;
		} else {
			rv = R::runif(0, 1) ;
			if (l1 == 1 || rv == 0) {
				Rho(i, 0) = srho ;
			} else {
				if (l2 - l1 >= log(rv)) {
					Rho(i, 0) = srho ;
				} else {
					Rho(i, 0) = Rho_old(i, 0) ;
				}
			}

		}



		/*while (l1 == 0 || l2 == 0) {
			l1 = l1 * 100000 ;
		    l2 = l2 * 100000 ;
		}

		if (l1 == 0) {
			Rho(i, 0) = srho ;
		} else {
			
			rate = l2 / l1 ;

			if (rate >= 1) {
				Rho(i, 0) = srho ;
				// update = 1 ;
			}
			else {
				rv = R::runif(0, 1) ;
				if (rv <= rate) {
					Rho(i, 0) = srho ;
				    // update = 1 ;
				}
				else {
					Rho(i, 0) = Rho_old(i, 0) ;
				}
			}

	    }*/

        // rate = l1 / l2 ;
        // rv = R::runif(0, 1) ;

        // if (rv <= rate) {
    	//    Rho(i, 0) = srho ;
        // }
        // else {
        //    Rho(i, 0) = Rho_old(i, 0) ;
        // }

        //if (l1 < 0) {
        //	Rho(i, 0) = Rho_old(i, 0) ;
        //} else {
        //	if (l2 == 0) {
	    //    	if (l1 == 0) {
	    //    		l1 = l1 * 1000000 ;
		//			l2 = l2 * 1000000 ;
					
		//			if (l1 > l2) {
		//			  	Rho(i, 0) = Rho_old(i, 0) ;
		//			} else {
		//			  	Rho(i, 0) = srho ;
		//			}
	        		
	    //    	} else {
	    //    		Rho(i, 0) = srho ;
	    //    	}
	    //    }
	    //    else {
	    //    	if (rv <= rate) {
	    //    	    Rho(i, 0) = srho ;
	    //        }
	    //        else {
	    //            Rho(i, 0) = Rho_old(i, 0) ;
	    //        }
	    //    }

        //}
        
        
		// while (srho >= 1/lmax || srho <= 1/lmin) {
		// 	srho = R::rnorm(phi, sqrt(Phi)) ;
		// }
		
	}

	return(Rho) ;
}

// [[Rcpp::export]]
arma::mat sampleRhoA(arma::mat Rho,
	                 arma::mat rhoZ,
	                 arma::mat RA0,
	                 arma::mat ra0,
	                 int type,
	                 double kappa,
	                 double sigma_n2) {
	
	int T = Rho.n_rows ;
	int p = rhoZ.n_cols ;

	arma::mat RA1(p, p, arma::fill::zeros) ;
	arma::mat ra1(p, 1, arma::fill::zeros) ;

	arma::mat rhoA(p, 1, arma::fill::zeros) ;

	arma::mat Rho1 = Rho ;
	arma::mat Rho2 = Rho ;
	arma::mat difRho(T-1, 1, arma::fill::zeros) ;
	arma::mat rhoZ2 = rhoZ ;

	arma::mat srho1 = Rho.row(0) ;

	if (type == 2) {
		RA1 = arma::inv(arma::inv(RA0) + rhoZ.t() * rhoZ / sigma_n2) ;
		ra1 = RA1 * (arma::inv(RA0) * ra0 + rhoZ.t() * Rho / sigma_n2) ;
	} else {
		Rho1.shed_row(T-1) ;
		Rho2.shed_row(0) ;
		difRho = Rho2 - Rho1 * kappa ;
		difRho.insert_rows(0, srho1) ;

		// rhoZ2.shed_row(0) ;

		RA1 = arma::inv(arma::inv(RA0) + rhoZ2.t() * rhoZ2 / sigma_n2) ;
		ra1 = RA1 * (arma::inv(RA0)*ra0 + rhoZ2.t() * difRho / sigma_n2) ;
	}

	if (p > 1) {
		rhoA = arma::mvnrnd(ra1.col(0), RA1) ;
	} else {
		rhoA(0, 0) = R::rnorm(ra1(0, 0), sqrt(RA1(0, 0))) ; 
	}

	return(rhoA) ;
}


// sample phi for each ar(1) process
// [[Rcpp::export]]
double sampleKappa(arma::mat Rho,  // a T*1 matrix
	               arma::mat rhoZ,
	               arma::mat rhoA,
	               double sigma_n2) { 
  
	int T = Rho.n_rows ;
	int p = rhoZ.n_cols ;
	  
	arma::mat lagRho = Rho ;
	lagRho.shed_row(T-1) ;

	arma::mat Rho1 = Rho ;
	Rho1.shed_row(0) ;

	arma::mat rhoZ2 = rhoZ ;

	arma::mat res = Rho1 ;
	
	if (p > 0) {
		rhoZ2.shed_row(0) ;
		res = res - rhoZ2 * rhoA ;
	}

	double kappa = 0 ;

	double muphi = 0 ;
	double varphi = 0 ;

	arma::mat xx = lagRho.t() * lagRho ;
	arma::mat xy = lagRho.t() * res ;

	varphi = sigma_n2 / xx(0, 0) ;
	muphi = varphi * xy(0, 0) / sigma_n2 ;

	kappa = R::rnorm(muphi, sqrt(varphi)) ;

	// while(kappa < -1 || kappa > 1) {
	//   	kappa = R::rnorm(muphi, sqrt(varphi)) ;
	// }

	return(kappa) ;

}


// sample variance of error term for each ar(1) process
// [[Rcpp::export]]
double sampleSigman2(arma::mat Rho,  // a T*1 matrix
	                 arma::mat rhoZ,
	                 arma::mat rhoA,
	                 int type,
	                 double kappa) {

	int T = Rho.n_rows ;
	int p = rhoZ.n_cols ;
	  
	arma::mat lagRho = Rho;
	lagRho.shed_row(T-1) ;

	arma::mat Rho1 = Rho ;
	Rho1.shed_row(0) ;

	//arma::mat rhoZ2 = rhoZ ;
	//if (p > 0) {
	//	rhoZ2.shed_row(0) ;
	//}

	arma::mat ss1 = Rho.row(0) ;


	arma::mat res1(T, 1, arma::fill::zeros) ;
	arma::mat res2(T - 1, 1, arma::fill::zeros) ;

	double a = 0 ;
	double b = 0 ;
	double c = 0 ;

	if (type == 2) {
		res1 = Rho ;
		if (p > 0) {
			res1 = res1 - rhoZ * rhoA ;
		}
		a = (T + 0.001) / 2 ;
	    b = 1 / ((arma::accu(pow(res1, 2)) + 0.001) / 2) ;
	} else {
		res2 = Rho1 - kappa * lagRho ;
		res2.insert_rows(0, ss1) ;

		if (p > 0) {
			res2 = res2 - rhoZ * rhoA ;
		}
		a = (T + 0.001) / 2 ;
	    b = 1 / ((arma::accu(pow(res2, 2)) + 0.001) / 2) ;
	}

	c = 1 / arma::randg<double>( arma::distr_param(a, b) ) ;

	return(c) ;
} 



// sample variance of error term 
// [[Rcpp::export]]
double sampleSigma2(arma::mat res) {

	int T = res.n_cols ;
	int N = res.n_rows ;


	double a = (N * T + 0.001) / 2 ;
	double b = 1 / ((arma::accu(pow(res, 2)) + 0.001) / 2) ;

	double c = 1 / arma::randg<double>( arma::distr_param(a, b) ) ;

	return(c) ;
} 


// permutation: generate binary value
// [[Rcpp::export]]
arma::mat genBi(int k) {
  arma::mat Bi(k, 1, arma::fill::ones) ;
  double rn = 0 ;
  for (int i = 0; i < k; i++) {
    rn = R::rnorm(0.0, 1.0) ;
    if (rn < 0) {
      Bi(i, 0) = -1 ;
    }
  }
  return(Bi) ;
}


// permutation result
// [[Rcpp::export]]
List permute(arma::mat omega,
             arma::mat xi) {
  
  int k = xi.n_cols ;
  int T = xi.n_rows ;
  arma::mat coef = genBi(k) ; // k*1

  arma::mat omega_p = omega % coef ;
  arma::mat xi_p = xi % repmat(coef.t(), T, 1) ;

  List result ;
  result["omega"] = omega_p ;
  result["Xi"] = xi_p ;
  return(result) ;
}

// [[Rcpp::export]]
double sampleG(double a, double b) {
  

  double c = arma::randg<double>( arma::distr_param(a, b) ) ;


  return(c) ; 
}

// sample inv-gaussian
// [[Rcpp::export]]
double rrinvgauss(double mu, double lambda) {
 
  double z, y, x, u, m ;
 
  z = R::rnorm(0,1) ;
  y = z * z ;
  x = mu + 0.5 * mu * mu * y / lambda - 0.5 * (mu / lambda) * sqrt(4 * mu * lambda * y + mu * mu * y * y) ;
  u = R::runif(0, 1) ;
  if(u <= mu / (mu + x)) {
    m = x ;
  } else {
    m = mu * mu / x ;
  } 

  return(m) ;
}













