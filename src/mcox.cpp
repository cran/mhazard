#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix calcF(NumericVector lambda10, NumericVector lambda01,
		     NumericMatrix lambda11) {
  int I = lambda10.size();
  int J = lambda01.size();
  NumericMatrix Fout(I, J);
  int i;
  int j;

  Fout(0, 0) = 1;
  for (i=1; i<I; i++) {
    Fout(i, 0) = Fout(i-1, 0)*(1 - lambda10[i]);
  }
  for (j=1; j<J; j++) {
    Fout(0, j) = Fout(0, j-1)*(1 - lambda01[j]);
  }
  for (i=1; i<I; i++) {
    for (j=1; j<J; j++) {
      Fout(i, j) = Fout(i, j-1) + Fout(i-1, j) - Fout(i-1, j-1)*
	(1 - lambda11(i, j));
    }
    Rcpp::checkUserInterrupt();
  }
  return(Fout);
}

// [[Rcpp::export]]
NumericMatrix calc_lambda11(NumericVector T1, NumericVector T2,
			    NumericVector Y1, NumericVector Y2,
			    NumericMatrix X, NumericMatrix dij11,
			    NumericVector beta11) {
  int I = T1.size();
  int J = T2.size();
  int K = Y1.size();
  int P = beta11.size();
  NumericMatrix lambda11(I+1, J+1);
  int i;
  int j;
  int k;
  int p;
  double rsum;
  double fsum;
  
  std::fill(lambda11.begin(), lambda11.end(), 0);

  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      if (dij11(i+1, j+1)>0) {
	fsum = 0;
	for (k=0; k<K; k++) {
	  if (Y1[k]>=T1[i] && Y2[k]>=T2[j]) {
	    rsum = 0;
	    for (p=0; p<P; p++) {
	      rsum += (X(k, p) * beta11[p]);
	    }
	    fsum += exp(rsum);
	  }
	}
	lambda11(i+1, j+1) = dij11(i+1, j+1) / fsum;
      }
    }
    Rcpp::checkUserInterrupt();
  }
  return(lambda11);
}

// [[Rcpp::export]]
NumericVector calcEE11(NumericVector beta11, NumericVector Y1,
		       NumericVector Y2, NumericVector T1,
		       NumericVector T2, NumericVector Delta1,
		       NumericVector Delta2, NumericMatrix X,
		       NumericMatrix dij11) {
  int I = T1.size();
  int J = T2.size();
  int K = Y1.size();
  int P = beta11.size();
  int i;
  int j;
  int k;
  int p;
  NumericVector temp1(K);
  NumericVector temp2(P);
  double rsum;
  double temp3;
  NumericVector temp4(P);
  NumericVector EE11out(P);

  std::fill(EE11out.begin(), EE11out.end(), 0);

  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      if (dij11(i+1, j+1)>0) {
	std::fill(temp1.begin(), temp1.end(), 0);
	temp3 = 0;
	for (k=0; k<K; k++) {
	  if (Y1[k]>=T1[i] && Y2[k]>=T2[j]) {
	    rsum = 0;
	    for (p=0; p<P; p++) {
	      rsum += (X(k, p) * beta11[p]);
	    }
	    temp1[k] = exp(rsum);
	  }
	  temp3 += temp1[k];
	}
	std::fill(temp2.begin(), temp2.end(), 0);
	std::fill(temp4.begin(), temp4.end(), 0);
	for (p=0; p<P; p++) {
	  for (k=0; k<K; k++) {
	    temp2[p] += X(k, p)*temp1[k];
	    if (Y1[k]==T1[i] && Y2[k]==T2[j] && Delta1[k]==1 && Delta2[k]==1) {
	      temp4[p] += X(k, p);
	    }
	  }
	  EE11out[p] += temp4[p] - dij11(i+1, j+1) * temp2(p) / temp3;
	}
      }
    }
    Rcpp::checkUserInterrupt();
  }
  return(EE11out);
}

// [[Rcpp::export]]
arma::cube calc_drij(NumericVector Y1, NumericVector Y2,
		     NumericVector T1, NumericVector T2,
		     NumericVector Delta1, NumericVector Delta2) {
  int I = T1.size();
  int J = T2.size();
  int K = Y1.size();
  int i;
  int j;
  int k;
  arma::cube drout(I, J, 2);

  std::fill(drout.begin(), drout.end(), 0);
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      for (k=0; k<K; k++) {
	if (Y1[k]>=T1[i] && Y2[k]>=T2[j]) {
	  drout(i, j, 0)++;
	}
	if (Y1[k]==T1[i] && Y2[k]==T2[j] && Delta1[k]==1 && Delta2[k]==1) {
	  drout(i, j, 1)++;
	}
      }
      Rcpp::checkUserInterrupt();
    }
  }
  return(drout);
}
