#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector margHaz(NumericVector Y, NumericVector T, NumericVector Delta) {
  int nY = Y.size();
  int nT = T.size();
  NumericVector Fhat(nT+1);
  int di;
  int ri;

  Fhat[0] = 1;
  for (int i=0; i<nT; i++) {
    di = 0;
    ri = 0;
    for (int j=0; j<nY; j++) {
      if (Y[j]==T[i] && Delta[j]==1) {
	di++;
      }
      if (Y[j]>=T[i]) {
	ri++;
      }
    }
    Fhat[i+1] = Fhat[i]*(1-double(di)/double(ri));
  }
  return(Fhat);
}

void calc_dij(NumericVector Y1, NumericVector Y2, NumericVector T1,
	      NumericVector T2, NumericVector Delta1, NumericVector Delta2,
	      NumericMatrix &dij00, NumericMatrix &dij01, NumericMatrix &dij10,
	      NumericMatrix &dij11, NumericMatrix &rij) {
  int nY = Y1.size();
  int nT1 = T1.size();
  int nT2 = T2.size();
  int i;
  int j;
  int k;

  for (i=0; i<nT1; i++) {
    for (j=0; j<nT2; j++) {
      for (k=0; k<nY; k++) {
	if (Y1[k]>=T1[i] && Y2[k]>=T2[j]) {
	  rij(i, j)++;
	}
	if (Y1[k]==T1[i] && Y2[k]==T2[j] && Delta1[k]==1 && Delta2[k]==1) {
	  dij11(i, j)++;
	}
	if (Y1[k]==T1[i] && Delta1[k]==1 &&
	    (Y2[k]>T2[j] ||
	     (Y2[k]==T2[j] && Delta2[k]==0))) {
	  dij10(i, j)++;
	}
	if (Y2[k]==T2[j] && Delta2[k]==1 &&
	    (Y1[k]>T1[i] ||
	     (Y1[k]==T1[i] && Delta1[k]==0))) {
	  dij01(i, j)++;
	}
      }
      dij00(i, j) = rij(i, j) - dij11(i, j) - dij10(i, j) - dij01(i, j);
    }
    Rcpp::checkUserInterrupt();
  }
  return;
}

NumericMatrix calc_Q(NumericMatrix rij, NumericMatrix dij11,
		     NumericMatrix dij10, NumericMatrix dij01,
		     NumericMatrix dij00) {
  int I = rij.nrow();
  int J = rij.ncol();
  NumericMatrix Q(I, J);
  int i;
  int j;

  for (i=0; i<I; i++) {
    Q(i, 0) = 1;
  }
  for (j=0; j<J; j++) {
    Q(0, j) = 1;
  }
  for (i=1; i<I; i++) {
    for (j=1; j<J; j++) {
      if (rij(i, 0)==0 || rij(0, j)==0 || rij(i, j)==0 ||
	  dij10(i, 0)==rij(i, 0) || dij01(0, j)==rij(0, j)) {
	Q(i, j) = 0;
      }
      else {
	Q(i, j) = Q(i, j-1)+Q(i-1, j)-Q(i-1, j-1)*
	  (1-(((dij11(i, j)/rij(i, j))-((dij11(i, j)+dij10(i, j))/rij(i, j))*
	       (dij01(0, j)/rij(0, j))-((dij11(i, j)+dij01(i, j))/rij(i, j))*
	       (dij10(i, 0)/rij(i, 0))+(dij10(i, 0)/rij(i, 0))*(dij01(0, j)/
								rij(0, j)))/
	      ((1-(dij10(i, 0)/rij(i, 0)))*(1-(dij01(0, j)/rij(0, j))))));
      }
    }
    Rcpp::checkUserInterrupt();
  }
  return(Q);
}

// [[Rcpp::export]]
arma::cube jointHazLambda(NumericVector Y1, NumericVector Y2, NumericVector T1,
			  NumericVector T2, NumericVector Delta1,
			  NumericVector Delta2, int method) {
  int nT1 = T1.size();
  int nT2 = T2.size();
  NumericMatrix rij(nT1+1, nT2+1);
  NumericMatrix dij11(nT1+1, nT2+1);
  NumericMatrix dij10(nT1+1, nT2+1);
  NumericMatrix dij01(nT1+1, nT2+1);
  NumericMatrix dij00(nT1+1, nT2+1);
  NumericVector T1new(nT1+1);
  NumericVector T2new(nT2+1);
  NumericMatrix Fhat(nT1+1, nT2+1);
  NumericMatrix Lambda10(nT1+1, nT2+1);
  NumericMatrix Lambda01(nT1+1, nT2+1);
  NumericMatrix Lambda11(nT1+1, nT2+1);
  NumericMatrix Q(nT1+1, nT2+1);
  arma::cube FLambda(nT1+1, nT2+1, 4);
  int i;
  int j;

  std::fill(rij.begin(), rij.end(), 0);
  std::fill(dij11.begin(), dij11.end(), 0);
  std::fill(dij10.begin(), dij10.end(), 0);
  std::fill(dij01.begin(), dij01.end(), 0);
  std::fill(dij00.begin(), dij00.end(), 0);
  std::fill(Fhat.begin(), Fhat.end(), 0);
  std::fill(Lambda10.begin(), Lambda10.end(), 0);
  std::fill(Lambda01.begin(), Lambda01.end(), 0);
  std::fill(Lambda11.begin(), Lambda11.end(), 0);

  T1new[0] = 0;
  for (i=0; i<nT1; i++) {
    T1new[i+1] = T1[i];
  }
  T2new[0] = 0;
  for (j=0; j<nT2; j++) {
    T2new[j+1] = T2[j];
  }
  calc_dij(Y1, Y2, T1new, T2new, Delta1, Delta2, dij00, dij01, dij10, dij11,
	   rij);

  Fhat(0, 0) = 1;
  for (i=0; i<=nT1; i++) {
    Lambda01(i, 0) = 0;
    Lambda11(i, 0) = 0;
    if (i>0) {
      Lambda10(i, 0) = dij10(i, 0) / rij(i, 0);
      Fhat(i, 0) = Fhat(i-1, 0)*(1 - Lambda10(i, 0));
    }
  }
  for (j=0; j<=nT2; j++) {
    Lambda10(0, j) = 0;
    Lambda11(0, j) = 0;
    if (j>0) {
      Lambda01(0, j) = dij01(0, j) / rij(0, j);
      Fhat(0, j) = Fhat(0, j-1)*(1 - Lambda01(0, j));
    }
  }

  if (method==1) {
    for (i=0; i<nT1; i++) {
      for (j=0; j<nT2; j++) {
	if (rij(i+1, j+1)>0) {
	  if (dij00(i+1, j+1)>0) {
	    Fhat(i+1, j+1) = (Fhat(i+1, j)*Fhat(i, j+1)/Fhat(i, j))*
	      dij00(i+1, j+1)*rij(i+1, j+1)/((dij00(i+1, j+1)+dij10(i+1, j+1))*
					     (dij00(i+1, j+1)+dij01(i+1, j+1)));
	  }
	  else {
	    Fhat(i+1, j+1) = 0;
	  }
	  Lambda10(i+1, j+1) = (Fhat(i, j+1)-Fhat(i+1, j+1))/Fhat(i, j+1);
	  Lambda01(i+1, j+1) = (Fhat(i+1, j)-Fhat(i+1, j+1))/Fhat(i+1, j);
	  Lambda11(i+1, j+1) = (Fhat(i, j)-Fhat(i+1, j)-Fhat(i, j+1)+
				Fhat(i+1, j+1))/Fhat(i, j);
	}
      }
      Rcpp::checkUserInterrupt();
    }
  }
  else if (method==2) {
    for (i=0; i<nT1; i++) {
      for (j=0; j<nT2; j++) {
	if (rij(i+1, j+1)>0) {
	  if (dij00(i+1, j+1)==0) {
	    Lambda11(i+1, j+1) = 0;
	  }
	  else {
	    Lambda11(i+1, j+1) = dij11(i+1, j+1)/rij(i+1, j+1);
	  }
	  Fhat(i+1, j+1) = Fhat(i+1, j)+Fhat(i, j+1)-Fhat(i, j)*
	    (1-Lambda11(i+1, j+1));
	  Lambda10(i+1, j+1) = (Fhat(i, j+1)-Fhat(i+1, j+1))/Fhat(i, j+1);
	  Lambda01(i+1, j+1) = (Fhat(i+1, j)-Fhat(i+1, j+1))/Fhat(i+1, j);
	}
      }
      Rcpp::checkUserInterrupt();
    }
  }
  if (method==3) {
    Q = calc_Q(rij, dij11, dij10, dij01, dij00);
    for (i=0; i<nT1; i++) {
      for (j=0; j<nT2; j++) {
	Fhat(i+1, j+1) = Fhat(i+1, 0)*Fhat(0, j+1)*Q(i+1, j+1);
	Lambda10(i+1, j+1) = (Fhat(i, j+1)-Fhat(i+1, j+1))/Fhat(i, j+1);
	Lambda01(i+1, j+1) = (Fhat(i+1, j)-Fhat(i+1, j+1))/Fhat(i+1, j);
	Lambda11(i+1, j+1) = (Fhat(i, j)-Fhat(i+1, j)-Fhat(i, j+1)+
			      Fhat(i+1, j+1))/Fhat(i, j);
      }
      Rcpp::checkUserInterrupt();
    }
  }

  for (i=0; i<=nT1; i++) {
    for (j=0; j<=nT2; j++) {
      FLambda(i, j, 0) = Fhat(i, j);
      FLambda(i, j, 1) = Lambda10(i, j);
      FLambda(i, j, 2) = Lambda01(i, j);
      FLambda(i, j, 3) = Lambda11(i, j);
    }
    Rcpp::checkUserInterrupt();
  }
  return(FLambda);
}

// [[Rcpp::export]]
arma::cube calcTemp2(NumericMatrix Fhat, NumericMatrix Lambda10,
		     NumericMatrix Lambda01, NumericMatrix Lambda11) {
  int I = Fhat.nrow();
  int J = Fhat.ncol();
  arma::cube temp(I, J, 4);
  int i;
  int j;

  std::fill(temp.begin(), temp.end(), 0);

  for(i=0; i<(I-1); i++) {
    for (j=0; j<(J-1); j++) {
      temp(i+1, j+1, 0) = temp(i+1, j, 0)+temp(i, j+1, 0)-temp(i, j, 0)+
	Lambda11(i+1, j+1)*Fhat(i, j);
      temp(i+1, j+1, 1) = temp(i+1, j, 1)+temp(i, j+1, 1)-temp(i, j, 1)+
	Lambda10(i+1, j)*Lambda01(i, j+1)*Fhat(i, j);
      temp(i+1, j+1, 2) = temp(i+1, j, 2)+temp(i, j+1, 2)-temp(i, j, 2)+
	(Lambda11(i+1, j+1)-Lambda10(i+1, j)*Lambda01(i, j+1))*
	Fhat(i, j)*Fhat(i, j);
      temp(i+1, j+1, 3) = temp(i+1, j, 3)+temp(i, j+1, 3)-temp(i, j, 3)+
	(Lambda11(i+1, j+1)+Lambda10(i+1, j)*Lambda01(i, j+1))*
	Fhat(i, j)*Fhat(i, j);
    }
    Rcpp::checkUserInterrupt();
  }
  return(temp);
}

// [[Rcpp::export]]
NumericMatrix jointHaz(NumericVector Y1, NumericVector Y2, NumericVector T1,
		       NumericVector T2, NumericVector Delta1,
		       NumericVector Delta2) {
  int nY = Y1.size();
  int nT1 = T1.size();
  int nT2 = T2.size();
  NumericMatrix Fhat(nT1+1, nT2+1);
  int rij;
  int dij11;
  int dij10;
  int dij01;
  int dij00;
  int i;
  int j;
  int k;

  NumericVector Fhat1 = margHaz(Y1, T1, Delta1);
  NumericVector Fhat2 = margHaz(Y2, T2, Delta2);

  for (i=0; i<=nT1; i++) {
    Fhat(i,0) = Fhat1[i];
  }
  for (j=0; j<=nT2; j++) {
    Fhat(0,j) = Fhat2[j];
  }

  for (i=0; i<nT1; i++) {
    for (j=0; j<nT2; j++) {
      if (Fhat(i,j)==0 || Fhat(i+1,j)==0 || Fhat(i,j+1)==0) {
	Fhat(i+1,j+1) = 0;
      }
      else {
	rij = 0;
	dij11 = 0;
	dij10 = 0;
	dij01 = 0;
	for (k=0; k<nY; k++) {
	  if (Y1[k]>=T1[i] && Y2[k]>=T2[j]) {
	    rij++;
	  }
	  if (Y1[k]==T1[i] && Y2[k]==T2[j] && Delta1[k]==1 && Delta2[k]==1) {
	    dij11++;
	  }
	  if (Y1[k]==T1[i] && Delta1[k]==1 &&
	      (Y2[k]>T2[j] ||
	       (Y2[k]==T2[j] && Delta2[k]==0))) {
	    dij10++;
	  }
	  if (Y2[k]==T2[j] && Delta2[k]==1 &&
	      (Y1[k]>T1[i] ||
	       (Y1[k]==T1[i] && Delta1[k]==0))) {
	    dij01++;
	  }
	}
	dij00 = rij - dij11 - dij10 - dij01;
	if (rij>0 && dij00>0) {
	  Fhat(i+1,j+1) = (Fhat(i+1,j)*Fhat(i,j+1)/Fhat(i,j)) *
	    (double(dij00)*double(rij))/((double(dij00)+double(dij10))*
					 (double(dij00)+double(dij01)));
	}
	else {
	  Fhat(i+1,j+1) = 0;
	}
      }
    }
    Rcpp::checkUserInterrupt();
  }
  return(Fhat);
}

// [[Rcpp::export]]
arma::cube triHaz(NumericVector Y1, NumericVector Y2, NumericVector Y3,
		  NumericVector T1, NumericVector T2, NumericVector T3,
		  NumericVector Delta1, NumericVector Delta2,
		  NumericVector Delta3) {
  int nY = Y1.size();
  int nT1 = T1.size();
  int nT2 = T2.size();
  int nT3 = T3.size();
  arma::cube Fhat(nT1+1, nT2+1, nT3+1);
  int rijk;
  int Ftilda111;
  int Ftilda110;
  int Ftilda101;
  int Ftilda011;
  int Ftilda100;
  int Ftilda010;
  int Ftilda001;
  double ftilda111;
  double ftilda110;
  double ftilda101;
  double ftilda011;
  double ftilda100;
  double ftilda010;
  double ftilda001;
  int i;
  int j;
  int k;
  int l;

  NumericMatrix Fhat12 = jointHaz(Y1, Y2, T1, T2, Delta1, Delta2);
  NumericMatrix Fhat13 = jointHaz(Y1, Y3, T1, T3, Delta1, Delta3);
  NumericMatrix Fhat23 = jointHaz(Y2, Y3, T2, T3, Delta2, Delta3);

  for (i=0; i<=nT1; i++) {
    for (j=0; j<=nT2; j++) {
      Fhat(i,j,0) = Fhat12(i,j);
    }
    Rcpp::checkUserInterrupt();
  }
  for (i=0; i<=nT1; i++) {
    for (k=0; k<=nT3; k++) {
      Fhat(i,0,k) = Fhat13(i,k);
    }
    Rcpp::checkUserInterrupt();
  }
  for (j=0; j<=nT2; j++) {
    for (k=0; k<=nT3; k++) {
      Fhat(0,j,k) = Fhat23(j,k);
    }
    Rcpp::checkUserInterrupt();
  }

  for (i=0; i<nT1; i++) {
    for (j=0; j<nT2; j++) {
      for (k=0; k<nT3; k++) {
	rijk = 0;
	for (l=0; l<nY; l++) {
	  if (Y1[l]>=T1[i] && Y2[l]>=T2[j] && Y3[l]>=T3[k]) {
	    rijk++;
	  }
	}
	if (rijk>0) {
	  Ftilda111 = 0;
	  Ftilda110 = 0;
	  Ftilda101 = 0;
	  Ftilda011 = 0;
	  Ftilda100 = 0;
	  Ftilda010 = 0;
	  Ftilda001 = 0;
	  for (l=0; l<nY; l++) {
	    if (Y1[l]>T1[i] && Y2[l]>T2[j] && Y3[l]>T3[k]) {
	      Ftilda111++;
	    }
	    if (Y1[l]>T1[i] && Y2[l]>T2[j] && Y3[l]>=T3[k]) {
	      Ftilda110++;
	    }
	    if (Y1[l]>T1[i] && Y2[l]>=T2[j] && Y3[l]>T3[k]) {
	      Ftilda101++;
	    }
	    if (Y1[l]>=T1[i] && Y2[l]>T2[j] && Y3[l]>T3[k]) {
	      Ftilda011++;
	    }
	    if (Y1[l]>T1[i] && Y2[l]>=T2[j] && Y3[l]>=T3[k]) {
	      Ftilda100++;
	    }
	    if (Y1[l]>=T1[i] && Y2[l]>T2[j] && Y3[l]>=T3[k]) {
	      Ftilda010++;
	    }
	    if (Y1[l]>=T1[i] && Y2[l]>=T2[j] && Y3[l]>T3[k]) {
	      Ftilda001++;
	    }
	  }
	  ftilda111 = double(Ftilda111) / double(rijk);
	  ftilda110 = double(Ftilda110) / double(rijk);
	  ftilda101 = double(Ftilda101) / double(rijk);
	  ftilda011 = double(Ftilda011) / double(rijk);
	  ftilda100 = double(Ftilda100) / double(rijk);
	  ftilda010 = double(Ftilda010) / double(rijk);
	  ftilda001 = double(Ftilda001) / double(rijk);
	  if (ftilda110>0 && ftilda101>0 && ftilda011>0) {
	    Fhat(i+1,j+1,k+1) = Fhat(i+1,j+1,k)*Fhat(i+1,j,k+1)*
	      Fhat(i,j+1,k+1)*Fhat(i,j,k)/(Fhat(i+1,j,k)*Fhat(i,j+1,k)*
					   Fhat(i,j,k+1))*
	      ftilda111*ftilda100*ftilda010*ftilda001/
	      (ftilda110*ftilda101*ftilda011);
	  }
	  else {
	    Fhat(i+1,j+1,k+1) = 0;
	  }
	}
	else {
	  Fhat(i+1,j+1,k+1) = 0;
	}
      }
      Rcpp::checkUserInterrupt();
    }
  }
  return(Fhat);
}

// [[Rcpp::export]]
arma::cube calcTemp3(arma::cube Fhat) {
  int I = Fhat.n_rows-1;
  int J = Fhat.n_cols-1;
  int K = Fhat.n_slices-1;
  arma::cube temp(I, J, K);
  int i;
  int j;
  int k;

  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      for (k=0; k<K; k++) {
	if (Fhat(i, j, k+1)>0 && Fhat(i, j+1, k)>0 && Fhat(i+1, j, k)>0) {
	  temp(i, j, k) = Fhat(i+1, j+1, k)*Fhat(i+1, j, k+1)*
	    Fhat(i, j+1, k+1)/(Fhat(i, j, k+1)*Fhat(i, j+1, k)*
			       Fhat(i+1, j, k));
	}
	else {
	  temp(i, j, k) = 0;
	}
      }
      Rcpp::checkUserInterrupt();
    }
  }
  return(temp);
}

NumericVector Quantile(NumericVector x, NumericVector probs) {
  // implementation of type 7
  const size_t n=x.size(), np=probs.size();
  if (n==0) return x;
  if (np==0) return probs;
  NumericVector index = (n-1.)*probs, y=x.sort(), x_hi(np), qs(np);
  NumericVector lo = Rcpp::floor(index), hi = Rcpp::ceiling(index);

  for (size_t i=0; i<np; ++i) { 
    qs[i] = y[lo[i]];
    x_hi[i] = y[hi[i]];
    if ((index[i]>lo[i]) && (x_hi[i] != qs[i])) {
      double h;
      h = index[i]-lo[i];
      qs[i] = (1.-h)*qs[i] + h*x_hi[i];
    }
  }

  return qs;
}

// [[Rcpp::export]]
NumericMatrix getBootCI(NumericMatrix T) {
  int J = T.ncol();
  NumericVector quants = NumericVector::create(0.025, 0.975);
  NumericVector curci(2);
  NumericMatrix outci(2, J); 
  int j;

  for (j=0; j<J; j++) {
    curci = Quantile(T(_,j), quants);
    outci(0,j) = curci[0];
    outci(1,j) = curci[1];
  }
  return(outci);
}
