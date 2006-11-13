#include "utils.h"
#include "llr_score.h"

//
// Initializes all the data members
// Q = Lattice size
//
llr_score_t::llr_score_t(int given_N, int given_K, 
			 double* given_pu, int given_Q) {

  N = given_N; K = given_K; pu = given_pu; Q = given_Q;

  double* log_pu = new double[K];
  for(int i = 0; i < K; i++) log_pu[i] = log(pu[i]);    
  
  double* logn = new double[N+1]; logn[0] = LOGZERO;
  for(int i = 1; i < N+1; i++) logn[i] = log((double) i);    
  
  int imin; min_array(pu, K, imin);
  max_ent = -N*log_pu[imin]; step = max_ent/(Q-1);
  
  //
  // Ip[k][n] = round(n*log(n/(N*pu[k]))/step)
  //
  Ip = new int*[K];
  for(int i = 0; i < K; i++) {
    
      Ip[i] = new int[N+1]; Ip[i][0] = 0;
      for(int j = 1; j < N+1; j++) 
	Ip[i][j] = NINT(j*(logn[j]-logn[N]-log_pu[i])/step);
  }
  
  delete[](log_pu); delete[](logn);
}

llr_score_t::~llr_score_t() {

  for(int i = 0; i < K; i++) delete[](Ip[i]);
  delete[](Ip);
}
