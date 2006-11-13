#include <math.h>
#include "utils.h"
#include "poiss_prob.h"

//
// Initializes all the data members
//
template<class p_t>
poiss_prob_t<p_t>::poiss_prob_t(int given_N, int given_K, double* pu) {

  N = given_N; K = given_K; double* log_pu = new double[K];
  for(int i = 0; i < K; i++) log_pu[i] = log(pu[i]);    
  
  log_factor = new double[N+1]; log_factor[0] = 0;
  for(int i = 1; i < N+1; i++) log_factor[i] = log_factor[i-1] + log((double) i);
  
  log_p = new double*[K]; p = new p_t*[K];
  double log_N = log((double) N); 
  for(int i = 0; i < K; i++) {
    
    log_p[i] = new double[N+1]; p[i] = new p_t[N+1];
    for(int j = 0; j < N+1; j++) 
      log_p[i][j] = -N*pu[i] + (log_N + log_pu[i])*j 
	- log_factor[j];		       
  }
  
  delete[](log_pu);
}

//
// Applies the theta1 shift to log_p
//
template<class p_t>
void poiss_prob_t<p_t>::shift1(double theta, int** Ip, double step) {

  theta1 = theta;
  for(int i = 0; i < K; i++) 
    for(int j = 0; j < N+1; j++)
      log_p[i][j] += theta1*step*Ip[i][j];
}

//
// Normalizes log_p after the theta1 shift
//
template<class p_t>
void poiss_prob_t<p_t>::norm1(double norm_factor) {

  for(int k = 0; k < K; k++)
    for(int n = 0; n < N+1; n++)
      log_p[k][n] -= norm_factor;
} 

//
// Applies the theta2 shift to log_p
//
template<class p_t>
void poiss_prob_t<p_t>::shift2(double theta) {
  
  theta2 = theta;
  for(int i = 0; i < K; i++)
    for(int j = 0; j < N+1; j++)
      log_p[i][j] -= theta2*j;
}

//
// Normalizes log_p[k] after the theta2 shift
//
template<class p_t>
void poiss_prob_t<p_t>::norm2(int k, double norm_factor) {

  for(int n = 0; n < N+1; n++)
    log_p[k][n] -= norm_factor;
} 

template<class p_t>
poiss_prob_t<p_t>::~poiss_prob_t() {
  
  delete[](log_factor);
  
  for(int i = 0; i < K; i++) {
    
    delete[](log_p[i]);
    delete[](p[i]);
  }
  
  delete[](log_p); delete[](p);
}
