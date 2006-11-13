#ifndef _POISS_PROB_H_
#define _POISS_PROB_H_

#include <math.h>

//
// Class for encapsulating the data and the methods for accessing them
// for the multinomial probability distribution (for use in the bagFFT
// algorithm) modelled as the sum of poissons. 
// If a different distribution needs to be used, then we can do so 
// by creating a class which has the same public interface as in this class.
//
// Note: Ideally, poiss_prob_t would be derived from a prob_t class, 
//       but we avoid this for efficiency reasons.
//
// p_t = Data type for storing the probability contributions.
//       The function exp_2Darray should be defined with an "out"
//       argument of this type (see utils.h).
//

template<class p_t>
class poiss_prob_t {

 private:

  //
  // N = Number of objects 
  // K = Number of bins 
  // p[k][n] = poisson distribution for the kth bin 
  //           having n objects in it = exp(-N*pu[k])*(N*pu[k])^n/n!
  // log_p[k][n] = log(p[k][n])
  // log_factor[n] = log(n!)
  //
  int N; int K; double** log_p; p_t** p; 
  double* log_factor;

  //
  // Thetas for the shifts
  //
  double theta1; double theta2;

 public:

  //
  // Initializes all the data members
  //
  poiss_prob_t(int given_N, int given_K, double* pu);

  //
  // Functions for applying the theta1 shift to log_p
  // and normalizing it.
  //
  void shift1(double theta, int** Ip, double step);
  void norm1(double norm_factor);

  //
  // Functions for applying the theta2 shift to log_p
  // and normalizing it.
  //
  void shift2(double theta);
  void norm2(int k, double norm_factor);

  void set_p() { exp_2Darray(log_p, K, N+1, p); }
  inline double get_log_p(int k, int n) { return log_p[k][n]; }
  inline double** get_log_p_copy() { return log_p; }
  inline p_t get_p(int k, int n) { return p[k][n]; }
  inline p_t** get_p_copy() { return p; }
  inline int get_N() { return N; }
  inline int get_K() { return K; }
  inline double get_theta1() { return theta1; }
  inline double get_theta2() { return theta2; }
  
  //
  // Returns log of the missing constant factor in the probability
  // i.e. log(exp(N)*N!/N^N)
  //
  inline double log_const_factor() { return N+log_factor[N]-N*log((double) N); }

  ~poiss_prob_t();

};

#include "poiss_prob.cpp"

#endif
