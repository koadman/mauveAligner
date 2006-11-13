#ifndef _THETA_OBJS_H_
#define _THETA_OBJS_H_

#include "extended_exponent.h"
#include "utils.h"
#include "convolution.h"

//
// Namespace for the variables and functions
// to use to compute theta1 in the bagFFT algorithm
//
template<class score_t, class prob_t>
class theta1Fn {
  
 private:

  static score_t* score; static prob_t* prob;

 public:

  static double s;

  //
  // To initialize the namespace
  //
  static void init(score_t* given_score, prob_t* given_prob, double given_s) 
    { score = given_score; prob = given_prob; s = given_s; }

  //
  // Returns -theta*s + log(M(theta))
  //
  static double minimizationFn(double theta) 
{

    int N = prob->get_N(); 

    double step = score->get_step();
    //ee_t op1[N+1]; 
    //ee_t op2[N+1];
    //ee_t* op1 = new ee_t[N+1];
    //ee_t* op2 = new ee_t[N+1];
    boost::scoped_array<ee_t> op1 (new ee_t[N+1]);
	boost::scoped_array<ee_t> op2 (new ee_t[N+1]);
    //
    // log(M(theta)) is computed by
    // convolving the shifted p[k]'s in extended exponent
    // arithmetic 
    //
    for(int i = 0; i < N+1; i++)
      op1[i].log_set(prob->get_log_p(0, i) + 
		     theta*step*score->get_Ip(0, i));
    
    for(int k = 1; k < prob->get_K(); k++) 
    {
      for(int i = 0; i < N+1; i++)      
	     op2[i].log_set(prob->get_log_p(k, i) + theta*step*score->get_Ip(k, i));
      
      conv(op2.get(), op1.get(), op1.get(), N+1, N+1);
    }
	double op1value = op1[N].log_get();
    //delete[](op1);
	//delete[](op2);
	op1.reset();
	op2.reset();
    return -s*theta + op1value + prob->log_const_factor();
  }

};

template<class score_t, class prob_t> score_t* theta1Fn<score_t, prob_t>::score; 
template<class score_t, class prob_t> prob_t* theta1Fn<score_t, prob_t>::prob;
template<class score_t, class prob_t> double theta1Fn<score_t, prob_t>::s;  

//
// Namespace for the variables and functions
// to use to compute theta2 in the bagFFT algorithm
//
template<class score_t, class prob_t>
class theta2Fn {

 private:

  static score_t* score; static prob_t* prob;
  
  //
  // FFT array size
  //
  static int N2;

 public:

  //
  // To initialize the namespace
  //  
  static void init(score_t* given_score, prob_t* given_prob, int given_N2) 
    { score = given_score; prob = given_prob; N2 = given_N2; }    

  //
  // Sets out[i] to exp(log_in[i]-theta*i)/C(theta) where
  // C(theta) is such that out sums to 1. Returns log(C(theta)).
  //
  static double shift_and_normalize(double* log_in, int size, double theta, double* out) {

    for(int i = 0; i < size; i++)
      out[i] = log_in[i] - theta*i; 
    
    double log_sum = log_array_sum(out, size);
    
    for(int i = 0; i < size; i++)
      out[i] = exp(out[i] - log_sum);
    
    return log_sum;
  }

  //
  // Returns the error bound for doing convolutions with
  // theta as described in the bagFFT paper
  //
  static double minimizationFn(double theta) 
{

    int N = prob->get_N();
    //double op1[N2]; 
    //double op2[N2];
    //double* op1 = new double[N2];
    //double* op2 = new double[N2];
    boost::scoped_array<double> op1 (new double[N2]);
	boost::scoped_array<double> op2 (new double[N2]);
    
	
	for(int i = 0; i < N2; i++)
      op1[i] = op2[i] = 0;

    double** log_p = prob->get_log_p_copy();
    double step = score->get_step();

    double total_log_sum = shift_and_normalize(log_p[0], N+1, theta, op1.get());
    //double op3[N+1];
    //double* op3 = new double[N+1];
    boost::scoped_array<double> op3 (new double[N+1]);
    for(int i = 0; i < N+1; i++)
      op3[i] = op1[i]*(1 + (i > 0 ? log((double) i): i) + 
      	       fabs((prob->get_theta1()-1)*score->get_Ip(0, i)*step) + 
      	       fabs((theta-1)*i));
    double err = two_norm(op3.get(), N+1);
  
    int dft_error_const = 10;
    for(int k = 1; k < prob->get_K(); k++) {
    
      total_log_sum += shift_and_normalize(log_p[k], N+1, theta, op2.get());
      for(int i = 0; i < N+1; i++)
	op3[i] = op2[i]*(1 + (i > 0 ? log((double) i): i) + 
			 fabs((prob->get_theta1()-1)*score->get_Ip(k, i)*step) + 
			 fabs((theta-1)*i));
      double op3_two_norm = two_norm(op3.get(), N+1);
      err += two_norm(op1.get(), N+1)*(2*dft_error_const*log((double) N2)/log(2.0) + 5) +
	one_norm(op1.get(), N+1)*(two_norm(op2.get(), N+1)*dft_error_const*log((double) N2)/log(2.0) + op3_two_norm);

#if NR
      real_conv_fftnr(op2.get(), op1.get(), op1.get(), N2);
#else
      real_conv_fftw(op2.get(), op1.get(), op1.get(), N2);
#endif

      for(int n = N+1; n < N2; n++) op1[n] = 0; 
    }
    
	//delete[](op1);
	//delete[](op2);
	//delete[](op3);
	op1.reset();
	op2.reset();
	op3.reset();

    return log(err*EPS) + theta*N + total_log_sum;
  }
};

template<class score_t, class prob_t> score_t* theta2Fn<score_t, prob_t>::score;
template<class score_t, class prob_t> prob_t* theta2Fn<score_t, prob_t>::prob;  
template<class score_t, class prob_t> int theta2Fn<score_t, prob_t>::N2;

#endif
