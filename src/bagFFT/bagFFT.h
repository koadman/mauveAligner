#ifndef _BAGFFT_H_
#define _BAGFFT_H_

struct bagFFT_return_value {

  //
  // An estimate for the p-value;
  //
  double log_pvalue;

  //
  // Upper and lower bounds for the p-value and
  // the corresponding bounds on the absolute error.
  //
  double log_max_pvalue;
  double log_min_pvalue;
  double log_max_error_bound;
  double log_min_error_bound;

  //
  // step = Step size for the lattice
  // theta1 = Optimal theta1 shift
  // log_mgf = Log of the mgf at theta1
  //
  double step;
  double theta1;
  double log_mgf;

  //
  // Lattice pmf and the corresponding absolute error bounds 
  // (log-values in an array of size Q).
  //
  int Q;
  double* log_pmf;
  double* log_error_bound;
};

//
// Function: Computes and returns a bagFFT_return_value structure 
//           for the multinomial llr score based on the bagFFT algorithm.
//
// Parameters: N = Number of objects 
//             K = Number of bins 
//             s = LLR score for which the p-value is to be computed
//             pu = Bin probabilities for the null multinomial distribution 
//             Q = Lattice size
// 
bagFFT_return_value bagFFT(int N, int K, double s, double* pu, int Q);

#endif
