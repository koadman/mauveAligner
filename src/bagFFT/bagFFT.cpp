//#define NR 1

#include "poiss_prob.h"
#include "llr_score.h"
#include "bagFFT-algorithm.h"
#include "utils.h"
#include "bagFFT.h"

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
bagFFT_return_value bagFFT(int N, int K, double s, double* pu, int Q) {

#if NR
  Q = RND_UP_TO_POWER_2(Q);
#endif

  //
  // The value to be returned.
  //
  bagFFT_return_value rv;

  //
  // Running the bagFFT algorithm to get the lattice pmf.
  //
  rv.log_pmf = bagFFTalgo_pmf<llr_score_t, poiss_prob_t<double> >
    (N, K, s, pu, Q, rv.step, rv.theta1, rv.log_mgf, rv.log_error_bound);

  //
  // Computing the p-value bounds and the error bounds.
  //
  rv.Q = Q;
  compute_theta_sum(rv.log_pmf, Q, s/rv.step, K, rv.theta1, 1, rv.log_min_pvalue, rv.log_max_pvalue);
  compute_theta_sum(rv.log_pmf, Q, s/rv.step, K, rv.theta1, 0, rv.log_min_error_bound, rv.log_max_error_bound);

  //
  // We estimate the p-value to be returned to be the geometic mean of the bounds.
  //
  rv.log_pvalue = (rv.log_max_pvalue + rv.log_min_pvalue)/2;

  return rv;
}
