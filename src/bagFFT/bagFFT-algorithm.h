#include <math.h>
#include "utils.h"
#include "theta_fns.h"
#include "minimize.h"
#include "complex.h"
#include "convolution.h"

//
// Function: Computes and returns the DFT of the shifted lattice pmf (normalized complex array of size Q)
//           based on the bagFFT algorithm.
//
// Parameters: N = Number of objects 
//             K = Number of bins 
//             Q = Lattice size 
//             score = Object encapsulating the partial score contributions
//             p[k][x] = probability contribution for x objects of the kth type
//             N2 = Array size for convolutions.
//             log_theta2_offset = Offset to compensate for normalization after shifting by theta2
//
template<class score_t>
complex* compute_phi_theta(int N, int K, int Q, score_t& score, double** p, int N2, double log_theta2_offset) {

  //
  // Omega
  //
  double omg = 2*pi/Q;
  
  //
  // exp_i_omg[x] = exp(i*omg*x)
  //
  complex *exp_i_omg = new complex[Q];
  //boost::shared_array<complex> exp_i_omg (new complex[Q]);
  for(int x = 0; x < Q; x++)
    exp_i_omg[x].set(cos(omg*x), sin(omg*x));

  //
  // phi_theta[l] = lth element of the DFT of the shifted lattice pmf.
  //
  complex *phi_theta = new complex[Q];
  //boost::shared_array<complex> phi_theta (new complex[Q]);
  //
  // Computing phi_theta[l].
  //
  
  for(int l = 0; l < 1+Q/2; l++) {

    //
    // Initializing curr_psi (= in LATEX notation $\Psi_{k, l, \theta, \theta_2}$ from the paper).
    //
    int k = 0; 
    //complex curr_psi[];
	
	//complex* curr_psi = new complex[N2];
	
	
	//vector<complex>::pointer complexp;
	//vector<complex> curr_psi(N2);

	//tjt: can boost help us?
	//boost::scoped_array<complex> curr_psi (new complex[N2]);
    boost::shared_ptr<complex> curr_psi (new complex[N2]);

	//vector<boost::shared_ptr<complex>> curr_psi;
	//curr_psi.push_back( new complex(N2));
	for(int x = 0; x < N2; x++)
      curr_psi.get()[x] = exp_i_omg[mod(score.get_Ip(k, x)+Q, l, Q)]*p[k][x];

	//vector<boost::shared_ptr<complex>> p_exp_k;
	//vector<boost::shared_ptr<complex>> next_psi;
    for(k = 1; k < K; k++) {

      //
      // Initializing p_exp_k (= in LATEX notation $p_{k, l, \theta, \theta_2}$ from the paper).
      //
      //complex p_exp_k[N2];    
	  boost::shared_ptr<complex> p_exp_k (new complex[N2]);
      //complex* p_exp_k = new complex[N2];
	
      for(int x = 0; x < N2; x++) 
		p_exp_k.get()[x] = exp_i_omg[mod(score.get_Ip(k, x)+Q, l, Q)]*p[k][x];

      //complex next_psi[N2];  
	  
	  boost::shared_ptr<complex> next_psi (new complex[N2]);
	  //vector<complex> next_psi(N2);
	  // complex* next_psi = new complex[N2];

#if NR

      conv_fftnr(curr_psi.get(), p_exp_k.get(), next_psi.get(), N2);
#else
	  //gets used in a strange way.. needs only curr_psi[0] ??
      //complex* curr_psi2 = new complex[3];
	  //curr_psi2[0] = curr_psi[0];
	  
      //conv_fftw((complex*)&curr_psi.at(0), p_exp_k, next_psi, N2);
	  conv_fftw(curr_psi.get(), p_exp_k.get(),next_psi.get(), N2);
#endif
      //delete[](p_exp_k);
      
	  //next_psi.reset();
      //
      // Updating curr_psi
      //
      for(int n = 0; n < N2; n++)
	  {
		if(n < N+1)  
		{
		  
		  curr_psi.get()[n] = next_psi.get()[n];
		}
		else
		  curr_psi.get()[n].set(0, 0);    

	  }
	  //next_psi = 0;
	  //delete[](next_psi);
	  //next_psi.reset();
    
    }
 
    //
    // Shifting curr_psi[N] to get phi_theta[l]
    //
    phi_theta[l] = curr_psi.get()[N].shift(log_theta2_offset);

		
    //tjt: delete is ok when NR=1..
	//delete[](curr_psi);
	
    //
    // Exploiting the symmetry for computing the DFT
    // of a real vector.
    //
    if(l > 0)
      phi_theta[Q-l] = phi_theta[l].conj();

	
	curr_psi.reset();
	
  }

  delete[](exp_i_omg);
  
  return phi_theta;
}

//
// Function: Computes and returns the entire lattice pmf (log-values
//           in an array of size Q) based on the bagFFT algorithm.
//
// Parameters: N = Number of objects 
//             K = Number of bins 
//             s = Score for which theta is optimized
//             pu = Null distribution parameters 
//             Q = Lattice size (Should be a power of 2 if NR = 1)
//             step = Step size for the lattice (set by the function)
//             theta1 = Optimal theta1 shift (set by the function)
//             log_mgf = Log of the mgf at theta1 (set by the function)
//             log_error_bound = Bound on the absolute error for the entries of 
//                              the pmf (set by the function)
//
template<class score_t, class prob_t>
double* bagFFTalgo_pmf(int N, int K, double s, double* pu, int Q, double& step, double& theta1, 
		       double& log_mgf, double*& log_error_bound) {

  //
  // Initializing the score and probability contributions.
  //
  score_t score(N, K, pu, Q);
  prob_t prob(N, K, pu);
  step = score.get_step();
  
  //
  // Minimizing theta1
  //
  theta1Fn<score_t, prob_t>::init(&score, &prob, s);
  log_mgf = minimize(theta1Fn<score_t, prob_t>::minimizationFn, theta1); 
  log_mgf += s*theta1;

  //
  // Shifting and normalizing log_p based on theta1
  //
  prob.shift1(theta1, score.get_Ip_copy(), score.get_step());
  prob.norm1(log_mgf/((double) K));

  //
  // N2 = 2*N+1
  //
  // Note that while rounding up to a power of 2 is only essential for the numerical recipes version,
  // in practice it was found that it makes the FFTW version more efficient too.
  //
  int N2 = RND_UP_TO_POWER_2(2*N+1);

  //
  // Minimizing theta2
  //
  double theta2; theta2Fn<score_t, prob_t>::init(&score, &prob, N2);
  double log_psi_error_bound = minimize(theta2Fn<score_t, prob_t>::minimizationFn, theta2);

  //
  // Shifting and normalizing log_p based on theta2
  //
  prob.shift2(theta2);
  double off_sum = 0; double** log_p = prob.get_log_p_copy();
  for(int k = 0; k < K; k++) {

    double off = log_array_sum(log_p[k], N+1); off_sum += off;
    prob.norm2(k, off);
  }
  double log_theta2_offset = theta2*N + off_sum + prob.log_const_factor();
  double log_phi_error_bound = log_psi_error_bound + prob.log_const_factor();

  //
  // Exponentiating log_p to get p
  //
  prob.set_p();
  double** p = prob.get_p_copy();

#if DEBUG 
  cout << "log(mgf): " << log_mgf << endl;
  cout << "theta1: " << theta1 << endl;
  cout << "theta2: " << theta2 << endl;
#endif

  //
  // Computing phi_theta.
  //
  complex* phi_theta = compute_phi_theta(N, K, Q, score, p, N2, log_theta2_offset);

  //
  // Applying the inverse DFT operator.
  //

#if NR
  fftnr(phi_theta, Q, -1);
#else
  fftw(phi_theta, Q, -1);
#endif

  //
  // log_pmf = log lattice pmf (return value)
  //
  double* log_pmf = new double[Q];
  //boost::scoped_array<double> log_pmf (new double[Q]);
  log_error_bound = new double[Q];
  //boost::scoped_array<double> log_error_bound (new double[Q]);

  //
  //s Setting the pmf and the bound arrays.
  //
  double log_p_theta_error_bound = log_sum(log_phi_error_bound, log(10*log((double) Q)/log(2.0)*EPS)); 
  for(int j = 0; j < Q; j++) {
    log_pmf[j] = (phi_theta[j].real() > 0 ? log(phi_theta[j].real()) : LOGZERO) - log((double) Q) - 
      theta1*step*j + log_mgf;
    log_error_bound[j] = log_p_theta_error_bound - theta1*step*j + log_mgf;
  }
  
  return log_pmf;
}
