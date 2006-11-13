#include <math.h>
#include "utils.h"

//
// Returns data[min_index] where data[min_index] = min(data[0..(size-1)])
//  
double min_array(double* data, int size, int& min_index) {

  double min_value = data[0]; min_index = 0;
  for(int i = 1; i < size; i++)
    if(min_value > data[i]) {
      min_value = data[i]; min_index = i;
    }

  return min_value;
}

//
// Sets out[0..(row_size-1)][0..(column_size-1)] to
// exp(data[0..(row_size-1)][0..(column_size-1)])
//
void exp_2Darray(double** data, int row_size, int column_size, double** out) {

  for(int i = 0; i < row_size; i++)
    for(int j = 0; j < column_size; j++)
      out[i][j] = exp(data[i][j]);
}

//
// Sets out[0..(row_size-1)][0..(column_size-1)] to
// exp(data[0..(row_size-1)][0..(column_size-1)])
//
void exp_2Darray(double** data, int row_size, int column_size, ee_t** out) {

  for(int i = 0; i < row_size; i++)
    for(int j = 0; j < column_size; j++)
      out[i][j].log_set(data[i][j]);
}

//
// Computes log(exp(log_a) + exp(log_b))
//
double log_sum(double log_a, double log_b) {
  
  return ((log_a > log_b) ? 
	  log_a + log(1+exp(log_b-log_a)) :
	  log_b + log(1+exp(log_a-log_b)));
}

//
// Computes log_sum for the entire array of size n
//    
double log_array_sum(double* log_x, int n) {

  double array_sum = log_x[0];
  for(int i = 1; i < n; i++)
    array_sum = log_sum(array_sum, log_x[i]);

  return array_sum;
}

//
// Computes the L1 norm of the vector x of size n
// 
double one_norm(double* x, int n) {

  double sum = 0;
  for(int i = 0; i < n; i++)
    sum += fabs(x[i]);

  return sum;
}

//
// Computes the L2 norm of the vector x of size n
// 
double two_norm(double* x, int n) {

  double sum = 0;
  for(int i = 0; i < n; i++)
    sum += x[i]*x[i];

  return sqrt(sum);
}

//
// Reads the parameter list for a program that computes
// the p-value for the multinomial llr statistic
//
void read_llr_parameters(int argc, char** argv, int& N, int& K, double& s, int& Q, double*& pu) {

  if(argc < 5) {
    cout << "Computes the p-value for the multinomial llr statistic" << endl << endl
	 << "Usage: <program name> N K s Q pu" << endl
	 << "N = Number of objects" << endl
	 << "K = Number of bins" << endl
	 << "s = LLR score" << endl
	 << "Q = Lattice size" << endl
	 << "pu = List of background probabilities (of size K)" << endl << endl
	 << "Example: <program name> 50 4 50 16384 0.25 0.25 0.25 0.25" << endl;
    exit(1);
  }

  int argno = 1;
  N = atoi(argv[argno++]); K = atoi(argv[argno++]); 
  s = atof(argv[argno++]); Q = atoi(argv[argno++]); 
  pu = new double[K];

  cout << "Parameters: N = " << N << " K = " << K << " s = " << s << " Q = " << Q << " pu =";
  for(int i = 0; i < K; i++) { pu[i] = atof(argv[argno++]); cout << " " << pu[i]; }
  cout << endl;
}

void compute_pvalues(double* log_pmf, double step, int Q, double s, int K, 
		     double& log_min_pval, double& log_max_pval) {


}

//
// Computes the log-sum of the entries of log_v, using theta to decide the
// tail that should be summed.
//
// size = size of the log_v array.
// start = location to start summing from.
// uncertainity = the margin used for computing log_min_sum and log_max_sum.
// invert = if set to 1, invert the sum when theta < 0.
//
void compute_theta_sum(double* log_v, int size, double start, int uncertainity, double theta, int invert, 
		       double& log_min_sum, double& log_max_sum ) { 
  
  int i; log_max_sum = LOGZERO; log_min_sum = LOGZERO;  
  for(i = (int) MAX(0, start - uncertainity/2.0); i < (int) MIN(start + uncertainity/2.0 + 1, size); i++)
    log_max_sum = log_sum(log_max_sum, log_v[i]);

  if(theta > 0)
    for(; i < size; i++)
      log_min_sum = log_sum(log_min_sum, log_v[i]);
  else
    for(i = 0; i < (int) MAX(0, start - uncertainity/2.0); i++)
      log_min_sum = log_sum(log_min_sum, log_v[i]);

  log_max_sum = log_sum(log_max_sum, log_min_sum);

  if(theta < 0) {
    double log_temp = (invert ? log(1-exp(log_max_sum)) : log_max_sum);
    log_max_sum = (invert ? log(1-exp(log_min_sum)) : log_min_sum);
    log_min_sum = log_temp;
  }
}
