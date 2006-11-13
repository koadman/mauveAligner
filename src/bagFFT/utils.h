#ifndef _UTILS_H_
#define _UTILS_H_

#define LOGZERO -1e100 
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define DEBUG 0
#define PRECISION 20

#include <iostream>
#include <iomanip>

using namespace std;

#include "extended_exponent.h"

//
// Double machine precision
//
#define EPS 2.220446049250313e-16

//
// Rounds to nearest integer
//
#define NINT(x) ((int) ((x)+((x) > 0 ? 0.5 : -0.5)))

#define LOG2(X) (log((double) X)/log(2.0))
#define RND_UP_TO_POWER_2(Y) (01 << (int) ceil(LOG2(Y)))

const double pi = 3.14159265358979;

//
// Computes log(exp(log_a) + exp(log_b))
//
double log_sum(double log_a, double log_b);

//
// Computes log_sum for the entire array of size n
//
double log_array_sum(double* log_x, int n);

//
// Returns data[min_index] where data[min_index] = min(data[0..(size-1)])
//  
double min_array(double* data, int size, int& min_index);

//
// Sets out[0..(row_size-1)][0..(column_size-1)] to
// exp(data[0..(row_size-1)][0..(column_size-1)]).
//
void exp_2Darray(double** data, int row_size, int column_size, double** out);
void exp_2Darray(double** data, int row_size, int column_size, ee_t** out);

inline int sign(double x) { return (x >= 0 ? 1 : -1); } 
//inline double log(int x) { return log((double) x); }

//
// Returns a*b mod c (computed carefully to avoid overflows).
//
inline unsigned long long mod(int a, int b, int c) 
{ return (((unsigned long long) a)*b)%c; }

//
// Computes the L1 norm of the vector x of size n
// 
double one_norm(double* x, int n);

//
// Computes the L2 norm of the vector x of size n
// 
double two_norm(double* x, int n);

//
// Reads the parameter list for a program that computes
// the p-value for the multinomial llr statistic
//
void read_llr_parameters(int argc, char** argv, int& N, int& K, double& s, int& Q, double*& pu);

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
		       double& log_min_sum, double& log_max_sum ); 

#endif 
 
