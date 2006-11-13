#include "minimize.h"
#include "mnbrak.h"
#include "brent.h"

//
// Numerically solves for the minimum of the function f.
// value is set to the argmin of f and f(value) is returned.
//
double minimize(double (*f)(double), double& value) {

  //
  // The search starts in the range [-1, 3] which seems
  // to be a reasonable choice for bagFFT.
  //
  double lower_bnd = -1; double middle = 0.5; double upper_bnd = 3.0;

  //
  // Desired level of precision.
  //
  double tolerance = 1e-6;

  //
  // The numerical recipe procedures mnbrak and brent are used
  // to find the minimum.
  //
  double fa, fb, fc;
  mnbrak(&lower_bnd, &middle, &upper_bnd, &fa, &fb, &fc, f);
  return brent(lower_bnd, middle, upper_bnd, f, tolerance, &value);
} 
