#ifndef _MINIMIZE_H_
#define _MINIMIZE_H_

//
// Numerically solves for the minimum of the function f.
// value is set to the argmin of f and f(value) is returned.
//
double minimize(double (*f)(double), double& value);

#endif
