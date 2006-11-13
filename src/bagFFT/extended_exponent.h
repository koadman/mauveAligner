#ifndef _EXTENDED_EXPONENT_H_
#define _EXTENDED_EXPONENT_H_

class ee_t;

#include <math.h>
#include "utils.h"

static const int max_shift = 60;
#define EXP1 2.71828182845905
#define ZEROEXP -10000

//
// shift[i] = exp(i) for i in [0..max_shift-1]
//
inline double* init() {

  double* shift = new double[max_shift];
  for(int i = 0; i < max_shift; i++)
    shift[i] = exp((double) i);
  return shift;
}
static double* shift = init();

//
// Class for representing values with very small
// or very large exponents
//
class ee_t {

 private:

  //
  // The value represented is base*exp(exponent)
  //
  double base;
  int exponent;

  inline void normalize(ee_t& v) 
    { if(v.base >= EXP1) { v.base /= EXP1; v.exponent++; } }

  inline ee_t(double given_base, int given_exponent) { base = given_base; exponent = given_exponent; }

 public:
  
  //
  // Functions for setting and getting log of the represented value.
  //
  inline double log_get() 
    { return (base == 0 ? ZEROEXP : exponent+log(base)); }

  inline void log_set(double log_v) {

    exponent = (log_v <= LOGZERO ? ZEROEXP : (int) log_v);
    base = (log_v <= LOGZERO ? 0 : exp(log_v-exponent));
  }

  inline ee_t() { base = 0; exponent = ZEROEXP; }
  inline ee_t(double log_v) { log_set(log_v); }

  inline double get() { return base*exp((double) exponent); }

  inline ee_t& operator=(const ee_t& a) 
    { base = a.base; exponent = a.exponent; return *this; }

  inline ee_t& operator=(double log_v) 
    { log_set(log_v); return *this; }

  inline int operator==(double v) {
    return (v == 0 ? (base == 0) : (get() == v));
  }

  inline int operator!=(double v) {
    return (v == 0 ? (base != 0) : (get() != v));
  }

  inline ee_t operator*(const ee_t& a) 
    { ee_t result(base*a.base, exponent+a.exponent); normalize(result); return result; }

  inline ee_t operator+(const ee_t& a) {

    //
    // Note: This special case is handled in this way
    // to optimize the runtime
    //
    if(base == 0)
      return a;

    ee_t result; 
    int diff = exponent - a.exponent;				

    //
    // Case 1: a is smaller => a.base needs to be shifted
    // 
    if(diff > 0)
      if(diff < max_shift) {
	result.exponent = exponent;
	result.base = base + a.base/shift[diff];
	normalize(result);
      }
      else
	result = *this;
    //
    // Case 2: a is larger => base needs to be shifted
    //
    else 
      if(diff > -max_shift) { 
	result.exponent = a.exponent;
	result.base = base/shift[-diff] + a.base; 
	normalize(result);
      }
      else 
	result = a;

    return result;
  }

  inline ee_t& operator+=(const ee_t& a) { *this = *this+a; return *this; }

};

inline double log(ee_t& a) { return a.log_get(); }

#endif
