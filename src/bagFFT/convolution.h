#ifndef _CONVOLUTION_H_
#define _CONVOLUTION_H_

#include <math.h>
#include <vector>
#include "complex.h"
#include "fft.h"
#include "boost/array.hpp"
#include "boost/scoped_array.hpp"
#include "boost/shared_ptr.hpp"
//
// Computes the convolution of a[0..(n-1)] and b[0..(to_fill-1)]
// and writes the result in result[0..(to_fill-1)].
// The computation is done "naively" i.e. without using FFT's.
// Note that b and result can be the same vector.
//
template<class T>
void conv(T* a, T* b, T* result, int n, int to_fill) {

  for(int i = to_fill-1; i >= 0; i--) {
    
    result[i] = a[0]*b[i];
    
    for(int j = 1; j <= MIN(i, n-1); j++)              
      result[i] += a[j]*b[i-j];
  }
}

//
// Computes the fft of a[0..(n-1)].
// Based on Numerical recipes code.
//
inline void fftnr(complex* a, int n, int sign) {

  fft(((double*) a)-1, n, sign);
}

//
// Sets result to the convolution of a and b where
// they are all vectors of size n and n is large
// enough to hold the result (also a power of 2).
// Computed using FFTs and based on Numerical recipes code.
//
inline void conv_fftnr(complex* a , complex* b, complex* result, int n) {

  fft(((double*) a)-1, n, -1); fft(((double*) b)-1, n, -1);
      
  for(int x = 0; x < n; x++)
    result[x] = a[x]*b[x]/((double) n);	

  fft(((double*) result)-1, n, 1);
}

//
// Sets result to the convolution of a and b where
// they are all real vectors of size n and n is large
// enough to hold the result.
// Computed using FFTs and based on Numerical recipes code.
//
inline void real_conv_fftnr(double* a , double* b, double* result, int n) {

  fft_real_conv(a-1, n/2, b-1, n/2, 1, result-1);
}

//
// FFTW code. To be compiled if NR = 0.
//

#if NR == 0

#include <fcntl.h>

//#include <unistd.h>
#include "fftw3.h"

static void my_fftw_write_char(char c, void *f) { fputc(c, (FILE *) f); }
#define fftw_export_wisdom_to_file(f) fftw_export_wisdom(my_fftw_write_char, (void*) (f))
#define fftwf_export_wisdom_to_file(f) fftwf_export_wisdom(my_fftw_write_char, (void*) (f))
#define fftwl_export_wisdom_to_file(f) fftwl_export_wisdom(my_fftw_write_char, (void*) (f))

static int my_fftw_read_char(void *f) { return fgetc((FILE *) f); }
#define fftw_import_wisdom_from_file(f) fftw_import_wisdom(my_fftw_read_char, (void*) (f))
#define fftwf_import_wisdom_from_file(f) fftwf_import_wisdom(my_fftw_read_char, (void*) (f))
#define fftwl_import_wisdom_from_file(f) fftwl_import_wisdom(my_fftw_read_char, (void*) (f))

//
// Returns an appropriate structure for
// read and write locks to pass to the function fcntl.
//
/*
struct flock* file_lock(short type, short whence) {

  static struct flock ret;
  ret.l_type = type; ret.l_start = 0; ret.l_whence = whence;
  ret.l_len = 0; ret.l_pid = getpid();
  return &ret;
}
*/
//
// For importing the wisdom file for FFTW
// before computing the first fft.
//
static int wisdom_init = 0;
inline void init_fftw() {

  if(!wisdom_init) {
    FILE* wisdom = fopen("wisdom", "r");
    //int w = open("wisdom", O_RDONLY);
    //fcntl(w, F_SETLKW, file_lock(F_RDLCK, SEEK_SET));
	
    fftw_import_wisdom_from_file(wisdom);
	
    //fcntl(w, F_SETLKW, file_lock(F_RDLCK, SEEK_SET));
    //close(w);
    fclose(wisdom);
    wisdom_init = 1;
  }
}

//
// To update the wisdom file with newly acquired wisdom.
//
inline void save_fftw() {

    FILE* wisdom = fopen("wisdom", "w");
    //int w = open("wisdom", O_WRONLY);
    //fcntl(w, F_SETLKW, file_lock(F_WRLCK, SEEK_SET));

	
    fftw_export_wisdom_to_file(wisdom);
	
    //fcntl(w, F_SETLKW, file_lock(F_WRLCK, SEEK_SET));
    //close(w);
    fclose(wisdom);
}

//
// Computes the fft of a[0..(n-1)].
// Based on FFTW code.
//
inline void fftw(complex* a, int n, int sign) {

  //
  // Getting a plan.
  //
  init_fftw();
  //fftw_complex out[n];
  //boost::array<fftw_complex,n> out;
  //vector<fftw_complex> out;
  //fftw_complex* out = new fftw_complex[n];
  boost::scoped_array<fftw_complex> out (new fftw_complex[n]);

  //cout << "allocating  fftw_complex[1000]" << endl;
  //boost::shared_array<fftw_complex> test (new fftw_complex[1000]);

  fftw_plan fft = fftw_plan_dft_1d(n, (fftw_complex*) a, out.get(), sign, FFTW_MEASURE);
  save_fftw();

  //
  // Executing it.
  //
  fftw_execute(fft);

  //
  // Copying the result back into a.
  //
  for(int i = 0; i < n; i++)
    a[i].set(out[i][0], out[i][1]);

  //delete[](out);
  out.reset();
}

//
// Sets result to the convolution of a and b where
// they are all vectors of size n and n is large
// enough to hold the result.
// Computed using FFTs and based on FFTW code.
//
inline void conv_fftw(complex* a , complex* b, complex* result, int n) {

  static int cur_n = 0;
  static fftw_complex* in = 0;
  static fftw_complex* out1 = 0; static fftw_complex* out2 = 0;
  static fftw_plan fft1 = 0; static fftw_plan fft2 = 0; static fftw_plan ifft = 0;

  vector<complex> a2(n);
  for(int i = 0; i < n; i++)
  {
	a2.push_back(a[i]);
  }

  vector<complex> b2(n);
  for(int i = 0; i < n; i++)
  {
	b2.push_back(b[i]);
  }

  vector<complex> result2(n);
  for(int i = 0; i < n; i++)
  {
	result2.push_back(result[i]);
  }
 
  //
  // Avoids creating new arrays and plans everytime
  // 
  if(n != cur_n) {

    init_fftw();
    if(cur_n != 0) {
      delete[](in); delete[](out1); delete[](out2);
      fftw_destroy_plan(fft1); fftw_destroy_plan(fft2); fftw_destroy_plan(ifft);
    }

    in = new fftw_complex[n]; out1 = new fftw_complex[n]; out2 = new fftw_complex[n];
    //fft1 = fftw_plan_dft_1d(n, (fftw_complex*) a, out1, FFTW_FORWARD, FFTW_MEASURE);
	fft1 = fftw_plan_dft_1d(n, (fftw_complex*) &a2[0], out1, FFTW_FORWARD, FFTW_MEASURE);
    fft2 = fftw_plan_dft_1d(n, (fftw_complex*) b, out2, FFTW_FORWARD, FFTW_MEASURE);
    ifft = fftw_plan_dft_1d(n, in, (fftw_complex*) result, FFTW_BACKWARD, FFTW_MEASURE);
    cur_n = n; save_fftw();
  }

  fftw_execute(fft1);
  fftw_execute(fft2);
  
  for(int i = 0; i < n; i++) {
    in[i][0] = (out1[i][0]*out2[i][0]-out1[i][1]*out2[i][1])/n;
    in[i][1] = (out1[i][1]*out2[i][0]+out1[i][0]*out2[i][1])/n;
  }  

  fftw_execute(ifft);
}

//
// Sets result to the convolution of a and b where
// they are all real vectors of size n and n is large
// enough to hold the result.
// Computed using FFTs and based on FFTW code.
//
inline void real_conv_fftw(double* a , double* b, double* result, int n) {

  static int cur_n = 0;
  static fftw_complex* in = 0; static fftw_complex* out1 = 0; static fftw_complex* out2 = 0; static double* out3 = 0;
  static fftw_plan fft1 = 0; static fftw_plan fft2 = 0; static fftw_plan ifft = 0;

  //
  // Avoids creating new arrays and plans everytime
  // 
  if(n != cur_n) {

    init_fftw();

    if(cur_n != 0) {
      delete[](in); delete[](out1); delete[](out2); delete[](out3);
      fftw_destroy_plan(fft1); fftw_destroy_plan(fft2); fftw_destroy_plan(ifft);
    }

    in = new fftw_complex[n]; out1 = new fftw_complex[n]; out2 = new fftw_complex[n]; out3 = new double[n];

    fft1 = fftw_plan_dft_r2c_1d(n, a, out1, FFTW_MEASURE);
    fft2 = fftw_plan_dft_r2c_1d(n, b, out2, FFTW_MEASURE);
    ifft = fftw_plan_dft_c2r_1d(n, in, out3, FFTW_MEASURE);
    cur_n = n; save_fftw();
  }

  fftw_execute(fft1);
  fftw_execute(fft2);
  
  for(int i = 0; i < n; i++) {
    in[i][0] = (out1[i][0]*out2[i][0]-out1[i][1]*out2[i][1])/n;
    in[i][1] = (out1[i][1]*out2[i][0]+out1[i][0]*out2[i][1])/n;
  }  

  fftw_execute(ifft);
  
  for(int i = 0; i < n; i++) 
    result[i] = out3[i];
}

#endif

#endif
