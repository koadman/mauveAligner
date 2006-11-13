#ifndef _FFT_H_
#define _FFT_H_

//
// Corresponds to fourl in numerical recipes.
//
void fft(double data[], unsigned long nn, int isign);

//
// Corresponds to realft in numerical recipes.
//
void fft_real(double data[], unsigned long n, int isign);

//
// Corresponds to convlv in numerical recipes.
//
void fft_real_conv(double data[], unsigned long n, double respns[], unsigned long m,
		   int isign, double ans[]);

#endif
