#ifndef FFT_H
#define FFT_H

#define _USE_MATH_DEFINES
#include <math.h>
#include "structure_utils.h"

class cFFT {
private:
	unsigned int N, which;
	unsigned int log_2_N;
	float pi2;
	unsigned int *reversed;
	complex_number **W;
	complex_number *c[2];
protected:
public:
	cFFT(unsigned int N);
	~cFFT();
	unsigned int reverse(unsigned int i);
	complex_number w(unsigned int x, unsigned int N);
	void fft(complex_number* input, complex_number* output, int stride, int offset);
};

#endif