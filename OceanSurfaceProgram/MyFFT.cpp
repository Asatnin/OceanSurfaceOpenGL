#include "MyFFT.h"
#define _USE_MATH_DEFINES // for M_PI constants
#include <math.h>

MyFFT::MyFFT(unsigned int N, unsigned M) : size_N(N), size_M(M) {
	bits_N = log2(N * 1.0f);
	bits_M = log2(M * 1.0f);
	binary_reverseN = new unsigned int[N];
	binary_reverseM = new unsigned int[M];
	calc_binary_reverse(binary_reverseN, N, bits_N);
	calc_binary_reverse(binary_reverseM, M, bits_M);
}

void MyFFT::calc_binary_reverse(unsigned int *binary_reverse, unsigned int N, unsigned int bits) {
	for (unsigned int i = 0; i < N; i++) {
		unsigned int nrev = i, n = i;
		for (unsigned int j = 1; j < bits; j++) {
			n >>= 1;
			nrev <<= 1;
			nrev |= n & 1;   // give LSB of n to nrev
		}
		nrev &= N - 1;

		binary_reverse[i] = nrev;
	}
}

// vertical for N-size arrays
void MyFFT::processVertical(complex_number *in, unsigned int stride, unsigned int offset) {
	// rearrange data in binary reverse order
	complex_number *data = new complex_number[size_N];
	for (unsigned int i = 0; i < size_N; i++) {
		data[i] = in[binary_reverseN[i] * stride + offset];
	}

	process(data, size_N);

	// copy back
	for (unsigned int i = 0; i < size_N; i++) {
		in[i  * stride + offset] = data[i];
	}
}

// horizontal for M-size arrays
void MyFFT::processHorizontal(complex_number *in, unsigned int stride, unsigned int offset) {
	// rearrange data in binary reverse order
	complex_number *data = new complex_number[size_M];
	for (unsigned int i = 0; i < size_M; i++) {
		data[i] = in[binary_reverseM[i] * stride + offset];
	}

	process(data, size_M);

	// copy back
	for (unsigned int i = 0; i < size_M; i++) {
		in[i  * stride + offset] = data[i];
	}
}

// data is already rearranged in binary reverse order
void MyFFT::process(complex_number *data, unsigned int size) {
	// steps has values: 1, 2, 4, 8...
	for (unsigned int step = 1; step < size; step <<= 1) { // that loop run log2(size) times. step is half period of exp
		// jump has values: 2, 4, 8, 16...
		unsigned int jump = step << 1;
		// delta has values: PI, PI / 2, PI / 4
		float delta = float(M_PI) / step; // exp has form: exp(i * PI / step * x) == exp(i * delta * x),
		// therefore angle of exp has stride 'delta' with increasing 'x' by one (x = 0..N - 1 as DFT)

		// mult is for following: exp(i * (ANGLE + delta)) = exp(i * ANGLE) + exp(i * ANGLE) * mult
		const complex_number mult(-2.0f * sinf(delta / 2) * sinf(delta / 2), sinf(delta));
		// therefore mult is used for moving angle from ANGLE by delta counterclockwise

		complex_number factor(1.0f, 0.0f);//twiddle factor exponent is initialized with exp(i * 0), i.e ANGLE = 0

		// we loop through half period of exponent
		// each group has exponent with the same angle but opposite signs
		for (unsigned int group = 0; group < step; group++) {
			// we iterate over pairs with identical exponent factor (twiddle factors)
			for (unsigned int pair = group; pair < size; pair += jump) {
				// pair and match is two buttefly sum terms
				unsigned int match = pair + step;
				// product is multiplication second term of butterfly on twiddle factor
				complex_number product = data[match] * factor;

				// process butterfly operation
				data[match] = data[pair] - product;
				data[pair] = data[pair] + product;
			}

			// increase ANGLE by delta counterclockwise
			factor = factor + factor * mult;
		}
	}
}