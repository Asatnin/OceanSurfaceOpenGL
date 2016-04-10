#pragma once

#include "structure_utils.h"

class MyFFT {
	private:
		unsigned int size_N, size_M; // input and output size, should be power of 2
		// size_N = vertical size, size_M = horizontal size
		unsigned int *binary_reverseN, *binary_reverseM; // binary reversed indices of number 0..N-1
		void calc_binary_reverse(unsigned int *binary_reverse, unsigned int N, unsigned int bits);
		unsigned int bits_N, bits_M;
		void process(complex_number *data, unsigned int size);

	public:
		MyFFT(unsigned int N, unsigned M);
		~MyFFT();
		// *in is 2-dim array, so we need 'stride' and 'offset' to iterate over rows and columns
		void processHorizontal(complex_number *in, unsigned int stride, unsigned int offset);
		void processVertical(complex_number *in, unsigned int stride, unsigned int offset);
};