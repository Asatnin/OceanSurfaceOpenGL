#pragma once

#include <GL\glew.h>

struct surface_vertex {
	GLfloat x, y, z; // x, z - horizontal position, y - height. We pass them to OpenGL pipeline, 
					//thats why GLfloat instead of float
};

struct complex_number {
	float re, im;

	complex_number() { re = 0.0f; im = 0.0f; }
	complex_number(float r, float i) { re = r; im = i; }

	complex_number operator *(const float b) const {
		return complex_number(b * re, b * im);
	}

	complex_number operator +(const complex_number& b) const {
		return complex_number(b.re + re, b.im + im);
	}

	complex_number operator *(const complex_number& b) const {
		return complex_number(re * b.re - im * b.im, re * b.im + im * b.re);
	}

	// return complex conjugate
	complex_number cc() {
		return complex_number(re, -im);
	}
};