#pragma once

#include <GL\glew.h>

struct surface_vertex {
	GLfloat x, y, z; // x, z - horizontal position, y - height. We pass them to OpenGL pipeline, 
					//thats why GLfloat instead of float
	GLfloat normal_x, normal_y, normal_z; // (x, y, z) components of normal vector
};

struct vertex_2d {
	GLfloat x, z;
};

struct tex_img_coord {
	GLfloat x, z;

	tex_img_coord() { x = 0; z = 0; }
	tex_img_coord(GLfloat _x, GLfloat _z) { x = _x; z = _z; }
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

	complex_number operator -(const complex_number& b) const {
		return complex_number(re - b.re, im - b.im);
	}

	complex_number operator *(const complex_number& b) const {
		return complex_number(re * b.re - im * b.im, re * b.im + im * b.re);
	}

	// return complex conjugate
	complex_number cc() {
		return complex_number(re, -im);
	}
};