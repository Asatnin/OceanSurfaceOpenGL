#include "OceanSurface.h"
#define _USE_MATH_DEFINES // for M_PI constants
#include <math.h>
#include <random>
#include <cassert>

// realization of OceanSurface class

OceanSurface::OceanSurface(int N, int M, float L_x, float L_z, float A, glm::vec2 wind, float g) : N(N), M(M), L_x(L_x),
L_z(L_z), A(A), wind(wind), g(g) {
	grid = new surface_vertex[N * M]; // construct wave height field grid
	num_indices = 0;
	indices = new GLuint[(N - 1) * (M - 1) * 6 + 2 * (N + M - 2)];

	// assign each point its original 3d coordinates
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++) {
			int pos = i * M + j; // one-dimensional pos from linear two-dimensional grid

			// we project our N*M grid coordinates to real world L_x*L_z coordinates
			grid[pos].x = (i - (N >> 1)) * L_x / N;
			grid[pos].y = 0; // height is 0
			grid[pos].z = j - (M >> 1) * L_z / M;

			// construct indices for index drawing
			if (i != N - 1 && j != M - 1) {
				indices[num_indices++] = pos;
				indices[num_indices++] = pos + 1;
				indices[num_indices++] = pos;
				indices[num_indices++] = pos + M;
				indices[num_indices++] = pos;
				indices[num_indices++] = pos + M + 1;
			}
			if (i != N - 1 && j == M - 1) {
				indices[num_indices++] = pos;
				indices[num_indices++] = pos + M;
			}
			if (i == N - 1 && j != M - 1) {
				indices[num_indices++] = pos;
				indices[num_indices++] = pos + 1;
			}
		}

	// check that its okay
	assert(num_indices == (N - 1) * (M - 1) * 6 + 2 * (N + M - 2));

	prepare_for_pipeline();
}

OceanSurface::~OceanSurface() {
	if (grid)
		delete[] grid;
	if (indices)
		delete[] indices;
}

void OceanSurface::prepare_for_pipeline() {
	vao = indices_vbo = points_vbo = 0;
	
	// points vbo
	glGenBuffers(1, &points_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(surface_vertex) * N * M, grid, GL_DYNAMIC_DRAW);

	// indices vbo
	glGenBuffers(1, &indices_vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indices_vbo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * num_indices, indices, GL_STATIC_DRAW);

	// vao
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(0);
}

float OceanSurface::phillips_spectrum(int n, int m) {
	glm::vec2 k = glm::vec2(M_PI * ((n << 1) - N) / L_x, M_PI * ((m << 1) - M) / L_z); // k is wavevector
	float len_k = k.length();
	// cut off very long waves
	if (len_k < 0.000001) {
		return 0.0;
	}

	float L = wind.length() * wind.length() / g;
	float l = L / 1000.0f; // coefficent to suppress small waves
	return A * exp(-1.0f / pow(len_k * L, 2)) / pow(len_k, 4) * pow(glm::dot(k, wind), 2) * exp(-pow(len_k * l, 2));
}

complex_number OceanSurface::h_t_0(int n, int m) {
	static std::default_random_engine generator;
	static std::normal_distribution<float> distribution(0.0f, 1.0f);

	float e_r = distribution(generator);
	float e_i = distribution(generator);

	return complex_number(e_r, e_i) * sqrt(phillips_spectrum(n, m) / 2.0f);
}

float OceanSurface::dispersion_relation(int n, int m) {
	glm::vec2 k = glm::vec2(M_PI * ((n << 1) - N) / L_x, M_PI * ((m << 1) - M) / L_z); // k is wavevector
	return sqrt(g * k.length());
}

complex_number OceanSurface::h_t(int n, int m, float t) {
	float w_freq = dispersion_relation(n, m);
	return h_t_0(n, m) * 2.0f * cos(w_freq * t);
}

float OceanSurface::h(float x, float z, float t) {
	complex_number height;
	glm::vec2 x_v = glm::vec2(x, z);

	for (int i = 0; i < N; i++) {
		float k_x = M_PI * ((i << 1) - N) / L_x;
		for (int j = 0; j < M; j++) {
			int pos = i * M + j;
			float k_z = M_PI * ((j << 1) - M) / L_z;

			glm::vec2 k = glm::vec2(k_x, k_z); // k is wavevector
			float d = glm::dot(k, x_v);
			complex_number h_0_t = h_t(i, j, t);
			complex_number h_0_t_cc = h_t(-i, -j, t).cc();
			complex_number c1(cos(d), sin(d));
			complex_number c2(cos(d), -sin(d));

			height = height + c1 * h_0_t + c2 * h_0_t_cc;
		}
	}

	return height.re;
}

void OceanSurface::render(float t) {
	glBindVertexArray(vao);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indices_vbo);
	glDrawElements(GL_LINES, num_indices, GL_UNSIGNED_INT, 0);
}