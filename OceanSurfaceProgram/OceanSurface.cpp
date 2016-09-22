#include "OceanSurface.h"
#define _USE_MATH_DEFINES // for M_PI constants
#include <math.h>
#include <random>
#include <cassert>

// realization of OceanSurface class

OceanSurface::OceanSurface(int N, int M, float L_x, float L_z, float A, glm::vec2 _wind, float g, bool _lines, bool _gpu) : N(N), M(M), L_x(L_x),
L_z(L_z), A(A), wind(_wind), g(g), lines(_lines), gpu(_gpu) {
	grid = new surface_vertex[N * M]; // construct wave height field grid
	init_positions = new vertex_2d[N * M];
	grid_tex_cord = new tex_img_coord[N * N];
	h_t0 = new complex_number[N * M];
	h_t0_cc = new complex_number[N * M];
	h_fft = new complex_number[N * M];
	dx_fft = new complex_number[N * M];
	dz_fft = new complex_number[N * M];
	gradient_x = new complex_number[N * M];
	gradient_z = new complex_number[N * M];
	lambda = -1.0f; // numeric constant for calculating displacement
	num_indices = 0;
	num_triangle_indices = 0;
	indices = new GLuint[(N - 1) * (M - 1) * 6 + 2 * (N + M - 2)];
	triangle_indices = new GLuint[(N - 1) * (M - 1) * 6 * 10];
	if (gpu) {
		wind = glm::vec2(wind.y, wind.x);
	}
	//fft
	myFFT = new MyFFT(N, M);

	// compute shader FFT stuff
	char *s_name = "compute_shaders/fft_compute.glsl";
	if (N == 32) {
		s_name = "compute_shaders/fft_compute_32.glsl";
	}
	if (N == 64) {
		s_name = "compute_shaders/fft_compute_64.glsl";
	}
	if (N == 128) {
		s_name = "compute_shaders/fft_compute_128.glsl";
	}
	if (N == 256) {
		s_name = "compute_shaders/fft_compute_256.glsl";
	}
	if (N == 512) {
		s_name = "compute_shaders/fft_compute_512.glsl";
	}
	fft_comp_program = create_comp_program_from_file(s_name);

	int bits = log2(N * 1.0f);
	calc_binary_reverse(bits);
	glGenBuffers(1, &rev_ind_buffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, rev_ind_buffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, N * sizeof(int), binary_reverse, GL_STATIC_DRAW);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	GLint steps_location = glGetUniformLocation(fft_comp_program, "steps");
	fft_column_location = glGetUniformLocation(fft_comp_program, "fft_column");

	glUseProgram(fft_comp_program);
	glUniform1i(steps_location, bits);

	// texture for inpute height map
	glGenTextures(1, &texture_H_t);
	glBindTexture(GL_TEXTURE_2D, texture_H_t);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, N, N, 0, GL_RG, GL_FLOAT, 0);
	/*glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);*/

	// texture for output height map per rows
	glGenTextures(1, &texture_H_fft_t_row);
	glBindTexture(GL_TEXTURE_2D, texture_H_fft_t_row);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);
	/*glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);*/

	// texture for output height map per columns
	glGenTextures(1, &texture_H_fft_t_col);
	glBindTexture(GL_TEXTURE_2D, texture_H_fft_t_col);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);
	/*glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);*/

	// texture for fft dx
	glGenTextures(1, &tex_dx_fft_row);
	glBindTexture(GL_TEXTURE_2D, tex_dx_fft_row);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);
	glGenTextures(1, &tex_dx_fft);
	glBindTexture(GL_TEXTURE_2D, tex_dx_fft);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);

	// texture for fft dz
	glGenTextures(1, &tex_dz_fft_row);
	glBindTexture(GL_TEXTURE_2D, tex_dz_fft_row);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);
	glGenTextures(1, &tex_dz_fft);
	glBindTexture(GL_TEXTURE_2D, tex_dz_fft);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);

	// texture for fft gradx
	glGenTextures(1, &tex_gradx_fft_row);
	glBindTexture(GL_TEXTURE_2D, tex_gradx_fft_row);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);
	glGenTextures(1, &tex_gradx_fft);
	glBindTexture(GL_TEXTURE_2D, tex_gradx_fft);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);

	// texture for fft gradz
	glGenTextures(1, &tex_gradz_fft_row);
	glBindTexture(GL_TEXTURE_2D, tex_gradz_fft_row);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);
	glGenTextures(1, &tex_gradz_fft);
	glBindTexture(GL_TEXTURE_2D, tex_gradz_fft);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);

	// unbind all textures
	glBindTexture(GL_TEXTURE_2D, 0);

	// assign each point its original 3d coordinates
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++) {
			int pos = i * M + j; // one-dimensional pos from linear two-dimensional grid

			grid_tex_cord[pos] = tex_img_coord(i, j);

			// precomputer amplitudes at each point
			h_t0[pos] = h_t_0(i, j);
			h_t0_cc[pos] = h_t_0(-i, -j).cc();

			// we project our N*M grid coordinates to real world L_x*L_z coordinates
			grid[pos].x = (i - (N >> 1)) * L_x / N;
			grid[pos].y = 0.0f; // height is 0
			grid[pos].z = (j - (M >> 1)) * L_z / M;

			// save initial positions
			init_positions[pos].x = grid[pos].x;
			init_positions[pos].z = grid[pos].z;

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

			// construct indices for triangle index drawing
			if (i != N - 1 && j != M - 1) {
				triangle_indices[num_triangle_indices++] = pos;
				triangle_indices[num_triangle_indices++] = pos + N + 1;
				triangle_indices[num_triangle_indices++] = pos + 1;

				triangle_indices[num_triangle_indices++] = pos;
				triangle_indices[num_triangle_indices++] = pos + N;
				triangle_indices[num_triangle_indices++] = pos + N + 1;
			}
		}

	// check that its okay
	assert(num_indices == (N - 1) * (M - 1) * 6 + 2 * (N + M - 2));
	//assert(num_triangle_indices, (N - 1) * (M - 1) * 6);

	prepare_for_pipeline();

	// update height map compute shader program
	upd_height_program = create_comp_program_from_file("compute_shaders/update_height_map.glsl");
	GLint len_location = glGetUniformLocation(upd_height_program, "L");;
	GLint n_location = glGetUniformLocation(upd_height_program, "N");
	total_time_location = glGetUniformLocation(upd_height_program, "totalTime");
	glUseProgram(upd_height_program);
	glUniform1f(len_location, L_x);
	glUniform1i(n_location, N);
	glUseProgram(0);

	// height map textures
	glGenTextures(1, &texture_H_t0);
	glGenTextures(1, &texture_H_t0_cc);
	glBindTexture(GL_TEXTURE_2D, texture_H_t0);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, N, N, GL_RG, GL_FLOAT, h_t0);
	glBindTexture(GL_TEXTURE_2D, texture_H_t0_cc);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, N, N, GL_RG, GL_FLOAT, h_t0_cc);
	glBindTexture(GL_TEXTURE_2D, 0);

	// displacement map textures
	glGenTextures(1, &texture_Dx);
	glGenTextures(1, &texture_Dz);
	glBindTexture(GL_TEXTURE_2D, texture_Dx);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);
	glBindTexture(GL_TEXTURE_2D, texture_Dz);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);
	glBindTexture(GL_TEXTURE_2D, 0);

	// gradient map textures
	glGenTextures(1, &texture_Gradx);
	glGenTextures(1, &texture_Gradz);
	glBindTexture(GL_TEXTURE_2D, texture_Gradx);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);
	glBindTexture(GL_TEXTURE_2D, texture_Gradz);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, N, N);
	glBindTexture(GL_TEXTURE_2D, 0);

	// Fresnel textures
	int fres_size = 128;
	float g_SkyBlending = 16.0f;
	unsigned int *buffer = new unsigned int[fres_size];
	for (int i = 0; i < fres_size; i++)
	{
		float cos_a = i / (float)fres_size;
		// Using water's refraction index 1.33
		unsigned int fresnel = (unsigned int)(D3DXFresnelTerm(cos_a, 1.33f) * 255);

		unsigned int sky_blend = (unsigned int)(powf(1 / (1 + cos_a), g_SkyBlending) * 255);

		buffer[i] = (sky_blend << 8) | fresnel;
	}
	glGenTextures(1, &tex_Fresnel);
	glBindTexture(GL_TEXTURE_1D, tex_Fresnel);
	glTexStorage1D(GL_TEXTURE_1D, 1, GL_RGBA8, fres_size);
	glTexSubImage1D(GL_TEXTURE_1D, 0, 0, fres_size, GL_RGBA, GL_UNSIGNED_INT, buffer);

	const GLfloat border[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glTexParameterfv(GL_TEXTURE_1D, GL_TEXTURE_BORDER_COLOR, border);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_COMPARE_FUNC, GL_NEVER);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

	glBindTexture(GL_TEXTURE_1D, 0);
	delete[] buffer;

	// Reflect cube textures
	tex_ReflectCube = SOIL_load_OGL_single_cubemap
		(
			"textures/sky_cube.dds",
			SOIL_DDS_CUBEMAP_FACE_ORDER,
			SOIL_LOAD_AUTO,
			SOIL_CREATE_NEW_ID,
			SOIL_FLAG_MIPMAPS | SOIL_FLAG_DDS_LOAD_DIRECT
		);
}

float OceanSurface::D3DXFresnelTerm(float costheta, float refractionindex) {
	float a, d, g, result;

	g = sqrtf(refractionindex * refractionindex + costheta * costheta - 1.0f);
	a = g + costheta;
	d = g - costheta;
	result = (costheta * a - 1.0f) * (costheta * a - 1.0f) / ((costheta * d + 1.0f) * (costheta * d + 1.0f)) + 1.0f;
	result *= 0.5f * d * d / (a * a);
	
	return result;
}

OceanSurface::~OceanSurface() {
	if (grid)
		delete[] grid;
	if (indices)
		delete[] indices;
}

void OceanSurface::calc_binary_reverse(int bits) {
	binary_reverse = new int[N];

	for (int i = 0; i < N; i++) {
		int nrev = i, n = i;
		for (int j = 1; j < bits; j++) {
			n >>= 1;
			nrev <<= 1;
			nrev |= n & 1;   // give LSB of n to nrev
		}
		nrev &= N - 1;

		binary_reverse[i] = nrev;
	}
}

void OceanSurface::prepare_for_pipeline() {
	vao = indices_vbo = points_vbo = triangle_indices_vbo = img_coord_vbo = 0;

	// img coords
	glGenBuffers(1, &img_coord_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, img_coord_vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(tex_img_coord) * N * N, grid_tex_cord, GL_STATIC_DRAW);
	
	// points vbo
	glGenBuffers(1, &points_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(surface_vertex) * N * M, grid, gpu ? GL_STATIC_DRAW : GL_DYNAMIC_DRAW);

	// indices vbo
	glGenBuffers(1, &indices_vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indices_vbo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * num_indices, indices, GL_STATIC_DRAW);

	// triangle indices vbo
	glGenBuffers(1, &triangle_indices_vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triangle_indices_vbo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)* num_triangle_indices, triangle_indices, GL_STATIC_DRAW);

	// vao
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
	// assign position to location = 0
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(surface_vertex), 0);
	// assign normal to location = 1
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(surface_vertex), (GLvoid *)12);
	// assign img coord to location = 2
	glBindBuffer(GL_ARRAY_BUFFER, img_coord_vbo);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(tex_img_coord), 0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

float OceanSurface::phillips_spectrum(int n, int m) {
	/*
	glm::vec2 k = glm::vec2(M_PI * (2.0f * n - N) / L_x, M_PI * (2.0f * m - M) / L_z); // k is wavevector
	if (n < 0 && m < 0) {
		k = glm::vec2(M_PI * (2.0f * n + N) / L_x, M_PI * (2.0f * m + M) / L_z); // k is wavevector
	}
	float len_k = glm::length(k);
	// cut off very long waves
	if (len_k < 0.000001) {
		return 0.0;
	}

	glm::vec2 k_norm = glm::normalize(k);
	glm::vec2 wind_norm = glm::normalize(wind);

	float L = glm::length(wind) * glm::length(wind) / g;
	float l = L / 1000.0f; // coefficent to suppress small waves
	return A * exp(-1.0f / pow(len_k * L, 2.0f)) / pow(len_k, 4.0f) * pow(glm::dot(k_norm, wind_norm), 2.0f) * exp(-pow(len_k * l, 2.0f));
	*/
	glm::vec2 k = glm::vec2(M_PI * (2.0f * n - N) / L_x, M_PI * (2.0f * m - M) / L_z); // k is wavevector
	if (n < 0 && m < 0) {
		k = glm::vec2(M_PI * (2.0f * n + N) / L_x, M_PI * (2.0f * m + M) / L_z); // k is wavevector
	}
	float k_length = glm::length(k);
	if (k_length < 0.000001) return 0.0;

	float k_length2 = k_length  * k_length;
	float k_length4 = k_length2 * k_length2;

	//float k_dot_w = glm::dot(k, wind);
	float k_dot_w = glm::dot(glm::normalize(k), glm::normalize(wind));
	float k_dot_w2 = k_dot_w * k_dot_w;

	float w_length = glm::length(wind);
	float L = w_length * w_length / g;
	float L2 = L * L;

	float damping = 0.001;
	float l2 = L2 * damping * damping;

	float res =  A * exp(-1.0f / (k_length2 * L2)) / k_length4 * k_dot_w2 * exp(-k_length2 * l2);
	return res;
}

complex_number OceanSurface::h_t_0(int n, int m) {
	static std::random_device rd;
	static std::mt19937 generator(rd());
	static std::normal_distribution<float> distribution(0.0f, 1.0f);

	float e_r = distribution(generator);
	float e_i = distribution(generator);
	//e_r = 0.02f; e_i = 0.02f;

	return complex_number(e_r, e_i) * sqrt(phillips_spectrum(n, m) / 2.0f);
}

float OceanSurface::dispersion_relation(int n, int m) {
	float w_0 = 2.0f * M_PI / 200.0f;
	glm::vec2 k = glm::vec2(M_PI * (2.0f * n - N) / L_x, M_PI * (2.0f * m - M) / L_z); // k is wavevector
	//return sqrt(g * glm::length(k));
	return floor(sqrt(g * glm::length(k)) / w_0) * w_0;
}

complex_number OceanSurface::h_t(int n, int m, float t) {
	float w_freq = dispersion_relation(n, m);
	float wt = w_freq * t;

	complex_number c1(cos(wt), sin(wt));
	complex_number c2(cos(wt), -sin(wt));

	int pos = n * M + m;

	return h_t0[pos] * c1 + h_t0_cc[pos] * c2;
	//return h_t0[pos] * c1 + h_t0[pos2].cc() * c2;
	//return h_t0[pos] * c1;
}

float OceanSurface::h(float x, float z, float t) {
	complex_number height;
	glm::vec2 x_v = glm::vec2(x, z);

	for (int i = 0; i < N; i++) {
		float k_x = M_PI * (2.0f * i - N) / L_x;
		for (int j = 0; j < M; j++) {
			int pos = i * M + j;
			float k_z = M_PI * (2.0f * j - M) / L_z;

			glm::vec2 k = glm::vec2(k_x, k_z); // k is wavevector
			float d = glm::dot(k, x_v);
			complex_number h_0_t = h_t(i, j, t);
			complex_number c1(cos(d), sin(d));

			height = height + c1 * h_0_t;
		}
	}

	return height.re;
}

void OceanSurface::updateOceanSlow(float t) {
	// update wave height
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			int pos = i * M + j;

			grid[pos].y = h(grid[pos].x, grid[pos].z, t);
			//grid[pos].y = t;
		}
	}
}

void OceanSurface::updateOcean(float t) {
	if (!gpu) {

		for (int i = 0; i < N; i++) {
			float k_x = float(M_PI) * (2.0f * i - N) / L_x; // x-coordinate of wavevector at 'pos'

			for (int j = 0; j < M; j++) {
				int pos = i * M + j;

				float k_z = M_PI * (2.0f * j - M) / L_z; // z-coordinate of wavevector at 'pos'
				float k_length = sqrtf(k_x * k_x + k_z * k_z); // length of wavevector 'k'

				// calc all complex amplitudes in time t
				h_fft[pos] = h_t(i, j, t);

				// calc all displacements in time t
				if (k_length < 0.000001) { // ignore small wavevectors, put (0, 0) displacement
					dx_fft[pos] = complex_number();
					dz_fft[pos] = complex_number();
				}
				else { // calculate displacement via formula
					dx_fft[pos] = complex_number(0.0f, -k_x / k_length) * h_fft[pos];
					dz_fft[pos] = complex_number(0.0f, -k_z / k_length) * h_fft[pos];
				}

				// calc x-, and z- component of gradient of h_fft
				gradient_x[pos] = complex_number(0.0f, k_x) * h_fft[pos];
				gradient_z[pos] = complex_number(0.0f, k_z) * h_fft[pos];
			}
		}

		// row flipping
		/*/complex_number *h_fft_2 = new complex_number[N * N];
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < M; j++) {
				int pos_1 = i * M + j;
				int pos_2 = (N - i - 1) * M + j;
				h_fft_2[pos_2] = h_fft[pos_1];
			}
		}*/
	} 
	else {

		// first update height map with time amplitudes (and displacement map with gradient map);
		glUseProgram(upd_height_program);

		glUniform1f(total_time_location, t);
		// height map
		glBindImageTexture(0, texture_H_t0, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(1, texture_H_t0_cc, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(2, texture_H_t, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);

		// displacement map
		glBindImageTexture(3, texture_Dx, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
		glBindImageTexture(4, texture_Dz, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);

		// gradient map
		glBindImageTexture(5, texture_Gradx, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
		glBindImageTexture(6, texture_Gradz, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);

		glDispatchCompute(N, N, 1);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		glUseProgram(0);

		// fft in compute shader
		glUseProgram(fft_comp_program);

		// update texture
		/*glBindTexture(GL_TEXTURE_2D, texture_H_t);
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, N, N, GL_RG, GL_FLOAT, h_fft);
		//glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, N, N, 0, GL_RG, GL_FLOAT, h_fft);
		glBindTexture(GL_TEXTURE_2D, 0);*/


		// bind reverse indices
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, rev_ind_buffer);

		/*______________HEIGHT MAP FFT________________*/

		// fft per rows
		glUniform1i(fft_column_location, 0);
		glBindImageTexture(0, texture_H_t, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(1, texture_H_fft_t_row, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
		glDispatchCompute(1, N, 1);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		// fft per columns
		glUniform1i(fft_column_location, 1);
		glBindImageTexture(0, texture_H_fft_t_row, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(1, texture_H_fft_t_col, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
		glDispatchCompute(1, N, 1);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		/*______________DISPLACEMENT-X MAP FFT________________*/

		// fft per rows
		glUniform1i(fft_column_location, 0);
		glBindImageTexture(0, texture_Dx, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(1, tex_dx_fft_row, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
		glDispatchCompute(1, N, 1);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		// fft per columns
		glUniform1i(fft_column_location, 1);
		glBindImageTexture(0, tex_dx_fft_row, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(1, tex_dx_fft, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
		glDispatchCompute(1, N, 1);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		/*______________DISPLACEMENT-Z MAP FFT________________*/

		// fft per rows
		glUniform1i(fft_column_location, 0);
		glBindImageTexture(0, texture_Dz, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(1, tex_dz_fft_row, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
		glDispatchCompute(1, N, 1);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		// fft per columns
		glUniform1i(fft_column_location, 1);
		glBindImageTexture(0, tex_dz_fft_row, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(1, tex_dz_fft, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
		glDispatchCompute(1, N, 1);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		/*______________GRADIENT-X MAP FFT________________*/

		// fft per rows
		glUniform1i(fft_column_location, 0);
		glBindImageTexture(0, texture_Gradx, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(1, tex_gradx_fft_row, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
		glDispatchCompute(1, N, 1);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		// fft per columns
		glUniform1i(fft_column_location, 1);
		glBindImageTexture(0, tex_gradx_fft_row, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(1, tex_gradx_fft, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
		glDispatchCompute(1, N, 1);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		/*______________GRADIENT-Z MAP FFT________________*/

		// fft per rows
		glUniform1i(fft_column_location, 0);
		glBindImageTexture(0, texture_Gradz, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(1, tex_gradz_fft_row, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
		glDispatchCompute(1, N, 1);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		// fft per columns
		glUniform1i(fft_column_location, 1);
		glBindImageTexture(0, tex_gradz_fft_row, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(1, tex_gradz_fft, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RG32F);
		glDispatchCompute(1, N, 1);
		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		glUseProgram(0);
	}

	int error = glGetError();
	if (error != GL_NO_ERROR) {
		int b = 4;
	}

	if (gpu) {
		return;
	}
	for (unsigned int i = 0; i < N; i++) { // horizontal fft for rows
		/*fft->fft(h_fft, h_fft, 1, i * M);
		fft->fft(dx_fft, dx_fft, 1, i * M);
		fft->fft(dz_fft, dz_fft, 1, i * M);
		fft->fft(gradient_x, gradient_x, 1, i * M);
		fft->fft(gradient_z, gradient_z, 1, i * M);*/
		myFFT->processVertical(h_fft, 1, i * M);
		myFFT->processVertical(dx_fft, 1, i * M);
		myFFT->processVertical(dz_fft, 1, i * M);
		myFFT->processVertical(gradient_x, 1, i * M);
		myFFT->processVertical(gradient_z, 1, i * M);
	}

	for (unsigned int j = 0; j < M; j++) { // vertical fft for columns
		/*fft->fft(h_fft, h_fft, M, j);
		fft->fft(dx_fft, dx_fft, M, j);
		fft->fft(dz_fft, dz_fft, M, j);
		fft->fft(gradient_x, gradient_x, M, j);
		fft->fft(gradient_z, gradient_z, M, j);*/
		myFFT->processHorizontal(h_fft, M, j);
		myFFT->processHorizontal(dx_fft, M, j);
		myFFT->processHorizontal(dz_fft, M, j);
		myFFT->processHorizontal(gradient_x, M, j);
		myFFT->processHorizontal(gradient_z, M, j);
	}

	
	glm::vec3 normal; // temp variable for normal calculation
	int sign;
	float signs[] = { 1.0f, -1.0f };
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			int pos = i * M + j;

			// flip sign for height
			sign = signs[(i + j) & 1];
			grid[pos].y = h_fft[pos].re * sign;

			// flip sign for displacement
			dx_fft[pos] = dx_fft[pos] * sign;
			dz_fft[pos] = dz_fft[pos] * sign;
			// then calc displacement for original position in grid
			//grid[pos].x = init_positions[pos].x + lambda * dx_fft[pos].re;
			//grid[pos].z = init_positions[pos].z + lambda * dz_fft[pos].re;

			// flip sign for gradient component;
			gradient_x[pos] = gradient_x[pos] * sign;
			gradient_z[pos] = gradient_z[pos] * sign;

			// then calculate normal of vertex and normalize it
			normal = glm::vec3(-gradient_x[pos].re, 1.0f, -gradient_z[pos].re);
			normal = glm::normalize(normal);
			// and save it to vertex structure
			grid[pos].normal_x = normal.x;
			grid[pos].normal_y = normal.y;
			grid[pos].normal_z = normal.z;
		}
	}
}

void OceanSurface::render() {
	// fresnel texture
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_1D, tex_Fresnel);
	
	// reflect texture
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_CUBE_MAP, tex_ReflectCube);

	// bind fft height map
	if (gpu) {
		glBindImageTexture(0, texture_H_fft_t_col, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(1, tex_dx_fft, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(2, tex_dz_fft, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(3, tex_gradx_fft, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
		glBindImageTexture(4, tex_gradz_fft, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RG32F);
	}
	//glActiveTexture(GL_TEXTURE0);
	//glBindTexture(GL_TEXTURE_2D, texture_H_t);
	//glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, N, N, GL_RG, GL_FLOAT, h_fft);
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, N, N, 0, GL_RG, GL_FLOAT, h_fft);
	//glBindTexture(GL_TEXTURE_2D, 0);

	if (!gpu) {
		// update vertices in gpu array buffer
		glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(surface_vertex)* N * M, grid);
	}

	// draw on screen
	glBindVertexArray(vao);

	if (lines) {
		// draw lines
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indices_vbo);
		glDrawElements(GL_LINES, num_indices, GL_UNSIGNED_INT, 0);
	}
	else {
		// draw triangles
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triangle_indices_vbo);
		glDrawElements(GL_TRIANGLES, num_triangle_indices, GL_UNSIGNED_INT, 0);
	}
}