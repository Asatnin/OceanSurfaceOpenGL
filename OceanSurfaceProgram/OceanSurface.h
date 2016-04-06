#pragma once

#include "structure_utils.h"
#include <glm\vec2.hpp>
#include <glm\geometric.hpp>

// class for ocean surface definition

class OceanSurface {
	private:
		int N, M; // grid size for wave height field (should be power of two both)
		float L_x, L_z; // real world dimensions for wave height field (in meters)
		surface_vertex *grid; // each point in grid consists of 3d-world coordinates (x, y, z)
		complex_number *h_t0, *h_t0_cc; // array of precomputed random complex amplitudes
		float A; // numeric constant for wave height in phillips spectrum
		glm::vec2 wind; // vector of wind direction
		float g; // gravitationl constant
		GLuint *indices; // indices for index drawing of grid
		int num_indices; // number of indices
		GLuint indices_vbo, points_vbo, vao; // vbo for indices and vertices, vao for vertices
		void prepare_for_pipeline();

	public:
		OceanSurface(int N, int M, float L_x, float L_z, float A, glm::vec2 wind, float g); // constructor
		~OceanSurface(); // destructor
		// here and further n = n' and m = m' (equally to formulas for n'=0..N-1, m'=0..M-1)
		float phillips_spectrum(int n, int m);
		complex_number h_t_0(int n, int m);
		float dispersion_relation(int n, int m);
		complex_number h_t(int n, int m, float t);
		float h(float x, float z, float t);
		void render(float t);
};