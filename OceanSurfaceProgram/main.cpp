#include "gl_utils.h"
#include <GL\glew.h>
#define GLFW_DLL
#include <GLFW\glfw3.h>
#include <cassert>
#include <cstdio>
#include <glm\glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm\vec3.hpp>
#include <glm\mat4x4.hpp>
#include <glm\gtc\matrix_transform.hpp>
#include "OceanSurface.h"

int g_width = 800;
int g_height = 600;
GLFWwindow *g_window = NULL;

int main() {
	// prepare logging
	assert(restart_gl_log());
	assert(start_gl());

	// create ocean
	OceanSurface ocean = OceanSurface(64, 64, 16.0f, 16.0f, 0.0005f, glm::vec2(15.0f, 0.0f), 9.8f);

	// create program
	GLuint shader_program = create_program_from_files("surface_vertex.glsl", "surface_fragment.glsl");

	// extract uniforms
	GLint model_location = glGetUniformLocation(shader_program, "model");
	GLint view_location = glGetUniformLocation(shader_program, "view");
	GLint projection_location = glGetUniformLocation(shader_program, "projection");

	// use shader program to get and push uniforms
	glUseProgram(shader_program);

	// prepare matrices
	glm::mat4 model = glm::mat4(1.0f);
	glm::mat4 view = glm::lookAt(glm::vec3(15.0f, 5.0f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f),
		glm::vec3(0.0f, 1.0f, 0.0f));
	glm::mat4 projection = glm::perspective(45.0f, 1.0f * g_width / g_height, 0.1f, 1000.0f);

	// and push it to shader program
	glUniformMatrix4fv(model_location, 1, GL_FALSE, &model[0][0]);
	glUniformMatrix4fv(view_location, 1, GL_FALSE, &view[0][0]);
	glUniformMatrix4fv(projection_location, 1, GL_FALSE, &projection[0][0]);

	//glClearColor(0.6f, 0.6f, 0.8f, 1.0f);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	// timer
	double old_time = glfwGetTime(), new_time;

	while (!glfwWindowShouldClose(g_window)) {
		update_fps_counter(g_window);

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glViewport(0, 0, g_width, g_height);

		new_time = glfwGetTime();

		/*-----render--------*/

		glUseProgram(shader_program);
		//ocean.updateOceanSlow((float)new_time);
		ocean.updateOcean((float)new_time);
		ocean.render();
		//printf("%f\n", new_time);

		/*-------------------*/

		glfwPollEvents();
		if (GLFW_PRESS == glfwGetKey(g_window, GLFW_KEY_ESCAPE)) {
			glfwSetWindowShouldClose(g_window, 1);
		}
		glfwSwapBuffers(g_window);
	}

	glfwTerminate();
	return 0;
}