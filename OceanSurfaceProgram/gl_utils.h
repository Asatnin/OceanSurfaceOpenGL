#pragma once

#include <GL\glew.h>
#include <GLFW\glfw3.h>
#include <cstdarg>

#define GL_LOG_FILE "gl.log"

extern int g_width;
extern int g_height;
extern GLFWwindow *g_window;

/*--------------------------------LOG FUNCTIONS------------------------------*/
bool restart_gl_log();
bool gl_log(const char *message, ...);
bool gl_log_err(const char *message, ...);
void log_gl_params();

/*--------------------------------GLFW3 and GLEW-----------------------------*/
bool start_gl();
void glfw_error_callback(int error, const char *description);
void glfw_window_size_callback(GLFWwindow *window, int width, int height);
void update_fps_counter(GLFWwindow *window);

/*--------------------------------Shaders and Programs-----------------------------*/
void print_shader_info_log(GLuint shader_index); // call it after unsuccessful compiling, for example
void print_program_info_log(GLuint program_index); // call it after unsuccessful linking, for example
const char* GL_type_to_string(GLenum type); // get string description of GLenum
void print_all(GLuint program); // print all information about program: attributes, uniforms, etc.
bool is_program_valid(GLuint program);

// next functions should be in ShaderManager, for example
bool parse_file_into_str(const char* file_name, char* shader_str, int max_len);
bool create_shader(const char* file_name, GLuint* shader, GLenum type);
bool create_program(GLuint vert, GLuint frag, GLuint* programme);
GLuint create_program_from_files(const char* vert_file_name, const char* frag_file_name);