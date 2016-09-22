#include "gl_utils.h"
#include <cstdio>
#include <ctime>
#include <cstring>
#include <cassert>

#define MAX_SHADER_LENGTH 262144

/*--------------------------------LOG FUNCTIONS------------------------------*/
bool restart_gl_log() {
	time_t now;
	char *date;
	FILE *file = NULL;
	fopen_s(&file, GL_LOG_FILE, "w");

	if (!file) {
		fprintf(stderr, "ERROR: could not open GL_LOG_FILE log file %s for writing\n", GL_LOG_FILE);
		return false;
	}
	now = time(NULL);
	date = ctime(&now);
	fprintf(file, "GL_LOG_FILE log. local time %s\n", date);
	fclose(file);
	return true;
}

bool gl_log(const char *message, ...) {
	va_list argptr;
	FILE *file = NULL;
	fopen_s(&file, GL_LOG_FILE, "a");

	if (!file) {
		fprintf(stderr, "ERROR: could not open GL_LOG_FILE %s file for appending\n", GL_LOG_FILE);
		return false;
	}
	va_start(argptr, message);
	vfprintf(file, message, argptr);
	va_end(argptr);
	fclose(file);
	return true;
}

bool gl_log_err(const char *message, ...) {
	va_list argptr;
	FILE *file = NULL;
	fopen_s(&file, GL_LOG_FILE, "a");

	if (!file) {
		fprintf(stderr, "ERROR: could not open GL_LOG_FILE %s file for appending\n", GL_LOG_FILE);
		return false;
	}
	va_start(argptr, message);
	vfprintf(file, message, argptr);
	va_end(argptr);
	va_start(argptr, message);
	vfprintf(stderr, message, argptr);
	va_end(argptr);
	fclose(file);
	return true;
}

void log_gl_params() {
	int i;
	int v[2];
	unsigned char s = 0;
	GLenum params[] = {
		GL_MAX_VERTEX_IMAGE_UNIFORMS,
		GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS,
		GL_MAX_CUBE_MAP_TEXTURE_SIZE,
		GL_MAX_DRAW_BUFFERS,
		GL_MAX_FRAGMENT_UNIFORM_COMPONENTS,
		GL_MAX_TEXTURE_IMAGE_UNITS,
		GL_MAX_TEXTURE_SIZE,
		GL_MAX_VARYING_FLOATS,
		GL_MAX_VERTEX_ATTRIBS,
		GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS,
		GL_MAX_VERTEX_UNIFORM_COMPONENTS,
		GL_MAX_VIEWPORT_DIMS,
		GL_STEREO,
	};
	const char* names[] = {
		"GL_MAX_VERTEX_IMAGE_UNIFORMS",
		"GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS",
		"GL_MAX_CUBE_MAP_TEXTURE_SIZE",
		"GL_MAX_DRAW_BUFFERS",
		"GL_MAX_FRAGMENT_UNIFORM_COMPONENTS",
		"GL_MAX_TEXTURE_IMAGE_UNITS",
		"GL_MAX_TEXTURE_SIZE",
		"GL_MAX_VARYING_FLOATS",
		"GL_MAX_VERTEX_ATTRIBS",
		"GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS",
		"GL_MAX_VERTEX_UNIFORM_COMPONENTS",
		"GL_MAX_VIEWPORT_DIMS",
		"GL_STEREO",
	};
	gl_log("GL Context Params:\n");
	// integers - only works if the order is 0-10 integer return types
	for (i = 0; i < 11; i++) {
		int v = 0;
		glGetIntegerv(params[i], &v);
		gl_log("%s %i\n", names[i], v);
	}
	// others
	v[0] = v[1] = 0;
	glGetIntegerv(params[11], v);
	gl_log("%s %i %i\n", names[10], v[0], v[1]);
	glGetBooleanv(params[12], &s);
	gl_log("%s %i\n", names[12], (unsigned int)s);
	gl_log("-----------------------------\n");
}

/*--------------------------------GLFW3 and GLEW-----------------------------*/
bool start_gl() {
	const GLubyte *renderer;
	const GLubyte *version;

	gl_log("starting GLFW %s\n", glfwGetVersionString());

	glfwSetErrorCallback(glfw_error_callback);
	if (!glfwInit()) {
		fprintf(stderr, "ERROR: could not start GLFW3\n");
		return false;
	}

	// Seting OpenGL version hints
	/*
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	*/

	g_window = glfwCreateWindow(g_width, g_height, "CG", NULL, NULL);
	if (!g_window) {
		fprintf(stderr, "ERROR: could not open window with GLFW3\n");
		glfwTerminate();
		return false;
	}

	glfwSetWindowSizeCallback(g_window, glfw_window_size_callback);
	glfwMakeContextCurrent(g_window);

	glfwWindowHint(GLFW_SAMPLES, 4);

	/* start GLEW extension handler */
	glewExperimental = GL_TRUE;
	glewInit();

	renderer = glGetString(GL_RENDERER);
	version = glGetString(GL_VERSION);
	printf("Renderer: %s\n", renderer);
	printf("OpenGL version supported: %s\n", version);
	gl_log("renderer: %s\nversion: %s\n", renderer, version);
	log_gl_params();

	return true;
}

void glfw_error_callback(int error, const char *description) {
	fputs(description, stderr);
	gl_log_err("GLFW ERROR: code %d msg: %s\n", error, description);
}

void glfw_window_size_callback(GLFWwindow *window, int width, int height) {
	g_width = width;
	g_height = height;
	printf("width: %d height: %d\n", width, height);
}

void update_fps_counter(GLFWwindow *window) {
	static double previous_seconds = glfwGetTime();
	static int frame_count;
	double current_seconds = glfwGetTime();
	double elapsed_seconds = current_seconds - previous_seconds;

	if (elapsed_seconds > 0.25) {
		previous_seconds = current_seconds;
		double fps = (double)frame_count / elapsed_seconds;
		char tmp[128];
		sprintf(tmp, "opengl @ fps: %.2f", fps);
		glfwSetWindowTitle(window, tmp);
		frame_count = 0;
	}
	frame_count++;
}

/*--------------------------------Shaders and Programs-----------------------------*/
void print_shader_info_log(GLuint shader_index) {
	int max_length = 2048, actual_length = 0;
	char log[2048];
	glGetShaderInfoLog(shader_index, max_length, &actual_length, log);
	printf("shader info log for GL shader index %d:\n%s\n", shader_index, log);
	gl_log("shader info log for GL shader index %d:\n%s\n", shader_index, log);
}

void print_program_info_log(GLuint program_index) {
	int max_length = 2048, actual_length = 0;
	char log[2048];
	glGetProgramInfoLog(program_index, max_length, &actual_length, log);
	printf("program info log for GL program index %d:\n%s\n", program_index, log);
	gl_log("program info log for GL program index %d:\n%s\n", program_index, log);
}

const char* GL_type_to_string(GLenum type) {
	switch (type) {
	case GL_BOOL: return "bool";
	case GL_INT: return "int";
	case GL_FLOAT: return "float";
	case GL_FLOAT_VEC2: return "vec2";
	case GL_FLOAT_VEC3: return "vec3";
	case GL_FLOAT_VEC4: return "vec4";
	case GL_FLOAT_MAT2: return "mat2";
	case GL_FLOAT_MAT3: return "mat3";
	case GL_FLOAT_MAT4: return "mat4";
	case GL_SAMPLER_2D: return "sampler2D";
	case GL_SAMPLER_3D: return "sampler3D";
	case GL_SAMPLER_CUBE: return "samplerCube";
	case GL_SAMPLER_2D_SHADOW: return "sampler2DShadow";
	default: return "other";
	}
	return "other";
}

void print_all(GLuint program) {
	printf("-------------\nshader program %d info:\n", program);
	gl_log("-------------\nshader program %d info:\n", program);
	int params = -1;

	glGetProgramiv(program, GL_LINK_STATUS, &params);
	printf("GL_LINK_STATUS = %d\n", params);
	gl_log("GL_LINK_STATUS = %d\n", params);

	glGetProgramiv(program, GL_ATTACHED_SHADERS, &params);
	printf("GL_ATTACHED_SHADERS = %d\n", params);
	gl_log("GL_ATTACHED_SHADERS = %d\n", params);

	// attributes info
	glGetProgramiv(program, GL_ACTIVE_ATTRIBUTES, &params);
	printf("GL_ACTIVE_ATTRIBUTES = %d\n", params);
	gl_log("GL_ACTIVE_ATTRIBUTES = %d\n", params);
	for (GLuint i = 0; i < (GLuint)params; i++) {
		char name[64];
		int max_length = 64, actual_length = 0, size = 0;
		GLenum type;
		glGetActiveAttrib(program, i, max_length, &actual_length, &size, &type, name);
		if (size > 1) {
			for (int j = 0; j < size; j++) {
				char long_name[64];
				sprintf(long_name, "%s[%d]", name, j);
				int location = glGetAttribLocation(program, long_name);
				printf("   %d) type: %s name: %s location: %d\n", i, GL_type_to_string(type), long_name, location);
				gl_log("   %d) type: %s name: %s location: %d\n", i, GL_type_to_string(type), long_name, location);
			}
		}
		else {
			int location = glGetAttribLocation(program, name);
			printf("   %d) type: %s name: %s location: %d\n", i, GL_type_to_string(type), name, location);
			gl_log("   %d) type: %s name: %s location: %d\n", i, GL_type_to_string(type), name, location);
		}
	}

	// uniforms info
	glGetProgramiv(program, GL_ACTIVE_UNIFORMS, &params);
	printf("GL_ACTIVE_UNIFORMS = %d\n", params);
	gl_log("GL_ACTIVE_UNIFORMS = %d\n", params);
	for (GLuint i = 0; i < (GLuint)params; i++) {
		char name[64];
		int max_length = 64, actual_length = 0, size = 0;
		GLenum type;
		glGetActiveUniform(program, i, max_length, &actual_length, &size, &type, name);
		if (size > 1) {
			for (int j = 0; j < size; j++) {
				char long_name[64];
				sprintf(long_name, "%s[%d]", name, j);
				int location = glGetUniformLocation(program, long_name);
				printf("   %d) type: %s name: %s location: %d\n", i, GL_type_to_string(type), long_name, location);
				gl_log("   %d) type: %s name: %s location: %d\n", i, GL_type_to_string(type), long_name, location);
			}
		}
		else {
			int location = glGetUniformLocation(program, name);
			printf("   %d) type: %s name: %s location: %d\n", i, GL_type_to_string(type), name, location);
			gl_log("   %d) type: %s name: %s location: %d\n", i, GL_type_to_string(type), name, location);
		}
	}

	print_program_info_log(program);
}

bool is_program_valid(GLuint program) {
	glValidateProgram(program);
	int params = -1;
	glGetProgramiv(program, GL_VALIDATE_STATUS, &params);
	if (params != GL_TRUE) {
		gl_log_err("program %i GL_VALIDATE_STATUS = GL_FALSE\n", program);
		print_program_info_log(program);
		return false;
	}
	gl_log("program %i GL_VALIDATE_STATUS = GL_TRUE\n", program);
	return true;
}

// shader manager functions
bool parse_file_into_str(const char* file_name, char* shader_str, int max_len) {
	shader_str[0] = '\0'; // reset string
	FILE* file = fopen(file_name, "r");
	if (!file) {
		gl_log_err("ERROR: opening file for reading: %s\n", file_name);
		return false;
	}
	int current_len = 0;
	char line[2048];
	strcpy(line, ""); // remember to clean up before using for first time!
	while (!feof(file)) {
		if (fgets(line, 2048, file) != NULL) {
			current_len += strlen(line); // +1 for \n at end
			if (current_len >= max_len) {
				gl_log_err("ERROR: shader length is longer than string buffer length %i\n", max_len);
			}
			strcat(shader_str, line);
		}
	}
	if (fclose(file) == EOF) { // probably unnecesssary validation
		gl_log_err("ERROR: closing file from reading %s\n", file_name);
		return false;
	}
	return true;
}

bool create_shader(const char* file_name, GLuint* shader, GLenum type) {
	gl_log("creating shader from %s...\n", file_name);
	char shader_string[MAX_SHADER_LENGTH];
	assert(parse_file_into_str(file_name, shader_string, MAX_SHADER_LENGTH));
	*shader = glCreateShader(type);
	const GLchar* p = (const GLchar*)shader_string;
	glShaderSource(*shader, 1, &p, NULL);
	glCompileShader(*shader);
	// check for compile errors
	int params = -1;
	glGetShaderiv(*shader, GL_COMPILE_STATUS, &params);
	if (GL_TRUE != params) {
		gl_log_err("ERROR: GL shader index %i did not compile\n", *shader);
		print_shader_info_log(*shader);
		return false; // or exit or something
	}
	gl_log("shader compiled. index %i\n", *shader);
	return true;
}

bool create_program(GLuint comp, GLuint* programme) {
	*programme = glCreateProgram();
	gl_log("created programme %u. attaching compute shader %u...\n", *programme, comp);
	glAttachShader(*programme, comp);
	// link the shader programme. if binding input attributes do that before link
	glLinkProgram(*programme);
	GLint params = -1;
	glGetProgramiv(*programme, GL_LINK_STATUS, &params);
	if (GL_TRUE != params) {
		gl_log_err("ERROR: could not link shader programme GL index %u\n", *programme);
		print_program_info_log(*programme);
		return false;
	}
	assert(is_program_valid(*programme));
	// delete shaders here to free memory
	glDeleteShader(comp);
	return true;
}

bool create_program(GLuint vert, GLuint geom, GLuint frag, GLuint* programme) {
	*programme = glCreateProgram();
	gl_log("created programme %u. attaching shaders %u, %u and %u...\n", *programme, vert, geom, frag);
	glAttachShader(*programme, vert);
	glAttachShader(*programme, geom);
	glAttachShader(*programme, frag);
	// link the shader programme. if binding input attributes do that before link
	glLinkProgram(*programme);
	GLint params = -1;
	glGetProgramiv(*programme, GL_LINK_STATUS, &params);
	if (GL_TRUE != params) {
		gl_log_err("ERROR: could not link shader programme GL index %u\n", *programme);
		print_program_info_log(*programme);
		return false;
	}
	assert(is_program_valid(*programme));
	// delete shaders here to free memory
	glDeleteShader(vert);
	glDeleteShader(geom);
	glDeleteShader(frag);
	return true;
}

bool create_program(GLuint vert, GLuint frag, GLuint* programme) {
	*programme = glCreateProgram();
	gl_log("created programme %u. attaching shaders %u and %u...\n", *programme, vert, frag);
	glAttachShader(*programme, vert);
	glAttachShader(*programme, frag);
	// link the shader programme. if binding input attributes do that before link
	glLinkProgram(*programme);
	GLint params = -1;
	glGetProgramiv(*programme, GL_LINK_STATUS, &params);
	if (GL_TRUE != params) {
		gl_log_err("ERROR: could not link shader programme GL index %u\n", *programme);
		print_program_info_log(*programme);
		return false;
	}
	assert(is_program_valid(*programme));
	// delete shaders here to free memory
	glDeleteShader(vert);
	glDeleteShader(frag);
	return true;
}

GLuint create_comp_program_from_file(const char* comp_file_name) {
	GLuint comp, programme;
	assert(create_shader(comp_file_name, &comp, GL_COMPUTE_SHADER));
	assert(create_program(comp, &programme));
	return programme;
}

GLuint create_program_from_files(const char* vert_file_name, const char* frag_file_name) {
	GLuint vert, frag, programme;
	assert(create_shader(vert_file_name, &vert, GL_VERTEX_SHADER));
	assert(create_shader(frag_file_name, &frag, GL_FRAGMENT_SHADER));
	assert(create_program(vert, frag, &programme));
	return programme;
}

GLuint create_program_from_files(const char* vert_file_name, const char* geom_file_name, const char* frag_file_name) {
	GLuint vert, geom, frag, programme;
	assert(create_shader(vert_file_name, &vert, GL_VERTEX_SHADER));
	assert(create_shader(geom_file_name, &geom, GL_GEOMETRY_SHADER));
	assert(create_shader(frag_file_name, &frag, GL_FRAGMENT_SHADER));
	assert(create_program(vert, geom, frag, &programme));
	return programme;
}