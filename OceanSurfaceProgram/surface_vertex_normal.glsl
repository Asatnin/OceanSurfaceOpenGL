#version 430

layout(location = 0) in vec3 vertex_position;
layout(location = 1) in vec3 normal;

out vec4 out_normal;

void main() {
	gl_Position = vec4(vertex_position, 1.0);
	out_normal = vec4(normal, 0.0);
}