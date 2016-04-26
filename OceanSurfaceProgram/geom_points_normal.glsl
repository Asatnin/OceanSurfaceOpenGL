#version 430
layout(lines) in;
layout(line_strip, max_vertices = 6) out;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

in vec4 out_normal[];

out vec4 vertex_color;

void main() {
	for (int i = 0; i < gl_in.length(); i++) {
		vec4 init_point = gl_in[i].gl_Position;

		gl_Position = projection * view * model * init_point;
		vertex_color = vec4(0.11, 0.42, 0.63, 1.0);
		EmitVertex();
	}

	EndPrimitive();

	for (int i = 0; i < gl_in.length(); i++) {
		vec4 init_point = gl_in[i].gl_Position;

		gl_Position = projection * view * model * init_point;
		vertex_color = vec4(1.0, 0.0, 0.0, 0.0);
		EmitVertex();

		gl_Position = projection * view * model * (init_point + out_normal[i] * 0.5);
		vertex_color = vec4(1.0, 0.0, 0.0, 0.0);
		EmitVertex();

		EndPrimitive();
	}
}