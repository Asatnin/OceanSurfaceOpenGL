#version 430

layout(location = 0) in vec3 vertex_position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec2 img_coord;

//layout(rg32f, binding = 3) readonly uniform image2D imageFFT;
layout(binding = 0) uniform sampler2D imageFFT;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main() {
	//vec2 height = imageLoad(imageFFT, ivec2(img_coord.x, img_coord.y)).xy;
	vec2 height = texture(imageFFT, vec2(0.5, 0.5)).xy;
	vec3 disp = vec3(0.0, 0.0, 0.0);
	if (height.x != 0.0 || height.y != 0.0) {
		disp = vec3(0.0, 5.0, 0.0);
	}

	gl_Position = projection * view * model * vec4(vertex_position + disp, 1.0);
}