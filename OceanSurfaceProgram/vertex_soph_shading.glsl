#version 430

layout(location = 0) in vec3 vertex_position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec2 img_coord;

layout(rg32f, binding = 0) readonly uniform image2D imageFFT;
layout(rg32f, binding = 1) readonly uniform image2D imageDx;
layout(rg32f, binding = 2) readonly uniform image2D imageDz;
layout(rg32f, binding = 3) readonly uniform image2D imageGradX;
layout(rg32f, binding = 4) readonly uniform image2D imageGradZ;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec3 Position;
out vec3 Normal;
out vec3 WorldPos;

void main() {
	vec2 height = imageLoad(imageFFT, ivec2(img_coord.x, img_coord.y)).xy;
	vec2 dx = imageLoad(imageDx, ivec2(img_coord.x, img_coord.y)).xy;
	vec2 dz = imageLoad(imageDz, ivec2(img_coord.x, img_coord.y)).xy;
	//vec3 disp = vec3(0.0, 0.0, 0.0);
	vec3 disp = vec3(-dx.x, height.x, -dz.x);
	//vec3 disp = vec3(0.0, height.x, 0.0);

	// normal
	vec2 grad_x = imageLoad(imageGradX, ivec2(img_coord.x, img_coord.y)).xy;
	vec2 grad_z = imageLoad(imageGradZ, ivec2(img_coord.x, img_coord.y)).xy;
	vec3 normal = normalize(vec3(-grad_x.x, 1.0, -grad_z.x));

	//Normal = (inverse(transpose(view * model)) * vec4(normal, 0.0)).xyz;
	Normal = (inverse(transpose(model)) * vec4(normal, 0.0)).xyz;
	//Normal = (model * vec4(normal, 0.0)).xyz;
	Position = vec3(view * model * vec4(vertex_position + disp, 1.0));
	WorldPos = vec3(model * vec4(vertex_position + disp, 1.0));
	gl_Position = projection * view * model * vec4(vertex_position + disp, 1.0);
}