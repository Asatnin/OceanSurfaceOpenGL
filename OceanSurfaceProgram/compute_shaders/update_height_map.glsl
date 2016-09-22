#version 430

#define M_PI 3.14159265358979323846
#define G 9.81

layout(rg32f, binding = 0) readonly uniform image2D imageIn;
layout(rg32f, binding = 1) readonly uniform image2D imageIn_CC;
layout(rg32f, binding = 2) writeonly uniform image2D imageOut;
layout(rg32f, binding = 3) writeonly uniform image2D dxOut;
layout(rg32f, binding = 4) writeonly uniform image2D dzOut;
layout(rg32f, binding = 5) writeonly uniform image2D gradxOut;
layout(rg32f, binding = 6) writeonly uniform image2D gradzOut;

uniform float totalTime;
uniform float L; // real length of surface
uniform int N; // size of surface

layout(local_size_x = 1, local_size_y = 1) in;

// complex number multiplication
vec2 multiply_complex(vec2 a, vec2 b) {
	vec2 res;

	res.x = a.x * b.x - a.y * b.y; // real part
	res.y = a.x * b.y + a.y * b.x; // imaginary part

	return res;
}

// dispersion relation
float dispersion_relation(int n, int m) {
	float w_0 = 2.0f * M_PI / 200.0f;
	vec2 k = vec2(M_PI * (2.0f * n - N) / L, M_PI * (2.0f * m - N) / L); // k is wavevector	
	//return sqrt(g * glm::length(k));
	return floor(sqrt(G * length(k)) / w_0) * w_0;
}

void main() {
	// gloval values of X and Y in computation grid
	int x = int(gl_GlobalInvocationID.x);
	int y = int(gl_GlobalInvocationID.y);

	ivec2 storePos = ivec2(x, y);

	vec2 h_t0 = imageLoad(imageIn, storePos).xy;
	vec2 h_t0_cc = imageLoad(imageIn_CC, storePos).xy;

	// coordinates of wavevector
	float k_x = M_PI * (2.0f * x - N) / L;
	float k_z = M_PI * (2.0f * y - N) / L;
	vec2 k = vec2(k_x, k_z); // k is wavevector
	float k_len = length(k);

	// dispersion relation
	float w_0 = 2.0f * M_PI / 200.0f;
	//float w_freq = dispersion_relation(x, y);
	float w_freq = floor(sqrt(G * k_len) / w_0) * w_0;
	float wt = w_freq * totalTime;

	vec2 c1 = vec2(cos(wt), sin(wt));
	vec2 c2 = vec2(cos(wt), -sin(wt));

	vec2 p1 = multiply_complex(h_t0, c1);
	vec2 p2 = multiply_complex(h_t0_cc, c2);

	// complex number addition
	vec2 res = vec2(p1.x + p2.x, p1.y + p2.y);

	// store height map
	imageStore(imageOut, storePos, vec4(res, 0.0, 0.0));

	// store displacement map
	if (k_len < 0.000001) { // ignore small wavevectors, put (0, 0) displacement
		imageStore(dxOut, storePos, vec4(0.0, 0.0, 0.0, 0.0));
		imageStore(dzOut, storePos, vec4(0.0, 0.0, 0.0, 0.0));
	}
	else { // calculate displacement via formula
		vec2 disp = multiply_complex(vec2(0.0, -k_x / k_len), res);
		imageStore(dxOut, storePos, vec4(disp, 0.0, 0.0));

		disp = multiply_complex(vec2(0.0, -k_z / k_len), res);
		imageStore(dzOut, storePos, vec4(disp, 0.0, 0.0));
	}

	// store gradient map
	vec2 grad = multiply_complex(vec2(0.0, k_x), res);
	imageStore(gradxOut, storePos, vec4(grad, 0.0, 0.0));

	grad = multiply_complex(vec2(0.0, k_z), res);
	imageStore(gradzOut, storePos, vec4(grad, 0.0, 0.0));
}