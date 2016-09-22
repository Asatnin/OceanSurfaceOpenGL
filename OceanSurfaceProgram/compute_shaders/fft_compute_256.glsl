#version 430

// compute shader for FFT

#define N 256 // N should be power of 2
#define M_PI 3.14159265358979323846

uniform int fft_column; // flag for processing FFT in rows or columns: 0 or 1
uniform int steps; // steps for FFT: steps = log2(N)

layout(rg32f, binding = 0) readonly uniform image2D imageIn; // FFT input
layout(rg32f, binding = 1) writeonly uniform image2D imageOut; // FFT output

layout(std430, binding = 2) buffer Ind {
	int indices[];
} bin_reversed; // binary reversed indices

shared vec2 sharedStore[N]; // here we write buttefly operation result on each step

layout(local_size_x = 128, local_size_y = 1) in; // local_size_x should be N / 2

// complex number multiplication
vec2 multiply_complex(vec2 a, vec2 b) {
	vec2 res;

	res.x = a.x * b.x - a.y * b.y; // real part
	res.y = a.x * b.y + a.y * b.x; // imaginary part

	return res;
}

// return exp(i * 2 * PI * k / n);
vec2 twiddleFactor(int n, int k) {
	vec2 res;

	float arg = 2.0 * M_PI * float(k) / float(n);

	res.x = cos(arg);
	res.y = sin(arg);

	return res;
}

void main() {
	// gloval values of X and Y in computation grid
	int x = int(gl_GlobalInvocationID.x);
	int y = int(gl_GlobalInvocationID.y);

	// we process two elements in one invocation: left element and right element
	int leftStoreIndex = 2 * x;
	int rightStoreIndex = 2 * x + 1;

	// extract binary reversed indices
	int leftLoadIndex = bin_reversed.indices[leftStoreIndex];
	int rightLoadIndex = bin_reversed.indices[rightStoreIndex];

	// position of left and right elemetn in FFT output imageOut
	ivec2 leftStorePos;
	ivec2 rightStorePos;

	// position of left and right elemetn in FFT input imageIn
	ivec2 leftLoadPos;
	ivec2 rightLoadPos;

	// make input and output positions for left and right elements due to processing FFT in rows or columns
	if (fft_column == 0) {
		// for rows
		leftLoadPos = ivec2(leftLoadIndex, y);
		rightLoadPos = ivec2(rightLoadIndex, y);

		leftStorePos = ivec2(leftStoreIndex, y);
		rightStorePos = ivec2(rightStoreIndex, y);
	}
	else {
		// for columns
		leftLoadPos = ivec2(y, leftLoadIndex);
		rightLoadPos = ivec2(y, rightLoadIndex);

		leftStorePos = ivec2(y, leftStoreIndex);
		rightStorePos = ivec2(y, rightStoreIndex);
	}

	// extract FFT input complex values
	vec2 leftValue = imageLoad(imageIn, leftLoadPos).xy;
	vec2 rightValue = imageLoad(imageIn, rightLoadPos).xy;

	// write to shared store to prepare for 'steps' passes
	sharedStore[leftStoreIndex] = leftValue;
	sharedStore[rightStoreIndex] = rightValue;

	// synchronization threads in work group
	memoryBarrierShared();
	barrier();

	// now all threads in work group are here

	int numberGroups = N / 2; // for N = 8: 4 groups with two elements in first pass
	int butterfliesInGroup = 1; // for N = 8: 1 butterfly operation with two elements in first pass
	// numberGroups: 4, 2, 1 for N = 8
	// butterfliesInGroup: 1, 2, 4 for N = 8

	int currentGroup = x;
	int currentButterfly = 0;

	// process steps = log2(N) passes
	for (int step = 1; step <= steps; step++) {
		int leftIndex = currentButterfly + currentGroup * butterfliesInGroup * 2;
		int rightIndex = currentButterfly + butterfliesInGroup + currentGroup * butterfliesInGroup * 2;

		leftValue = sharedStore[leftIndex];
		rightValue = sharedStore[rightIndex];

		// process butterfly operation
		vec2 twiddle = twiddleFactor(2 * butterfliesInGroup, currentButterfly);

		vec2 t = multiply_complex(twiddle, rightValue);
		
		sharedStore[leftIndex] = leftValue + t;
		sharedStore[rightIndex] = leftValue - t;

		// synchronize shared memory operations by threads in work group
		memoryBarrierShared();

		numberGroups /= 2;
		butterfliesInGroup *= 2;

		currentGroup /= 2;
		currentButterfly = x % butterfliesInGroup;

		barrier();

		// now all threads in work group are here
	}

	// make sign flipping after columns passes
	if (fft_column == 1) {
		if ((leftStorePos.x + leftStorePos.y) % 2 == 1) {
			sharedStore[leftStoreIndex] *= -1.0;
		}
		if ((rightStorePos.x + rightStorePos.y) % 2 == 1) {
			sharedStore[rightStoreIndex] *= -1.0;
		}
		// synchronize shared memory operations by threads in work group
		memoryBarrierShared();
	}

	// write FFT output
	imageStore(imageOut, leftStorePos, vec4(sharedStore[leftStoreIndex], 0.0, 0.0));
	imageStore(imageOut, rightStorePos, vec4(sharedStore[rightStoreIndex], 0.0, 0.0));
}