#version 430

layout(binding = 0) uniform sampler1D g_texFresnel;
layout(binding = 1) uniform samplerCube g_texReflectCube;

in vec3 Position;
in vec3 Normal;
in vec3 WorldPos;

vec3 g_SkyColor = vec3(0.38f, 0.45f, 0.56f);
vec3 g_WaterbodyColor = vec3(0.07f, 0.15f, 0.2f);
//vec3 g_WaterbodyColor = vec3(0.11f, 0.42f, 0.63f);

vec3 g_SunDir = vec3(0.936016f, -0.343206f, 0.0780013f);
vec3 g_SunColor = vec3(1.0f, 1.0f, 0.6f);
float g_Shineness = 400.0f;

vec3 g_BendParam = vec3(0.1f, -0.4f, 0.2f);

vec3 g_WorldCameraPos = vec3(15.0f, 3.0f, 0.0f);

out vec4 fragment_color;

void main() {
	// --------------- Water body color
	vec3 eye_dir = normalize(g_WorldCameraPos - WorldPos);
	//vec3 eye_dir = normalize(-Position);
	// Calculate normal here.
	vec3 normal = normalize(Normal);
	// Reflected ray
	vec3 reflect_vec = reflect(-eye_dir, normal);
	// dot(N, V)
	float cos_angle = dot(normal, eye_dir);

	// A coarse way to handle transmitted light
	vec3 body_color = g_WaterbodyColor;

	// --------------- Reflected color

	// ramp.x for fresnel term. ramp.y for sky blending
	//!float4 ramp = g_texFresnel.Sample(g_samplerFresnel, cos_angle).xyzw;
	vec4 ramp = texture(g_texFresnel, cos_angle).xyzw;
	// A workaround to deal with "indirect reflection vectors" (which are rays requiring multiple
	// reflections to reach the sky).
	if (reflect_vec.z < g_BendParam.x) {
		ramp = mix(ramp, vec4(g_BendParam, 1.0), (g_BendParam.x - reflect_vec.z) / (g_BendParam.x - g_BendParam.y));
		//ramp.x = mix(ramp.x, g_BendParam.z, (g_BendParam.x - reflect_vec.z) / (g_BendParam.x - g_BendParam.y));
	}
	reflect_vec.z = max(0, reflect_vec.z);

	//!float3 reflection = g_texReflectCube.Sample(g_samplerCube, reflect_vec).xyz;
	vec3 reflection = texture(g_texReflectCube, reflect_vec).xyz;
	// Hack bit: making higher contrast
	reflection = reflection * reflection * 2.5f;

	// Blend with predefined sky color
	vec3 reflected_color = mix(g_SkyColor, reflection, ramp.y);

	// Combine waterbody color and reflected color
	vec3 water_color = mix(body_color, reflected_color, ramp.x);

	// --------------- Sun spots
	float cos_spec = clamp(dot(reflect_vec, g_SunDir), 0, 1);
	float sun_spot = pow(cos_spec, g_Shineness);
	water_color += g_SunColor * sun_spot;

	fragment_color = vec4(water_color, 1.0);
}