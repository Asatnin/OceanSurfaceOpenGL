#version 430

in vec3 Position;
in vec3 Normal;

vec3 LightPosition = vec3(1.0, 1.0, 1.0); // vector towards directional light source
//vec3 LightIntensity = vec3(0.3, 0.3, 1.8);
vec3 LightIntensity = vec3(0.07, 0.15, 0.2);
vec3 Kd = vec3(0.5, 0.65, 0.75); // diffuse reflectivity
vec3 Ka = vec3(0.0, 0.65, 0.75); // ambient reflectivity
vec3 Ks = vec3(1.0, 0.25, 0.0); // specular reflectivity
float Shininess = 120.0; // specular shineness factor

out vec4 fragment_color;

vec3 ads() {
	vec3 n = normalize(Normal);
	vec3 s = normalize(LightPosition);
	vec3 v = normalize(-Position);
	vec3 h = normalize(v + s);

	float sDotN = max(dot(s, n), 0.0);
	vec3 spec = vec3(0.0);
	if (sDotN > 0.0) {
		spec = Ks * pow(max(dot(h, n), 0.0), Shininess);
	}

	return LightIntensity * (Ka + Kd * sDotN + spec);
}

void main() {
	//fragment_color = vec4(0.11, 0.42, 0.63, 1.0);
	fragment_color = vec4(ads(), 1.0);
}