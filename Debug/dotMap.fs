#version 330

in vec3 Normal0;
in vec3 vp_pos;
uniform vec3 light_direction;
layout (location = 0) out vec4 FragColor;

void main()
{
	vec3 toCam = normalize(vec3(0,0,0)-vp_pos);
	float strength = dot(normalize(Normal0),toCam);
	FragColor = vec4(strength,strength,strength,1);
}