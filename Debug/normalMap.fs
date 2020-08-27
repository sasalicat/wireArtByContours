#version 330

in vec3 Normal0;
uniform vec3 light_direction;
layout (location = 0) out vec4 FragColor;

void main()
{
	//vec3 light_color = vec3(1.0, 1.0, 1.0);
	//vec3 ambient = vec3(0.2,0.2,0.2);
	//FragColor = vec4(max(0,dot(light_direction,Normal0))*light_color + ambient,1);
	//FragColor = vec4(normalize(Normal0),1.0)*0.5+vec4(0.5,0.5,0.5,0.5);
	FragColor = vec4(Normal0,1.0);
}