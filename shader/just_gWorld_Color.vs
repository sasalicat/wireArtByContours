#version 330

layout (location = 0) in vec3 Position;
layout (location = 1) in vec4 Color;
out vec4 color;
uniform mat4 gWorld;


void main()
{
	color = Color;
    gl_Position = gWorld * vec4(Position, 1.0);
}
