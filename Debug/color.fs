#version 330

uniform sampler2D ourTexture;

in vec2 TexCoord;
out vec4 FragColor;
uniform float zNear;
uniform float zFar;

void main()
{

	vec4 Color = texture(ourTexture,TexCoord);
	FragColor = Color;
}