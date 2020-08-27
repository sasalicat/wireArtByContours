#version 330

in vec2 UV;

out vec4 FragColor;

uniform sampler2D myTextureSampler;

void main()
{
    FragColor = vec4(texture( myTextureSampler, UV ).rgb,1.0f); //Color;
}
