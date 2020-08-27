#version 330

uniform sampler2D ourTexture;

in vec2 TexCoord;
out vec4 FragColor;
uniform float zNear;
uniform float zFar;

void main()
{

	//vec4 Color = texture(ourTexture,TexCoord);
	//float z_b = Color.x;
	//float z_n = 2.0*z_b - 1.0;
	//float z_e = (z_n*(zNear-zFar)+2*zNear*zFar)/(zNear+zFar);
	//float z_final = (z_e - zNear)/(zFar - zNear);
	FragColor = texture(ourTexture,TexCoord);
}