#version 330

uniform sampler2D ourTexture;

in vec2 TexCoord;
out vec4 FragColor;
uniform float zNear;
uniform float zFar;

void main()
{

	vec4 Color = texture(ourTexture,TexCoord);
	//FragColor = Color;
	vec3 real_normal =(Color.xyz-vec3(0.5,0.5,0.5))*2;
	
	if(abs(real_normal.z)<0.1){
		FragColor = vec4(1.0);
	}
	else{
		FragColor = vec4(vec3(0.0),1.0);
	}
	
	/*
	
	if(FragColor.x<0){
		FragColor = vec4(1.0,0.0,0.0,1.0);
	}
	else{
		FragColor = vec4(color.x);
	}
	*/
}