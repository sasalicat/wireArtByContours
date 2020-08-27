#version 330

uniform sampler2D ourTexture;
uniform sampler2D depthTexture;

in vec2 TexCoord;
out vec4 FragColor;
uniform float zNear;
uniform float zFar;
float threshold = 0.2;


void main()
{

	//if(texture(ourTexture,TexCoord).x)
	
	vec4 Color = texture(ourTexture,TexCoord);
	//FragColor = Color;
	vec3 real_normal =(Color.xyz-vec3(0.5,0.5,0.5))*2;
	
	if(abs(real_normal.z)<threshold){
		FragColor = vec4(1.0);
		return;
	}
	else{
		FragColor = vec4(vec3(0.0),1.0);
	}
	
	//增加深度测试来补齐边
	vec2 size = textureSize(depthTexture,0);
	float interval_x = 1/(size.x);
	float interval_y = 1/(size.y);
	vec2 offsets[4];//= {vec2(1,0),vec2(0,-1),vec2(-1.0),vec2(0,1)};
	offsets[0] = vec2(interval_x,0);
	offsets[1] = vec2(0,-interval_y);
	offsets[2] = vec2(-interval_x,0);
	offsets[3] = vec2(0,interval_y);
	
	vec4 dColor = texture(depthTexture,TexCoord);
	float z_b = dColor.x;
	
	if(z_b<1.0){
		for(int i=0;i<4;i++){
			vec2 pos = TexCoord+offsets[i];
			if(pos.x>=0 && pos.x<1 && pos.y>=0 && pos.y<1){
				if(texture(depthTexture,pos).x == 1.0){
					//FragColor = vec4(1.0);
					FragColor = vec4(1.0);
					return;
				}
			}
		}
		FragColor = vec4(vec3(0.0),1.0);
	}
	else{
		FragColor = vec4(vec3(0.0),1.0);
	}
	//FragColor= vec4(z_b,0,0,1);
}