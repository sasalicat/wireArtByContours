#version 330

uniform sampler2D ourTexture;

in vec2 TexCoord;
out vec4 FragColor;
uniform float zNear;
uniform float zFar;

float getRealDepth(float z_b){
	float z_n = 2.0*z_b - 1.0;
	float z_e = (z_b*(zNear - zFar)+2*zNear*zFar)/(zNear+zFar);
	float z_final = (z_e-zNear)/(zFar-zNear);
	return z_final;
}
void main()
{
	vec2 size = textureSize(ourTexture,0);
	float interval_x = 1/(size.x);
	float interval_y = 1/(size.y);
	vec2 offsets[4];//= {vec2(1,0),vec2(0,-1),vec2(-1.0),vec2(0,1)};
	offsets[0] = vec2(interval_x,0);
	offsets[1] = vec2(0,-interval_y);
	offsets[2] = vec2(-interval_x,0);
	offsets[3] = vec2(0,interval_y);
	
	vec4 Color = texture(ourTexture,TexCoord);
	float z_b = Color.x;
	float fc = getRealDepth(z_b);
	if(fc>=1){
		FragColor=vec4(0,0,1,1);
	}
	else if(fc<0.5 && fc>0){
		FragColor=vec4(0,1,0,1);
	}
	else if(fc<1 && fc>0.5){
		FragColor=vec4(fc,fc,fc,1);
	}
	else{
		FragColor=vec4(1,0,0,1);
	}
	
	/*
	if(z_b<1.0){
		for(int i=0;i<4;i++){
			vec2 pos = TexCoord+offsets[i];
			if(pos.x>=0 && pos.x<1 && pos.y>=0 && pos.y<1){
				if(texture(ourTexture,pos).x == 1.0){
					//FragColor = vec4(1.0);
					if(i==0){
						FragColor = vec4(1.0,0.0,0.0,1.0);
					}
					else if(i==1){
						FragColor = vec4(0.0,1.0,0.0,1.0);
					}
					else if(i==3){
						FragColor = vec4(0.0,0.0,1.0,1.0);
					}
					else{
						FragColor = vec4(0.5,0.5,0.5,1.0);
					}
					return;
				}
			}
		}
		FragColor = vec4(vec3(0.0),1.0);
	}
	else{
		FragColor = vec4(vec3(0.0),1.0);
	}
	*/
	//FragColor = vec4(TexCoord,0.0,1.0);
	//FragColor = vec4(size.x/10000,size.y/10000,0,1.0);
}
