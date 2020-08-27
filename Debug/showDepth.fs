#version 330

uniform sampler2D depthTexture;
uniform sampler2D ourTexture;

in vec2 TexCoord;
out vec4 FragColor;

uniform float zNear;
uniform float zFar;
uniform float threshold=0.02;
uniform vec3  viewPos;

float getRealDepth(float z_b){
	float z_n = 2.0*z_b - 1.0;
	float z_e = (z_b*(zNear - zFar)+2*zNear*zFar)/(zNear+zFar);
	float z_final = (z_e-zNear)/(zFar-zNear);
	return z_final;
}
void main(){
	//獲得pixel上下左右的uv位置
	vec2 size = textureSize(depthTexture,0);
	float interval_x = 1/(size.x);
	float interval_y = 1/(size.y);
	vec2 offsets[8];//= {vec2(1,0),vec2(0,-1),vec2(-1.0),vec2(0,1)};
	offsets[0] = vec2(-interval_x,-interval_y);//左上
	offsets[1] = vec2(0,-interval_y);//上
	offsets[2] = vec2(interval_x,-interval_y);//右上
	offsets[3] = vec2(interval_x,0);//右
	offsets[4] = vec2(interval_x,interval_y);//右下
	offsets[5] = vec2(0,interval_y);//下
	offsets[6] = vec2(-interval_x,interval_y);//左下
	offsets[7] = vec2(-interval_x,0);//左
	//索伯算子矩陣x
	float sobel_x[8];
	sobel_x[0] = 1.0;
	sobel_x[1] = 0.0;
	sobel_x[2] = -1.0;
	sobel_x[3] = -2.0;
	sobel_x[4] = -1.0;
	sobel_x[5] = 0.0;
	sobel_x[6] = 1.0;
	sobel_x[7] = 2.0;
	//索伯算子矩陣y
	float sobel_y[8];
	sobel_y[0] = 1.0;
	sobel_y[1] = 2.0;
	sobel_y[2] = 1.0;
	sobel_y[3] = 0.0;
	sobel_y[4] = -1.0;
	sobel_y[5] = -2.0;
	sobel_y[6] = -1.0;
	sobel_y[7] = 0.0;
	/*
	float Gx =0.0;
	float Gy =0.0;
	for(int i=0;i<8;i++){
		float z_b =texture(depthTexture,TexCoord+offsets[i]).x;
		Gx += z_b*sobel_x[i];
		Gy += z_b*sobel_y[i];
	}
	float ans = sqrt(pow(Gx,2)+pow(Gy,2));
	//float depth = getRealDepth(z_b);
	//FragColor = vec4(depth,depth,depth,1);
	if(ans >threshold){
		FragColor = vec4(1.0,1.0,1.0,1.0);
	}
	else{
		FragColor = vec4(0.0,0.0,0.0,1);
	}*/
	float Gx =0.0;
	float Gy =0.0;
	for(int i=0;i<8;i++){
		vec4 normal = texture(ourTexture,TexCoord+offsets[i]);
		vec3 real_normal = (normal.xyz-vec3(0.5,0.5,0.5))*2;
		float mag= dot(real_normal,viewPos)/(length(real_normal)*length(viewPos));
		Gx += mag*sobel_x[i];
		Gy += mag*sobel_y[i];
	}
	float ans = sqrt(pow(Gx,2)+pow(Gy,2));
	FragColor = vec4(ans,ans,ans,1);
}