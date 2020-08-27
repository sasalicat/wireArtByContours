#version 330

uniform sampler2D depthTexture;
uniform sampler2D ourTexture;

in vec2 TexCoord;
out vec4 FragColor;

uniform float zNear;
uniform float zFar;
uniform float threshold=0.35;
uniform vec3  viewPos;
uniform mat4 gMV;
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
	
	//用depthmap找到強輪廓---------------------------------------------------
	
	int edge = 0;
	if(texture(depthTexture,TexCoord).x<1){
		if(texture(depthTexture,TexCoord+offsets[1]).x==1){
			edge=1;
		}
		else if(texture(depthTexture,TexCoord+offsets[3]).x==1){
			edge =1;
		}
		else if(texture(depthTexture,TexCoord+offsets[5]).x==1){
			edge =1;
		}
		else if(texture(depthTexture,TexCoord+offsets[7]).x==1){
			edge =1;
		}
	}
	if(edge ==1){
		FragColor = vec4(1,1,1,1);
	}
	//用dotmap做索伯邊緣偵測-------------------------------------------------
	
	else{
		
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
		float Gx =0.0;
		float Gy =0.0;
		for(int i=0;i<8;i++){
			vec4 pixel = texture(ourTexture,TexCoord+offsets[i]);
			float mag = pixel[0];
			Gx += mag*sobel_x[i];
			Gy += mag*sobel_y[i];
		}
		float ans_ori = sqrt(pow(Gx,2)+pow(Gy,2));
		float ans = ans_ori/4;
		
		if(ans>threshold){
			ans = 1;
		}
		else{
			ans = 0;
		}
		FragColor = vec4(ans,ans,ans,1);
		
	}
	
	//float depth = texture(depthTexture,TexCoord).x;
	//FragColor = vec4(depth,depth,depth,1);
}