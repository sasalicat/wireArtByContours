#pragma once
#include "polarMap.h"
//用來計算tensorMap的代碼,並沒有寫完,實際上也沒有使用
struct mat22{
	OpenMesh::Vec2f num[2];
	mat22(float e00, float e10, float e01, float e11){
		num[0] = OpenMesh::Vec2f(e00, e01);
		num[1] = OpenMesh::Vec2f(e10, e11);
	}
	mat22(){
		num[0] = OpenMesh::Vec2f(1, 0);
		num[1] = OpenMesh::Vec2f(0, 1);
	}
	mat22 operator*(float c){
		return mat22(num[0][0]*c,num[0][1]*c,num[1][0]*c,num[1][1]*c);
	}
	OpenMesh::Vec2f* eigens(){
		float a = 0;
		float b = 0;
		float c = 0;
		cout << "matrix:" << endl << num[0][0] << "," << num[0][1] << endl << num[1][0] << "," << num[1][1] << endl;
		c += num[1][1]*num[0][0];
		cout << "c += " << num[1][1] * num[0][0];
		b -= num[0][0];
		b -= num[1][1];
		a += 1;
		c -= num[0][1] * num[1][0];
		cout << "c -= " << num[0][1] * num[1][0];
		cout << "expr:" << a << "x^2 +" << b << "x +" << c << endl;
		float landa1 = (-b + sqrtf(powf(b,2)-4*a*c))/(2 * a);
		float landa2 = (-b - sqrtf(powf(b, 2) - 4 * a*c)) / (2 * a);
		//for landa1
		OpenMesh::Vec2f row1 = num[0]- OpenMesh::Vec2f(landa1,0);
		OpenMesh::Vec2f row2 = num[1] - OpenMesh::Vec2f(0, landa1);
		
		cout << "landa1:" << landa1 << " row1:" << row1[0] << "," << row1[1] << " row2:" << row2[0] << "," << row2[1] << endl;
		//
		return NULL;
	}
};
class tensorMap :public polarMap
{
public:
	tensorMap(MyMesh mesh, int poleIndex, mat22(*kernel)(float U, float theta));
	~tensorMap();
	void createTensor();
protected:
	mat22(*kernel)(float U, float theta);
	mat22 *tensor;
};

mat22 kernel_node(float U, float theta);