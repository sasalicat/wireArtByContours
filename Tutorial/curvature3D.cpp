#include "stdafx.h"
#include "curvature3D.h"
#include <stdio.h>
#include <iostream>
glm::vec2 solve2x3(glm::mat2x3 mat,bool &fail){
	if ((mat[0][0] == 0 && mat[1][0] == 0) || (mat[1][1] == 0 && mat[0][1] == 0)){
		fail = true;
		return glm::vec2(0, 0);
	}
	else
	{
		fail = false;
		if (mat[0][0] == 0){
			glm::vec3 temp = mat[0];
			mat[0] = mat[1];
			mat[1] = temp;
		}
		mat[0] /= mat[0][0];//歸一化第1行
		mat[1] -= mat[0] * mat[1][0];//將第二行第一個數字歸零
		mat[1] /= mat[1][1];//歸一化第2行
		mat[0] -= mat[1] * mat[0][1];//將第一行第二個數字歸零
		return glm::vec2(mat[0][2], mat[1][2]);
	}
}
curvature3D::curvature3D(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	float l21 = glm::length(p1 - p2);
	float l23 = glm::length(p3 - p2);
	float cosValue = glm::dot((p1 - p2), (p3 - p2))/(l21*l23);
	float sinValue = sqrtf(1 - cosValue*cosValue);
	point2Ds[0] = glm::vec2(glm::length(p1 - p2),0);
	point2Ds[1] = glm::vec2(0, 0);
	point2Ds[2] = glm::vec2(l23*cosValue, l23*sinValue);
}
curvature3D::curvature3D(glm::vec2 p1, glm::vec2 p2, glm::vec2 p3){
	point2Ds[0] = p1;
	point2Ds[1] = p2;
	point2Ds[2] = p3;
}
curvature3D::curvature3D(){
	
}

curvature3D::~curvature3D()
{
}
float curvature3D::getResult(){

	bool fail = true;
	float ta = glm::length(point2Ds[0] - point2Ds[1]);
	float tb = glm::length(point2Ds[2] - point2Ds[1]);
	glm::mat2x3 matrix_x;
	matrix_x[0] = glm::vec3(-ta, ta*ta, point2Ds[0][0] - point2Ds[1][0]);
	matrix_x[1] = glm::vec3(tb, tb*tb, point2Ds[2][0] - point2Ds[1][0]);

	glm::vec2 ans= solve2x3(matrix_x,fail);
	float a1 = point2Ds[1][0];
	float a2 = ans[0];
	float a3 = ans[1];
	glm::mat2x3 matrix_y;
	matrix_y[0] = glm::vec3(-ta, ta*ta, point2Ds[0][1] - point2Ds[1][1]);
	matrix_y[1] = glm::vec3(tb, tb*tb, point2Ds[2][1] - point2Ds[1][1]);

	ans = solve2x3(matrix_y, fail);
	float b1 = point2Ds[1][1];
	float b2 = ans[0];
	float b3 = ans[1];
	return 2 * (a3*b2 - a2*b3) / powf(a2*a2 + b2*b2, 3 / 2);
}