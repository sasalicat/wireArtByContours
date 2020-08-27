#include "stdafx.h"
#include "tensorMap.h"

mat22 kernel_node(float U, float theta){
	float x = U*cosf(theta);
	float y = U*sinf(theta);
	return mat22(x*x-y*y,-2*x*y,-2*x*y,-(x*x-y*y));
}
tensorMap::tensorMap(MyMesh mesh, int poleIndex,mat22 (*kernel)(float U,float theta)) :polarMap(mesh,poleIndex)
{
	tensor = new mat22[mesh.n_vertices()];
	this->kernel = kernel;
	mat22 m1 = mat22(1, 2, 3, 2);
	m1.eigens();
	mat22 m2 = mat22(9, 1, 2, 3);
	m2.eigens();
	mat22 m3 = mat22(15, 7, 47, 21);
	m3.eigens();
}


tensorMap::~tensorMap()
{
}

void tensorMap::createTensor(){

	for (int i = 0; i < mesh.n_vertices(); i++){
			tensor[i] = kernel(U[i], theta[i]);
			float weight = exp2f(-U[i]);
			tensor[i] = tensor[i] * weight;
	}
}