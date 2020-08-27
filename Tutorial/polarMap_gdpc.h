#pragma once
#include "polarMap.h"
#include "GDPC\Mesh.h";
#include "GDPC\Generator.h"
typedef DGPC::Vector3<double> Point;
typedef DGPC::MeshOM<Point> omMesh;
typedef DGPC::Generator<omMesh> DGPCgenerator;
//用來測試polar map原方法的代碼,並不能拿來用,基本上的目的就就是drawMeshToImage看看效果
class polarMap_gdpc:public polarMap
{
public:
	omMesh mesh_om;
	DGPCgenerator *generator=NULL;
	polarMap_gdpc(string path,int polarIdx);
	int pointNum();
	void createMap();
	void drawMeshToImage();
	int getFasthestIdx();
	polarPoint getPointFrom(int vertexIndex);
	void findExtremumVertexs();
	~polarMap_gdpc();
private:
	int endIdx;
};

