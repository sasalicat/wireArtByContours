#include "stdafx.h"
#include "meshCreater.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

meshCreater::meshCreater()
{
}


meshCreater::~meshCreater()
{
}

bool contain(vector<pathNode*> list, pathNode* elem){
	for each (pathNode *item in list)
	{
		if (item == elem)
		{
			return true;
		}
	}
	return false;
}
void meshCreater::createFaceBtw(OpenMesh::Vec3f p1, OpenMesh::Vec3f p2, OpenMesh::Vec3f normal){
	const float faceWidth = 0.005f;
	OpenMesh::Vec3f VecX = (p2 - p1).normalize();
	OpenMesh::Vec3f tempZ = normal.normalize();
	if (VecX[0] == tempZ[0] && VecX[1] == tempZ[1] && VecX[2] == tempZ[2]){
		cout << "p2-p1與normal的方向相同,無法生成表面!";
		return;
	}
	OpenMesh::Vec3f VecY = OpenMesh::cross(VecX,tempZ).normalize();
	OpenMesh::Vec3f vertex_pos[4];
	vertex_pos[0] = p1 + VecY*(faceWidth / 2);
	vertex_pos[1] = p1 - VecY*(faceWidth / 2);
	vertex_pos[2] = p2 + VecY*(faceWidth / 2);
	vertex_pos[3] = p2 - VecY*(faceWidth / 2);
	MyMesh::VertexHandle vhandles[4];
	vhandles[0] = mesh.add_vertex(vertex_pos[0]);
	vhandles[1] = mesh.add_vertex(vertex_pos[1]);
	vhandles[2] = mesh.add_vertex(vertex_pos[2]);
	vhandles[3] = mesh.add_vertex(vertex_pos[3]);

	vector<MyMesh::VertexHandle> face1(3);
	face1[0] = vhandles[0];
	face1[1] = vhandles[1];
	face1[2] = vhandles[2];
	mesh.add_face(face1);

	vector<MyMesh::VertexHandle> face2(3);
	face2[0] = vhandles[1];
	face2[1] = vhandles[3];
	face2[2] = vhandles[2];
	mesh.add_face(face2);
}

void meshCreater::createCylinderBtw(OpenMesh::Vec3f p1, OpenMesh::Vec3f p2, OpenMesh::Vec3f zAxis){
	OpenMesh::Vec3f AxisX = (p2 - p1).normalize();
	OpenMesh::Vec3f tempZ = zAxis.normalize();
	OpenMesh::Vec3f AxisY = OpenMesh::cross(AxisX, tempZ).normalize();
	OpenMesh::Vec3f AxisZ = OpenMesh::cross(AxisY,AxisX).normalize();
	
	float interval = 1 /(float)cylinderPointNum;
	vector<MyMesh::VertexHandle> vhandles(cylinderPointNum * 2);
	for (size_t i = 0; i < cylinderPointNum; i++)
	{
		float y1 = cos(interval*i * 2 * M_PI)*cylinderRadius;
		float z1 = sin(interval*i * 2 * M_PI)*cylinderRadius;
		OpenMesh::Vec3f v1= p1+(AxisY*y1 + AxisZ*z1);
		OpenMesh::Vec3f v2 =p2+(AxisY*y1 + AxisZ*z1);
		vhandles[i * 2] = mesh.add_vertex(v1);
		vhandles[i * 2 + 1] = mesh.add_vertex(v2);
	}
	for (size_t i = 0; i < cylinderPointNum; i++)
	{
		OpenMesh::VertexHandle v0 = vhandles[i * 2];
		OpenMesh::VertexHandle v1 = vhandles[i * 2 + 1];
		OpenMesh::VertexHandle v2 = vhandles[(i * 2 + 2) % vhandles.size()];
		OpenMesh::VertexHandle v3 = vhandles[(i * 2 + 3) % vhandles.size()];
		vector<OpenMesh::VertexHandle> face1(3);
		face1[0] = v0;
		face1[1] = v1;
		face1[2] = v2;
		mesh.add_face(face1);

		vector<OpenMesh::VertexHandle> face2(3);
		face2[0] = v1;
		face2[1] = v3;
		face2[2] = v2;
		mesh.add_face(face2);
	}

}
void meshCreater::create(GraphArray nodeGraph, string savePath){
	const OpenMesh::Vec3f globalY(0, 1, 0);
	vector<pathNode*> haveBuild;
	vector<pathNode*> candidates;
	candidates.push_back(nodeGraph.array[0][0]);
	while (candidates.size() > 0){
		pathNode* nowNode = candidates[0];
		haveBuild.push_back(nowNode);
		candidates.erase(candidates.begin());
		for (int i = 0; i < nowNode->connects.size(); i++){
			pathNode *next = nowNode->connects[i];
			if (!contain(haveBuild, next)){
				if (nowNode->path[i].size() == 0){
					//createFaceBtw(nowNode->position, next->position, globalY);
					createCylinderBtw(nowNode->position, next->position, globalY);
				}
				else
				{
					//createFaceBtw(nowNode->position, nowNode->path[i].front(), globalY);
					createCylinderBtw(nowNode->position, nowNode->path[i].front(), globalY);
					for (int j = 0; j < nowNode->path[i].size() - 1; j++){
						//createFaceBtw(nowNode->path[i][j], nowNode->path[i][j + 1], globalY);
						createCylinderBtw(nowNode->path[i][j], nowNode->path[i][j + 1], globalY);
					}
					//createFaceBtw(nowNode->path[i].back(), next->position, globalY);
					createCylinderBtw(nowNode->path[i].back(), next->position, globalY);
				}
				if (!contain(candidates, next)){
					candidates.push_back(next);
				}
			}
		}
	}
	try
	{
		if (!OpenMesh::IO::write_mesh(mesh, savePath))
		{
			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
		}
	}
	catch (std::exception& x)
	{
		std::cerr << x.what() << std::endl;
	}
}