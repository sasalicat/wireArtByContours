#pragma once
#include "stdlib.h"
#include "polarMap.h"
#include "iostream"
#include <vector>
#include "OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh"
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>

using namespace std;
//用於創建3d mesh

struct pathNode
{
	int id=-1;
	float x_2d=-1;
	float y_2d=-1;//x_2d,y_2d實際上是一組float的id,用於確定兩個pathNode是不是同一個,並沒有實際的坐標意義
	OpenMesh::Vec3f position;
	OpenMesh::Vec3f normal;//沒有實際用途
	vector<pathNode*> connects;
	vector<int> pixelIds;
	map<pathNode*, pathNode*> directedConnects;//key是connect中的值的id, value都是connects中的的值,用途是在線重組的時候確定順序,從key出發經過次node 后下一個node會是value.但如果沒有從key出發經過node的線條的話,directedConnects就不會有value
	//理論上directedConnects的elem數量會是1/2 connect
	vector<vector<OpenMesh::Vec3f>> path;
	pathNode(int nid,float x, float y, OpenMesh::Vec3f pos,int pid){
		id = nid;
		x_2d = x;
		y_2d = y;
		position = pos;
		pixelIds.push_back(pid);
	}
	bool hasConnect(pathNode* other){
		for each (pathNode* node in connects)
		{
			if (node == other){
				return true;
			}
		}
		return false;
	}
};
struct path
{
	float start_2d[2];
	float end_2d[2];
	int pixelId[2];
	OpenMesh::Vec3f normals[2];
	vector<OpenMesh::Vec3f> point3ds;

};
struct GraphArray
{
	int nodeNum = 0;
	vector<vector<pathNode*>> array;
	const float mergeThreshold_x = 0.001f;
	GraphArray(){

	}
	GraphArray(int contourNum){
		array.resize(contourNum);
	}
	pathNode* getNodeBy(int id){
		for (int i = 0; i < array.size(); i++){
			for each (pathNode* node in array[i])
			{
				if (node->id == id){
					return node;
				}
			}
		}
		return NULL;
	}
	void addAPath(path newPath){
		
		int arrayIndex1 = floorf(newPath.start_2d[1]);
		int pixelId1 = newPath.pixelId[0];
		pathNode* node1 = NULL;
		for each (pathNode* node in array[arrayIndex1])
		{
			if (abs(node->x_2d - loopClampTo01(newPath.start_2d[0]))<mergeThreshold_x && node->y_2d == newPath.start_2d[1]){//找到同一個node
				node1 = node;
				if(!contain(node1->pixelIds,pixelId1))
				{
					node1->pixelIds.push_back(pixelId1);
				}
				//cout << "连接原node1:(" << node->x_2d << "," << node->y_2d << ")";
				break;
			}
		}
		if (node1 == NULL){
			
			pathNode* newNode = new pathNode(nodeNum++,loopClampTo01(newPath.start_2d[0]),newPath.start_2d[1],newPath.point3ds.front(),pixelId1);
			//cout << "创建新node1:(" <<newNode->x_2d<<","<<newNode->y_2d<<")";
			array[arrayIndex1].push_back(newNode);
			node1 = newNode;
		}

		int arrayIndex2 = floorf(newPath.end_2d[1]);
		int pixelId2 = newPath.pixelId[1];
		//cout << "newPath.end_2d[1]:" << newPath.end_2d[1] << "arrayIndex2:" << arrayIndex2 << ";";
		pathNode* node2 = NULL;
		for each (pathNode* node in array[arrayIndex2])
		{
			if (abs(node->x_2d- loopClampTo01(newPath.end_2d[0]))<mergeThreshold_x && node->y_2d == newPath.end_2d[1]){//找到同一個node
				node2 = node;
				if (!contain(node2->pixelIds, pixelId2)){
					node2->pixelIds.push_back(pixelId2);
				}
				//cout << "连接原node2:(" << node->x_2d << "," << node->y_2d << ")";
				break;
			}
		}
		if (node2 == NULL){
			
			pathNode* newNode = new pathNode(nodeNum++, loopClampTo01(newPath.end_2d[0]), newPath.end_2d[1], newPath.point3ds.back(),pixelId2);
			//cout << "创建新node2:(" << newNode->x_2d << "," << newNode->y_2d << ")";
			array[arrayIndex2].push_back(newNode);
			node2 = newNode;
		}
		//cout << "arrayIndex1:" << arrayIndex1 << " arrayIndex2:" << arrayIndex2;
		if (!node1->hasConnect(node2))
		{
			//cout << "创建node1 connect to node2";
			node1->connects.push_back(node2);
			if (newPath.point3ds.size() > 2){
				vector<OpenMesh::Vec3f> newone;
				for (size_t i = 1; i < newPath.point3ds.size()-1; i++)
				{
					newone.push_back(newPath.point3ds[i]);
				}
				node1->path.push_back(newone);
			}
			else{
				node1->path.push_back(vector<OpenMesh::Vec3f>());
			}
		}
		if (!node2->hasConnect(node1)){
			//cout << "创建node2 connect to node1";
			node2->connects.push_back(node1);
			if (newPath.point3ds.size()>2){
				vector<OpenMesh::Vec3f> newone;
				for (size_t i = newPath.point3ds.size()-2; i >= 1; i--)
				{
					newone.push_back(newPath.point3ds[i]);
				}
				node2->path.push_back(newone);
			}
			else{
				node2->path.push_back(vector<OpenMesh::Vec3f>());
			}
		}
	}
};
class meshCreater
{
public:
	const float cylinderRadius = 0.01f;//0.01f;
	const float cylinderPointNum = 6;
	meshCreater();
	~meshCreater();
	void create(GraphArray nodeGraph,string savePath);
	MyMesh mesh;
	void createFaceBtw(OpenMesh::Vec3f p1,OpenMesh::Vec3f p2,OpenMesh::Vec3f normal);
	void createCylinderBtw(OpenMesh::Vec3f p1, OpenMesh::Vec3f p2, OpenMesh::Vec3f zAxis);
};
bool contain(vector<pathNode*> list, pathNode* elem);
