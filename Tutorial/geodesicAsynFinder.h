#pragma once
#include "stdlib.h"
#include <fstream>
#include <stdio.h>
#include <vector>
#include <sstream>
#include "OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh"
using namespace std;
//用於異步尋路,
//為何要使用異步尋路的原因是geodesic_algorithm_exact這個代碼會導致持續佔用內存,所以geodesicAsynFinder的做法是先把請求存到mission.txt里,在用system call調用geodesicGroupFinder.exe找到路徑后存到pathResults.txt,然後再讀取,
//這樣做雖然在geodesicGroupFinder.exe運行的時候仍然會持續佔用內存,但是當geodesicGroupFinder.exe結束時佔用的內存就會被施放
struct gdc_pair
{
	int fIndex1;
	float point1[3];
	int fIndex2;
	float point2[3];
	gdc_pair(){

	}
	gdc_pair(int fIdx1, float x1, float y1, float z1, int fIdx2, float x2, float y2, float z2){
		fIndex1 = fIdx1;
		point1[0] = x1;
		point1[1] = y1;
		point1[2] = z1;

		fIndex2 = fIdx2;
		point2[0] = x2;
		point2[1] = y2;
		point2[2] = z2;
	}
	string to_string(){
		ostringstream oss;
		oss << fIndex1 << " " << point1[0] << " " <<point1[1] << " " << point1[2] << " " <<fIndex2 << " " << point2[0] << " " << point2[1] << " " << point2[2] << ";";
		return string(oss.str());
	}
};
struct  marker
{
	int cindex = 0;
	int lindex = 0;
};
const int READ_PATH_BUFFER_SIZE = 2048;
class geodesicAsynFinder
{
protected:

	const string PATH_FILE_RESULT = "pathResults.txt";
	const string PATH_FILE_BOUNDARY = "boundaryIndexs.txt";
	const string PATH_FILE_MISSION = "mission.txt";
	vector<vector<int>> boundarys;
	vector<vector<gdc_pair>> missionList;
	marker nowPointer;
	string fileName;
public:

	geodesicAsynFinder();
	~geodesicAsynFinder();
	void init(int contourNum,string fname);
	void initBoundaryRecord(int index,vector<int> boundary);
	void writeBoundaryData();
	vector<vector<OpenMesh::Vec3f>> readPathResult();
	void requestPath(int cindex, int faceIdx_s, float* pos_s, int faceIdx_e, float* pos_e);
	bool writeToMissionFile(int missionNum);
	void doIt();
	void printMissionNum();
};

