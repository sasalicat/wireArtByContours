#pragma once
#include "stdlib.h";
#include <geodesic/geodesic_algorithm_exact.h>
using namespace std;
//比較老的geodesic尋路版本,有持續累積佔用記憶體的問題,當前版本使用geodesicAsynFinder作為代替

class geodesicPathFinder
{
public:
	geodesicPathFinder();
	geodesicPathFinder(vector<float> points,vector<int> faces);
	geodesicPathFinder(vector<float> points, vector<int> faces,int subFinderNum);
	void initMesh(vector<float> points, vector<int> faces, int subFinderNum);
	void initSubFinder(int index, vector<int> boundary_vIdxs);
	~geodesicPathFinder();
	void testMem(int start_vIdx, int end_vIdx);
	vector<float*> getPathBtw(int start_vIdx, int end_vIdx);
	vector<float*> getPathBtw_face(int start_fIdx, int end_fIdx);
	vector<float*> getPathInSubArea(int areaIndex,float* point1,int face1_idx,float* point2,int face2_idx);
	void logMesh();
	void testMultipleInput();
protected:
	geodesic::Mesh mesh;
	vector<geodesic::GeodesicAlgorithmExact*> subFinders;
	vector<vector<geodesic::SurfacePoint>> boundaries;
};

