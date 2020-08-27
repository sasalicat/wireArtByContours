#include "stdafx.h"
#include "geodesicPathFinder.h"
#include <geodesic/geodesic_algorithm_exact.h>

geodesicPathFinder::geodesicPathFinder(){

}
geodesicPathFinder::geodesicPathFinder(vector<float> points, vector<int> faces)
{
	mesh.initialize_mesh_data(points, faces);
}
geodesicPathFinder::geodesicPathFinder(vector<float> points, vector<int> faces,int num)
{
	mesh.initialize_mesh_data(points, faces);
	subFinders.resize(num);
	boundaries.resize(num);
}

geodesicPathFinder::~geodesicPathFinder()
{
}
void geodesicPathFinder::initMesh(vector<float> points, vector<int> faces, int num){
	mesh.initialize_mesh_data(points, faces);
	subFinders.resize(num);
	boundaries.resize(num);
}
void geodesicPathFinder::logMesh(){

}
void geodesicPathFinder::initSubFinder(int index, vector<int> boundary_vIdxs){
	subFinders[index] = new geodesic::GeodesicAlgorithmExact(&mesh);
	boundaries[index].resize(boundary_vIdxs.size());
	for (size_t i = 0; i < boundary_vIdxs.size(); i++)
	{
		boundaries[index][i] = geodesic::SurfacePoint(&mesh.vertices()[boundary_vIdxs[i]]);
	}
}
void geodesicPathFinder::testMem(int start_vIdx, int end_vIdx){
	geodesic::GeodesicAlgorithmExact algorithm(&mesh);
	vector<geodesic::SurfacePoint> path;

	geodesic::SurfacePoint source(&mesh.vertices()[start_vIdx]);
	vector<geodesic::SurfacePoint> all_source(1, source);
	geodesic::SurfacePoint traget(&mesh.vertices()[end_vIdx]);
	vector<geodesic::SurfacePoint> stop_points(1, traget);
	algorithm.propagate(all_source, geodesic::GEODESIC_INF, &stop_points);
	algorithm.trace_back(traget, path);
	algorithm.clear();
}
vector<float*> geodesicPathFinder::getPathBtw(int start_vIdx, int end_vIdx){
	geodesic::GeodesicAlgorithmExact algorithm(&mesh);
	vector<geodesic::SurfacePoint> path;
	vector<float*> result;

	geodesic::SurfacePoint source(&mesh.vertices()[start_vIdx]);
	vector<geodesic::SurfacePoint> all_source(1, source);
	geodesic::SurfacePoint traget(&mesh.vertices()[end_vIdx]);
	vector<geodesic::SurfacePoint> stop_points(1, traget);
	algorithm.propagate(all_source, geodesic::GEODESIC_INF, &stop_points);
	algorithm.trace_back(traget, path);
	for (int i = 0; i < path.size(); i++){
		//cout << "path[" << i << "]:(" << path[i].x()<<","<<path[i].y()<<","<<path[i].z()<<") ";
		float* point=new float[3]{(float)path[i].x(),(float)path[i].y(),(float)path[i].z()};
		result.push_back(point);
	}
	cout << endl;
	return result;
}
vector<float*> geodesicPathFinder::getPathBtw_face(int start_fIdx, int end_fIdx){
	float time_start = clock();
	geodesic::GeodesicAlgorithmExact algorithm(&mesh);
	vector<geodesic::SurfacePoint> path;
	vector<float*> result;

	geodesic::SurfacePoint source(&mesh.faces()[start_fIdx]);
	vector<geodesic::SurfacePoint> all_source(1, source);
	geodesic::SurfacePoint traget(&mesh.faces()[end_fIdx]);
	vector<geodesic::SurfacePoint> stop_points(1, traget);
	algorithm.propagate(all_source, geodesic::GEODESIC_INF, &stop_points);
	algorithm.trace_back(traget, path);
	for (int i = 0; i < path.size(); i++){
		float* point = new float[3]{(float)path[i].x(), (float)path[i].y(), (float)path[i].z()};
		result.push_back(point);
	}
	cout << endl;
	float cost_time = (clock() - time_start) / CLOCKS_PER_SEC;
	cout << "Œ¤ÕÒœyµØ¾€ºÄ•r:" << cost_time << endl;
	return result;
}
vector<float*> geodesicPathFinder::getPathInSubArea(int areaIndex, float* point1, int face1_idx, float* point2, int face2_idx){
	float time_start = clock();
	vector<geodesic::SurfacePoint> path;
	vector<float*> result;
	geodesic::SurfacePoint source(&mesh.faces()[face1_idx]);
	source.set(point1[0],point1[1],point1[2]);
	vector<geodesic::SurfacePoint> all_source(1, source);

	geodesic::SurfacePoint traget(&mesh.faces()[face2_idx]);
	traget.set(point2[0], point2[1], point2[2]);
	subFinders[areaIndex]->propagate(all_source, geodesic::GEODESIC_INF, &boundaries[areaIndex]);
	subFinders[areaIndex]->trace_back(traget, path);
	for (int i = 0; i < path.size(); i++){
		float* point = new float[3]{(float)path[i].x(), (float)path[i].y(), (float)path[i].z()};
		result.push_back(point);
	}
	subFinders[areaIndex]->clear();
	float cost_time = (clock() - time_start) / CLOCKS_PER_SEC;
	return result;
}
