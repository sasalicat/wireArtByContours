#pragma once

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <ogldev_math_3d.h>
#include "stdlib.h"
#include "OpenMesh\Core\IO\MeshIO.hh"
#include "OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh"
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
#include "iostream"
#include "drawMap.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "geodesicPathFinder.h"
#include "geodesicAsynFinder.h"
#include "curvature3D.h";
#include "IS_FeatureMatrix.h";
//#include "Tutorial.cpp"
using namespace std;

//主要使用的類別,用於計算polarMap和創建等高線,計算角度對應關係.以及尋找等高線之間的路徑,包括投影路徑和geodesic
typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;
float loopFindNearestCycle(float traget, float testAngle);
float loopFindNearestBiggerValue(float traget, float testAngle);
float loopFindNearestSmallerValue(float traget, float testAngle);
float loopTransToLayer(float value, int layer);
float loopClampTo01(float oriAngle);
struct Vec3f_angle{
	OpenMesh::Vec3f pos;
	float angle;
	Vec3f_angle(OpenMesh::Vec3f positon,float theta){
		pos = positon;
		angle = theta;
	}
	Vec3f_angle(){
		angle = -1;
	}
};
struct float_pair{
	float values[2];
	float operator[] (int index){
		return values[index];
	}
	float_pair(){

	}
	float_pair(float v1, float v2){
		values[0] = v1;
		values[1] = v2;
	}
};
struct divergenceGridInf
{
	int gridIndex=-1;
	float angleCorresponds[2];
};
struct fp_mapping
{
	float_pair group[2];
	vector<divergenceGridInf> gridMarkers;
	//vector<float_pair> intersectionMapping;
	float_pair operator[] (int index){
		return group[index];
	}
	fp_mapping(){

	}
	fp_mapping(float_pair p1, float_pair p2){
		group[0] = p1;
		group[1] = p2;
	}
};

struct PointCurve{
	vector<Vec3f_angle> points;
	//vector<PointCurve*> nextCurve;
	void Add(Vec3f_angle va){
		//cout << "PointCurve Add" << va.angle<<" pos:"<<va.pos[0]<<","<<va.pos[1]<<","<<va.pos[2]<<" ";
		bool find = false;
		for (int i = 0; i < points.size(); i++){
			if (va.angle < points[i].angle){
				find = true;
				points.insert(points.begin()+i,va);
				break;
			}
		}
		if (!find){
			points.push_back(va);
		}
	}
	OpenMesh::Vec3f getPointAt(float anglePercentage){
		for (int p = 0; p < points.size()-1; p++){
			float lastAngle = points[p].angle;
			float nextAngle = points[p+1].angle;
			if (points[p].angle <= anglePercentage&&points[p + 1].angle>anglePercentage){
				OpenMesh::Vec3f toNext = points[p + 1].pos - points[p].pos;
				float percentage = (anglePercentage - lastAngle) / (nextAngle - lastAngle);
				OpenMesh::Vec3f now_toNext = toNext*percentage;
				return points[p].pos + now_toNext;
			}
		}
		int lastIdx = points.size() - 1;
		float lastAngle = points[lastIdx].angle;
		float nextAngle = 1;
		OpenMesh::Vec3f toNext = points[0].pos - points[lastIdx].pos;//从上个点到下个点的向量
		float percentage = (anglePercentage - lastAngle) / (nextAngle - lastAngle);
		OpenMesh::Vec3f now_toNext = toNext*percentage;
		return points[lastIdx].pos + now_toNext;
	}
};
struct pointFromTwoSource{
	int sourceIdx_Small;
	int sourceIdx_Large;
	OpenMesh::Vec3f pointPos;
	float angleByCal = 0;

	pointFromTwoSource(int v1, int v2, OpenMesh::Vec3f pos){
		sourceIdx_Small = v1;
		sourceIdx_Large = v2;
		pointPos = pos;
		//cout << "pfts pos:(" << pointPos[0] << "," << pointPos[1] << "," << pointPos[2] << ");";
	}
	pointFromTwoSource(){
		sourceIdx_Large = -1;
		sourceIdx_Large = -1;
	}
	static OpenMesh::Vec3f Interpolation(pointFromTwoSource p_small, pointFromTwoSource p_large, float angle){
		OpenMesh::Vec3f result(0,0,0);
		//cout << "p_small angleByCal:" << p_small.angleByCal << "p_large angleByCal:" << p_large.angleByCal << endl;
		if (p_small.angleByCal>p_large.angleByCal){//最大值點和最小值點之間
			float length = (1 - p_small.angleByCal) + p_large.angleByCal;
			float subLength = 0;
			if (angle > p_small.angleByCal&&angle > p_large.angleByCal){//此時angle應該在p_small.angleByCal到1之間
				subLength= angle - p_small.angleByCal;
			}
			else if (angle < p_small.angleByCal&&angle < p_large.angleByCal){
				subLength = (1 - p_small.angleByCal) + angle;
			}
			return p_small.pointPos + (p_large.pointPos - p_small.pointPos)*(subLength / length);
		}
		else if (angle > p_large.angleByCal||angle<p_small.angleByCal){
			cout << "pointFromTwoSource 錯誤的Interpolartion! traget angle:"<<angle<<" p_small:"<<p_small.angleByCal<<"p_large:"<<p_large.angleByCal << endl;
			return  result;
		}
		float percent = (angle - p_small.angleByCal) / (p_large.angleByCal - p_small.angleByCal);
		return p_small.pointPos + (p_large.pointPos - p_small.pointPos)*percent;
	}
	static pointFromTwoSource Interpolation_newpfts(pointFromTwoSource p_small, pointFromTwoSource p_large, float angle){
		pointFromTwoSource result;
		//cout << "p_small angleByCal:" << p_small.angleByCal << "p_large angleByCal:" << p_large.angleByCal << endl;
		if (p_small.angleByCal>p_large.angleByCal){//最大值點和最小值點之間
			float length = (1 - p_small.angleByCal) + p_large.angleByCal;
			float subLength = 0;
			if (angle > p_small.angleByCal&&angle > p_large.angleByCal){//此時angle應該在p_small.angleByCal到1之間
				subLength = angle - p_small.angleByCal;
			}
			else if (angle < p_small.angleByCal&&angle < p_large.angleByCal){
				subLength = (1 - p_small.angleByCal) + angle;
			}
			result.angleByCal = angle;
			result.pointPos = p_small.pointPos + (p_large.pointPos - p_small.pointPos)*(subLength / length);
			return result;
		}
		else if (angle > p_large.angleByCal || angle<p_small.angleByCal){
			cout << "pointFromTwoSource 錯誤的Interpolartion!" << endl;
			return  result;
		}
		float percent = (angle - p_small.angleByCal) / (p_large.angleByCal - p_small.angleByCal);
		result.angleByCal = angle;
		result.pointPos = p_small.pointPos + (p_large.pointPos - p_small.pointPos)*percent;
		return result;
	}
};
struct pointAtContour
{
	OpenMesh::Vec3f pos;
	int contourIdx=-1;
	float atAngle = -1;
	pointAtContour(int idx, float angle, OpenMesh::Vec3f pt){
		contourIdx = idx;
		atAngle = angle;
		pos = pt;
	}
	pointAtContour(){

	}
};
float calEdgeTime(OpenMesh::Vec3f startPoint, OpenMesh::Vec3f edgePoint1, OpenMesh::Vec3f edgePoint2);
OpenMesh::Vec3f calCloestPoint(OpenMesh::Vec3f startPoint, OpenMesh::Vec3f edgePoint1, OpenMesh::Vec3f edgePoint2);
struct int_pair{
	int values[2];
	int_pair(){};
	int_pair(int v1, int v2){
		values[0] = v1;
		values[1] = v2;
	}
	int_pair(int begin){
		values[0] = begin;
	}
	int operator[](int index){
		return values[index];
	}
};
struct contour
{

	int id = -1;
	float distance = 0;
	float length = 0;
	vector<contour*> lastContours;
	vector<contour*> nextContours;
	vector<pointFromTwoSource> pfts;
	vector<int> pointsBelong;//數量和pfts相同,標記每一個pointFromTwoSource對應于哪一個nextContours
	vector<vector<fp_mapping>> divergenceSegm;
	vector<OpenMesh::Vec3f> closetNextPoint;
	vector<float> curvatures;//記錄每一點的曲率
	contour(){
	}
	contour(vector<pointFromTwoSource> traget){
		pfts = traget;
		for (int i = 1; i < traget.size(); i++){
			length += (traget[i].pointPos - traget[i - 1].pointPos).length();
		}
		length += (traget[traget.size() - 1].pointPos - traget[0].pointPos).length();
	}
	contour(vector<pointFromTwoSource> traget,float distance2OrigenPoint){
		pfts = traget;
		distance = distance2OrigenPoint;
		for (int i = 1; i < traget.size(); i++){
			length += (traget[i].pointPos - traget[i - 1].pointPos).length();
		}
		length += (traget[traget.size() - 1].pointPos - traget[0].pointPos).length();
	}
	bool contact(int vidx){
		for each (pointFromTwoSource point in pfts)
		{
			if (point.sourceIdx_Small == vidx){
				return true;
			}
		    if (point.sourceIdx_Large == vidx)
			{
				return true;
			}
		}
		return false;
	}
	void debug_contact(int vidx){
		bool found = false;
		for each (pointFromTwoSource point in pfts)
		{
			if (point.sourceIdx_Small == vidx){
				cout << vidx << "是sorceSmall large是:" << point.sourceIdx_Large<<",";
				found = true;
			}
			if (point.sourceIdx_Large == vidx)
			{
				cout << vidx << "是sorceLarge small是:" << point.sourceIdx_Small<<",";
				found = true;
			}
		}
		if (found){
			cout << endl;
		}
	}
	bool contact(int vidx1, int vidx2){
		for each (pointFromTwoSource point in pfts)
		{
			if (point.sourceIdx_Small == vidx1 && point.sourceIdx_Large == vidx2){
				return true;
			}
			if (point.sourceIdx_Large == vidx1 && point.sourceIdx_Small == vidx2){
				return true;
			}
		}
		return false;
	}
	bool contact_asSmall(int vidx){
		bool found = false;
		for each (pointFromTwoSource point in pfts)
		{
			if (point.sourceIdx_Small == vidx){
				return true;
			}
		}
		return false;
	}
	pointFromTwoSource get_contact_asSmall_pfts(int vidx){
		for each (pointFromTwoSource point in pfts)
		{
			if (point.sourceIdx_Small == vidx){
				return point;
			}
		}
		return pointFromTwoSource(-1,-1,OpenMesh::Vec3f());
	}
	bool contact_asLarge(int vidx){
		for each (pointFromTwoSource point in pfts)
		{
			if (point.sourceIdx_Large == vidx){
				return true;
			}
		}
		return false;
	}
	int indexOfTwoVidx(int vidx1, int vidx2){
		for (int i = 0; i < pfts.size(); i++){
			if ((pfts[i].sourceIdx_Small == vidx1&&pfts[i].sourceIdx_Large == vidx2) || (pfts[i].sourceIdx_Small == vidx2&&pfts[i].sourceIdx_Large == vidx1)){
				return i;
			}
		}
		return -1;
	}
	PointCurve toCurve(){
		PointCurve curve;
		for each (pointFromTwoSource point in pfts)
		{
			curve.Add(Vec3f_angle(point.pointPos, 0));
		}
		return curve;
	}
	vector<OpenMesh::Vec3f> toVec3fList(){
		vector<OpenMesh::Vec3f> list(pfts.size());
		for (int i = 0; i < pfts.size(); i++){
			list[i] = pfts[i].pointPos;
		}
		return list;
	}
	vector<pointFromTwoSource> searchRange(float start_angle, float end_angle){
		vector<pointFromTwoSource> pointList;
		pointFromTwoSource p11;
		pointFromTwoSource p12;
		int nowIndex = -1;
		bool crossZeroPoint = end_angle < start_angle;//因為end_angle<start_angle說明起點到終點跨過了0點
		cout << "search Range 開始:";
		if (start_angle <= pfts[0].angleByCal&&start_angle>=pfts[0].angleByCal){
			p11 = pfts.back();
			p12 = pfts[0];
			nowIndex = 0;
		}
		else
		{
			
			for (int i = 0; i < pfts.size() - 1; i++){
				if (start_angle >= pfts[i].angleByCal&&start_angle <= pfts[i + 1].angleByCal){
					p11 = pfts[i];
					p12 = pfts[i + 1];
					nowIndex = i + 1;
					break;
				}
				else if (start_angle>=pfts[i].angleByCal&&start_angle>=pfts[i+1].angleByCal&&pfts[i+1].angleByCal<pfts[i].angleByCal)
				{
					p11 = pfts[i];
					p12 = pfts[i + 1];
					nowIndex = i + 1;
					break;
				}
			}
		}
		cout << "p11 angle:" << p11.angleByCal << "p12 angle:" << p12.angleByCal << endl;
		pointList.push_back(pointFromTwoSource::Interpolation_newpfts(p11, p12, start_angle));
		cout << "起始點:" << pointList[0].angleByCal;
		while (true){
			pointList.push_back(pfts[nowIndex]);
			if (nowIndex == 0){
				if (pfts.back().angleByCal <= end_angle&&end_angle <= pfts[0].angleByCal){
					pointList.push_back(pointFromTwoSource::Interpolation_newpfts(pfts.back(), pfts[0], end_angle));
				}
			}
			else if (pfts[nowIndex - 1].angleByCal <= end_angle&&end_angle <= pfts[nowIndex].angleByCal){
				pointList.push_back(pointFromTwoSource::Interpolation_newpfts(pfts[nowIndex - 1], pfts[nowIndex], end_angle));
				break;
			}
			nowIndex++;
			if (nowIndex >= pfts.size()){
				nowIndex = 0;
			}
		}
		return pointList;
	}
	pointFromTwoSource* searchAngle(float angle){
		if (false){
			cout << "contour" << id<<" pfts size:"<<pfts.size()<< "searchAngle:" << angle<<endl;
			for (int i = 0; i < pfts.size(); i++){
				cout << "i=" << i << " angle:" << pfts[i].angleByCal << "; ";
			}
		}
		pointFromTwoSource *result= new pointFromTwoSource[2];
		if ((pfts[0].angleByCal >= angle && pfts.back().angleByCal<= angle)||(pfts[0].angleByCal>=angle&&pfts.back().angleByCal>=angle&&pfts.back().angleByCal>pfts[0].angleByCal)){
			result[0] = pfts.back();
			result[1] = pfts[0];
			//cout << "angle 位於index 0之前";

		}
		else{
			for (int i = 1; i < pfts.size(); i++){
				//cout << "i=" << i << " angleByCal:" << pfts[i].angleByCal << "|"<<endl;
				if (pfts[i].angleByCal >= angle&&pfts[i - 1].angleByCal <= angle){//比前一個大比后一個小,則找到了正確位置
					result[0] = pfts[i - 1];
					result[1] = pfts[i];
					//cout << "比較pfts[i].angleByCal:" << pfts[i].angleByCal << "pfts[i-1].angleByCal:" << pfts[i - 1].angleByCal;
					break;
				}
				else if ((pfts[i].angleByCal >= angle&&pfts[i - 1].angleByCal >= angle) && (pfts[i - 1].angleByCal>pfts[i].angleByCal)){
					result[0] = pfts[i - 1];
					result[1] = pfts[i];
					break;
				}
			}
		}
		//cout << "result[0].angle:" << result[0].angleByCal << " result[1].angle:" << result[1].angleByCal<<endl;
		return result;
	}
	OpenMesh::Vec3f getPosAtAngle(float angle){
		angle = loopClampTo01(angle);
		pointFromTwoSource *pair = searchAngle(angle);
		return pointFromTwoSource::Interpolation(pair[0], pair[1], angle);
	}
	OpenMesh::FaceHandle getFaceAtAngle(MyMesh &mesh, float angle){
		angle = loopClampTo01(angle);
		pointFromTwoSource* pair = searchAngle(angle);
		MyMesh::VertexHandle vh(pair[0].sourceIdx_Small);
		for (MyMesh::VertexFaceCWIter vf_cwit = mesh.vf_cwbegin(vh); vf_cwit != mesh.vf_cwend(vh); vf_cwit++){
			//cout << "搜尋面" << (*vf_cwit).idx()<<":";
			int matchNum = 0;
			for (MyMesh::FaceVertexIter fv_it = mesh.fv_begin(*vf_cwit); fv_it != mesh.fv_end(*vf_cwit); fv_it++){
				int vidx = (*fv_it).idx();
				//cout << "vidx:" << vidx<<" ";
				if (pair[0].sourceIdx_Small == vidx || pair[0].sourceIdx_Large == vidx || pair[1].sourceIdx_Small == vidx || pair[1].sourceIdx_Large == vidx){
					matchNum++;
					//cout << "match" << matchNum << "次. ";
				}
			}
			if(matchNum==3)
			{
				return *vf_cwit;
			}
		}
		return OpenMesh::FaceHandle(-1);
	}

	Vec3f_angle* sample(int sampleNumber){
		int nowIndex = -1;
		for (int i = 0; i < pfts.size(); i++){
			if (abs(pfts[i].angleByCal - 0.f) <0.000001f || abs(1.f - pfts[i].angleByCal) <0.000001f){//這裡寫 pfts[i].angleByCal ==1沒有用,我也不知道為什么 1!=1
				nowIndex = i;
				break;
			}
		}
		if (nowIndex<0)
		{
			cout << "沒有找到角度為0的點,錯誤的contour!" << endl;
			return NULL;
		}
		Vec3f_angle* result = new Vec3f_angle[sampleNumber];
		result[0] = Vec3f_angle(pfts[nowIndex].pointPos, 0);
		//cout << "起始index為:" << nowIndex<<endl;
		int counter = 0;
		float interval = 1.0f / (float)sampleNumber;
		while (counter<sampleNumber-1)
		{
			//cout << "nowIndex:" << nowIndex << " counter:" << counter;
			int nextIndex = (nowIndex + 1) % pfts.size();//使超出範圍的index回到頭部,打個比方pfts.size()=10,nowIndex=9,則nowIndex+1=10,10%10=0,所以nextIndex變為0
			while ((counter+1)*interval<= pfts[nextIndex].angleByCal){//下一個點人在nextIndex的pfts的角度內
				counter++;
				OpenMesh::Vec3f point = pointFromTwoSource::Interpolation(pfts[nowIndex], pfts[nextIndex], counter*interval);
				result[counter] = Vec3f_angle(point,counter*interval);
			}
			nowIndex = nextIndex;

		}
		return result;
	}
	vector<OpenMesh::Vec3f> sample_withoutAngle(int sampleNumber){

		int nowIndex = -1;
		for (int i = 0; i < pfts.size(); i++){
			if (abs(pfts[i].angleByCal - 0.f) <0.000001f || abs(1.f - pfts[i].angleByCal) <0.000001f){//這裡寫 pfts[i].angleByCal ==1沒有用,我也不知道為什么 1!=1
				nowIndex = i;
				break;
			}
		}
		if (nowIndex<0)
		{
			cout << "沒有找到角度為0的點,錯誤的contour!" << endl;
			return vector<OpenMesh::Vec3f>();
		}
		//Vec3f_angle* result = new Vec3f_angle[sampleNumber];
		vector<OpenMesh::Vec3f> result(sampleNumber);
		result[0] =pfts[nowIndex].pointPos;
		//cout << "起始index為:" << nowIndex<<endl;
		int counter = 0;
		float interval = 1.0f / (float)sampleNumber;
		while (counter<sampleNumber - 1)
		{
			//cout << "nowIndex:" << nowIndex << " counter:" << counter;
			int nextIndex = (nowIndex + 1) % pfts.size();//使超出範圍的index回到頭部,打個比方pfts.size()=10,nowIndex=9,則nowIndex+1=10,10%10=0,所以nextIndex變為0
			while ((counter + 1)*interval <= pfts[nextIndex].angleByCal){//下一個點人在nextIndex的pfts的角度內
				counter++;
				OpenMesh::Vec3f point = pointFromTwoSource::Interpolation(pfts[nowIndex], pfts[nextIndex], counter*interval);
				result[counter] = point;
			}
			nowIndex = nextIndex;

		}
		return result;
	}
	float getContourLengthBtw(float angle1, float angle2){
		int minIndex = 0;
		int lastMinIndex = pfts.size() - 1;
		float lengthCount = 0;
		bool foundAngle1 = false;
		for (int i = 1; i < pfts.size(); i++){
			if (pfts[i].angleByCal < pfts[minIndex].angleByCal){
				minIndex = i;
				lastMinIndex = i - 1;
			}
		}
		if (pfts[minIndex].angleByCal < angle1)
		{
			foundAngle1 = true;
			float vecLength = (pfts[lastMinIndex].pointPos - pfts[minIndex].pointPos).length();
			float angleOffset = (1 - pfts[lastMinIndex].angleByCal) + pfts[minIndex].angleByCal;
			float subAngle= pfts[minIndex].angleByCal - angle1;
			lengthCount += vecLength*(subAngle / angleOffset);
		}
		for (int i = 1; i < pfts.size(); i++){
			int nowIndex = (minIndex + i) % pfts.size();
			int lastIndex = (minIndex + i - 1) % pfts.size();
			if (!foundAngle1){//如果還沒有找到angle1位於哪條邊
				if (pfts[nowIndex].angleByCal >= angle1){
					foundAngle1 = true;
					float vecLength = (pfts[lastIndex].pointPos - pfts[nowIndex].pointPos).length();
					float angleOffset = pfts[nowIndex].angleByCal - pfts[lastIndex].angleByCal;
					float subAngle = pfts[nowIndex].angleByCal - angle1;
					lengthCount += vecLength*(subAngle / angleOffset);
				}
			}
			else
			{
				if (pfts[nowIndex].angleByCal >= angle2){//如果找到一個點的角度比angle2大,說明angle2存在于這段邊之中
					float vecLength = (pfts[lastIndex].pointPos - pfts[nowIndex].pointPos).length();
					float angleOffset = pfts[nowIndex].angleByCal - pfts[lastIndex].angleByCal;
					float subAngle = angle2- pfts[lastIndex].angleByCal;
					break;
				}
				else//如果不存在這段邊之中則將邊的長度累加上去
				{
					float vecLength = (pfts[lastIndex].pointPos - pfts[nowIndex].pointPos).length();
					lengthCount += vecLength;
				}
			}
		}
		return lengthCount;
	}
	void calDivergence(){
		//cout << "contour" << id << " calDivergence";
		if (nextContours.size() <= 1){
			//cout << "nextContour.size()<=1!"<<endl;
			return;
		}
		//cout << "is divergence"<<endl;
		vector<OpenMesh::Vec3f> pointList(pfts.size());
		for (int i = 0; i < pfts.size(); i++){
			pointList[i] = pfts[i].pointPos;
		}
		
		pointsBelong.resize(pfts.size());
		closetNextPoint.resize(pfts.size());
		for (size_t i = 0; i < pfts.size(); i++)
		{
			int minContourId = -1;
			float minDistance = 1000000;
			pointFromTwoSource nowPoint = pfts[i];
			OpenMesh::Vec3f closetPos;
			//cout << "i = " << i << ":";
			for each (contour *next in nextContours)
			{
				for (size_t idx = 0; idx<next->pfts.size(); idx++)
				{
					OpenMesh::Vec3f p1 = nowPoint.pointPos;
					OpenMesh::Vec3f p2 = next->pfts[idx].pointPos;
					OpenMesh::Vec3f p3 = next->pfts[(idx + 1) % next->pfts.size()].pointPos;
					float time= calEdgeTime(p1, p2, p3);
					float distance = (p2 - p1).length();
					OpenMesh::Vec3f nowPos = p2;
					if (time >= 0 && time <= 1){
						nowPos = p2 + (p3 - p2)*time;
						distance = (nowPos - p1).length();
					}
					
					if (distance<minDistance)
					{
						minContourId = next->id;
						minDistance = distance;
						closetPos = nowPos;
					}
				}
			}
			pointsBelong[i] = minContourId;
			closetNextPoint[i] = pfts[i].pointPos;//closetPos;
		}
	}
	void calCurvature(){
		curvatures.resize(pfts.size());
		for (size_t i = 0; i < pfts.size(); i++)
		{
			int last_i = i - 1;
			if (last_i < 0){
				last_i = last_i + pfts.size();
			}
			int next_i = (i + 1) % pfts.size();
			OpenMesh::Vec3f p1 = pfts[last_i].pointPos;
			OpenMesh::Vec3f p2 = pfts[i].pointPos;
			OpenMesh::Vec3f p3 = pfts[next_i].pointPos;
			curvature3D curvature(glm::vec3(p1[0], p1[1], p1[2]), glm::vec3(p2[0], p2[1], p2[2]), glm::vec3(p3[0], p3[1], p3[2]));
			curvatures[i] = curvature.getResult();
			cout << "算出來曲率為:" << curvatures[i] << endl;
		}
	}
};
vector<vector<Vec3f_angle>> combineSample(vector<contour*> contourUnion, int sampleNum);
vector<vector<OpenMesh::Vec3f>> combineSample_withoutAngle(vector<contour*> contourUnion, int sampleNum);
struct polarPoint{
	int vertexIdx;
	OpenMesh::Vec3f euclideanPos;
	float distance;
	float angle;
	polarPoint(int idx,OpenMesh::Vec3f ecdpos,float U,float theta){
		vertexIdx = idx;
		euclideanPos = ecdpos;
		distance = U;
		angle = theta;
		
	}
	OpenMesh::Vec2f toEuclid2d(){
		return OpenMesh::Vec2f(distance*cosf(angle),distance*sinf(angle));
	}
};
struct adapter_3dto2d{
	glm::mat2x3 matrix;
	static glm::mat3x4 solveMat(glm::mat3x4 mat){
		//cout << "輸出solveMat開始:" << endl;
		for (int row = 0; row < 2; row++){//不可靠的排序,通過0值的說少為權重決定運算順序
			vector<int> weights;//記錄權重,由於每個迴圈會封掉一個row,要計算權重數量會逐漸減少
			for (size_t i = row; i < 3; i++)
			{
				int weight = 0;
				if (mat[i][row] == 0){//如果關鍵數值為0那就沒什麼好談了
					weight = -1;
				}
				else{
					for (size_t j = 0; j < 3; j++)
					{
						if (mat[i][j] == 0)
							weight++;
					}
				}
				weights.push_back(weight);
			}
			int max_weight_index = 0;//找到權重最大的row
			for (int i = 0; i < weights.size(); i++){
				if (weights[i]>weights[max_weight_index])
				{
					max_weight_index = i;
				}
			}
			if (row + max_weight_index != row){//如果最大權重的row不是當前row,則交換這兩個row的值
				for (size_t i = 0; i < 4; i++)
				{
					float temp = mat[row][i];
					mat[row][i] = mat[row + max_weight_index][i];
					mat[row + max_weight_index][i] = temp;
				}
			}
		}
		//前處理完成,開始解矩陣
		for (size_t i = 0; i < 3; i++)
		{
			mat[i] /= mat[i][i];//歸一化第i行,例如 mat[1]: 4 2 0 1=>mat[1]: 2 1 0 1/2
			glm::vec4 current = mat[i];
			//使其他行第i個元素歸零
			for (size_t j = 0; j < 3; j++)
			{
				if (j!=i)
					mat[j] -= current*mat[j][i];
			}
		}
		return mat;
	}
	void init(glm::mat3x4 mat_x, glm::mat3x4 mat_y){
		//cout << "solve result_x:";
		glm::mat3x4 result_x = solveMat(mat_x);
		glm::vec3 row1(result_x[0][3] / result_x[0][0], result_x[1][3] / result_x[1][1], result_x[2][3] / result_x[2][2]);
		//cout << "solve result_y:";
		glm::mat3x4 result_y = solveMat(mat_y);
		glm::vec3 row2(result_y[0][3] / result_y[0][0], result_y[1][3] / result_y[1][1], result_y[2][3] / result_y[2][2]);
		matrix = glm::mat2x3(row1, row2);
	}
	adapter_3dto2d(glm::vec3 point3dx3[3], glm::vec2 point2dx3[3]){
		glm::mat3x4 mat_x(point3dx3[0][0], point3dx3[0][1], point3dx3[0][2], point2dx3[0][0],
			point3dx3[1][0], point3dx3[1][1], point3dx3[1][2], point2dx3[1][0],
			point3dx3[2][0], point3dx3[2][1], point3dx3[2][2], point2dx3[2][0]);
		glm::mat3x4 mat_y(point3dx3[0][0], point3dx3[0][1], point3dx3[0][2], point2dx3[0][1],
			point3dx3[1][0], point3dx3[1][1], point3dx3[1][2], point2dx3[1][1],
			point3dx3[2][0], point3dx3[2][1], point3dx3[2][2], point2dx3[2][1]);
		init(mat_x, mat_y);
	}
	adapter_3dto2d(float point3dx3[3][3], float point2dx3[3][2]){
		glm::mat3x4 mat_x(point3dx3[0][0], point3dx3[0][1], point3dx3[0][2], point2dx3[0][0],
						  point3dx3[1][0], point3dx3[1][1], point3dx3[1][2], point2dx3[1][0],
						  point3dx3[2][0], point3dx3[2][1], point3dx3[2][2], point2dx3[2][0]);
		glm::mat3x4 mat_y(point3dx3[0][0], point3dx3[0][1], point3dx3[0][2], point2dx3[0][1],
						  point3dx3[1][0], point3dx3[1][1], point3dx3[1][2], point2dx3[1][1],
						  point3dx3[2][0], point3dx3[2][1], point3dx3[2][2], point2dx3[2][1]);
		init(mat_x, mat_y);
	}

	float* to2dPoint(float* point3d){
		glm::vec3 origen(point3d[0], point3d[1], point3d[2]);
		glm::vec2 result = origen*matrix;
		float arr[2] = {result[0],result[1]};
		return arr;
	}
	glm::vec2 to2dPoint(glm::vec3 point3d){
		return point3d*matrix;
	}
};
struct coverAxis{
	glm::vec3 axis[3];
	glm::vec3 tempAxis;
	glm::vec3 center;
	coverAxis(){

	}
	coverAxis(glm::vec3 oripoint,glm::vec3 axis_x, glm::vec3 axis_y, glm::vec3 axis_z){
		axis[0] = glm::normalize(axis_x);
		axis[1] = glm::normalize(axis_y);
		axis[2] = glm::normalize(axis_z);
		center = oripoint;
	}
	glm::vec3 cover(glm::vec3 input){
		glm::vec3 offset = input - center;
		float x = glm::dot(offset, axis[0]);
		float y = glm::dot(offset, axis[1]);
		float z = glm::dot(offset, axis[2]);

		return glm::vec3(x, y, z);
	}
};
struct pathRequestPoint{
	int contour_index[2];
	float contour_angle[2];
	float y;
	OpenMesh::Vec3f pointPos;
	int pixelId = -1;
	pathRequestPoint(int contourIndex1, float angle1, int contourIndex2, float angle2, float at_y,int pId){
		contour_index[0] = contourIndex1;
		contour_index[1] = contourIndex2;
		contour_angle[0] = angle1;
		contour_angle[1] = angle2;
		y = at_y;
		pixelId = pId;
	}

};
class polarMap
{
public:
	vector<contour*> *contours = NULL;
	set<int> debugVertexIndex;
protected:
	MyMesh mesh;

	virtual OpenMesh::ArrayKernel getMesh();
	int poleIdx;
	int zero_anix_vertexIdx;

	void initStartRing();
	bool inited = false;
	OpenMesh::Vec3f baseDirection;
	void calU(int vIdx, vector<MyMesh::VertexHandle> &candidates,bool debug);
	void calU(int vidx, vector<MyMesh::VertexHandle> &candidates, int debugIdx);
	float computeDistance(int vIdx);
	//void calTheta(int vIdx);
	MyMesh::VertexHandle popSmallestNode(vector<MyMesh::VertexHandle> &list);
	vector<MyMesh::VertexHandle> neighbour(MyMesh::VertexHandle vertex);
	draw2dConnectMap *painter;
	int farthestIdx = -1;
	//找出0度軸需要用到的資料
	virtual void findExtremumVertexs();//初始化LocalExtremum_idx的方法
	void removeIncorrectExtremumVertexs();//此方法需要用到contours,所以必須要等getContourLine被呼叫后才能使用
private: void exploreVertex(int nowVIdx, vector<int> &regionVertexIdx, set<int> &boundaryContourIdx,bool debug);//removeIncorrectExtremumVertexs使用的遞迴方法,目的是找到被等高線圍起來的區域(region)內的所有頂點(vertex),boundaryContourIdx是在擴張過程中遇到的等高線的索引值,只有被一拳等高線環繞的局部極值點才是真的
		 vector<contour*> explorePossibleBoundary;
		 int exploreOrigenVIdx = -1;
		 
public: vector<int> localExtremum_idx;//局部極值點,包含最遠的點
		void exploreFace(int contourIdx, int nowFaceIdx, MyMesh::HalfedgeHandle *fromEdge, vector<int> &allFaceIdx);
	//-------------------------------
public:
	//點的資料,包括距離(U)和角度(theta)
	virtual polarPoint getPointFrom(int vertexIndex);
	virtual int pointNum();
	int totalLevel=-1;
	float intervalBtwLevels=-1;
	float *U;
	float *theta;
	float assign_maxPercent=0.9f;
	float assign_minPercent = 0.1f;
	int assign_segmentNum = 16;
	float assigh_distEnergy_weight = 8.0f;
	const float getPathYAsZeroThreshold = 0.00001f;
	const float matchEndPointThreshold = 0.001f;
	//-------------------------------
	const float BIG_NUMBER = 1000000;
	//debug方法,用於印出點資料
	polarMap(MyMesh mesh,int poleIndex);
	string meshPath;
	void showPoints(int range);
	void showPoints(vector<int> vList);
	void debugDist(int rings);
	//-------------------------------
	~polarMap();
	int debugCounter = 0; 
	//把參數化坐標畫在圖片上的方法,同樣是用來除錯的方法
	virtual void createMap();
	virtual void drawMeshToImage();
	static void drawContourToImage(vector<PointCurve> contours);
	//-------------------------------
	virtual int getFasthestIdx();//最遠的點的索引值
	const float dropContourThreshold = 0;//0.1;//總長度(length())小於dropContourThreshold的contour會被捨棄掉
	vector<PointCurve> getContourLine(int line_number,vector<int> axis);//已被廢棄
	vector<PointCurve> getContourLine(int line_number);
	float AssignIgnorePartAngleThreshold = 0.01;//如果做mapping的起始段太過短,則mapping的結果不可信,當該段所佔的角度百分比小於AssignIgnorePartAngleThreshold時,該段不納入考慮
	void assignDivergence();
	vector<vector<OpenMesh::Vec3f>> dividePointSet(MyMesh mesh, vector<contour> pointSet);//已被廢棄
	vector<vector<OpenMesh::Vec3f>> dividePointSet(MyMesh mesh, vector<pointFromTwoSource> pointSet);//已被廢棄
	vector<vector<Vec3f_angle>> dividePointSet(MyMesh mesh, vector<pointFromTwoSource> pointSet, vector<int> axis);//已被廢棄
	vector<vector<OpenMesh::Vec3f>> getIntersections(vector<PointCurve> contours,vector<int> pointNum);//已被廢棄
	//void calAngleUseAxis(vector<PointCurve> &lines, vector<int> baseAxis);//
	bool Image_show_vertex_idx=false;
	//vector<OpenMesh::Vec3f> pathHistory
	
protected:
	void tryMeetNextContour(MyMesh::VertexHandle now_vh, contour current, vector<contour*> candiate, set<int> &meetContourIndex,set<int> &vertexHistory);//遞迴方法,試著從current中出發在,遇到candiate中的一個contour,回傳值為是遇到的所有的contour在candiate中的索引值
	void calRingVertex(vector<OpenMesh::VertexHandle> &candiates,float* &distArray,int* &lastVertexArray);
	void calRingVertex(vector<OpenMesh::VertexHandle> &candiates, float* &distArray, int* &lastVertexArray, vector<int> vidx_domain);
	float farestFromVertex = 0;
	//int getFaceIdxAtContourAngle();
	float* distFormVertex;
	geodesicPathFinder pFinder;

	void initPathFinder();


public:	
	geodesicAsynFinder asynFinder;
	vector<vector<pathRequestPoint>> requestTable;
	float pathLengthThreshold = 0.5f;
	virtual vector<vector<int>> calBaseAxis();
	float getPercentageFromVertex(int vidx);
	void showContours();
	static OpenMesh::Vec3f calPercentagetPosAtPath(vector<OpenMesh::Vec3f> path, float percentage);
	vector<contour*> getPossibleBoundary(float U);
	vector<contour*>  getPossibleBoundary(float U,int extraLayer);//額外增加一個等級的boundary
	vector<OpenMesh::Vec3f> getPointBetween(int contour1_index,float at_angle, float percentage);//已經廢棄
	vector<OpenMesh::Vec3f> getPointsBtw(int contour1_index, float at_angle);
	vector<OpenMesh::Vec3f> getPointsBtw(int contour1_index, float at_angle, bool &fail);
	vector<OpenMesh::Vec3f> getPointsBtw(int contour1_index, float at_angle1, float at_angle2, bool &fail);
	vector<OpenMesh::Vec3f> getPointsBtw(int contour1_index, int next_contour_index, float at_angle1, float at_angle2, bool &fail);
	vector<OpenMesh::Vec3f> getPointsBtw_gdcPath(int contour1_index, float at_angle1, float at_angle2, vector<int> &shortestPath);
	vector<OpenMesh::Vec3f> getgdcPointsBtw(int contour1_index, float at_angle1, float at_angle2);
	vector<OpenMesh::Vec3f> getgdcPointsBtw(int contour1_index, int next_contour_index, float at_angle1, float at_angle2);
	
	vector<OpenMesh::Vec3f> pathRecord;
	OpenMesh::Vec3f getPointAtBtw(int contour1_index, float at_angle,float percent);//已經棄用
	OpenMesh::Vec3f getPointAtBtw(int contour1_index, float at_angle, float percent,bool &fail);//已經棄用
	OpenMesh::Vec3f getPointAtBtw(int contour1_index, float at_angle1, float at_angle2,float percent, bool &fail);
	OpenMesh::Vec3f getPointAtBtw(int contour1_index, int next_contour_index, float at_angle1, float at_angle2, float percent, bool &fail);
	OpenMesh::Vec3f getgdcPointAtBtw(int contour1_index, float at_angle1, float at_angle2, float percent);
	void recordGdcPointLineBtw(vector<pathRequestPoint> requstLine);
	vector<OpenMesh::Vec3f> getAxisAt(int contour1_index, float at_angle,float vecLen);//已經棄用
	pointAtContour closetPointInNext(int contourIdx, float angle);//已經棄用
	vector<vector<pathRequestPoint>> RespondsGdcPoints();
	const int RESPOND_CYCLE = 100;//多少次尋路之後重新運行一次asynFinder.doit();
	coverAxis lastCover;

};

struct anglepair{
	float ph_ki;
	float ph_ij;
	anglepair(float a1, float a2){
		ph_ki = a1;
		ph_ij = a2;
	}
};

 float angleOf(OpenMesh::Vec3f v1, OpenMesh::Vec3f v2);
 anglepair cal_phki_phij(OpenMesh::Vec3f ek, OpenMesh::Vec3f ej, float Uk, float Uj);
 OpenMesh::Vec2f normalizePos(float* border, OpenMesh::Vec2f ori);

 struct PlaneFit
 {
	 float A, B, C, D;
	 bool is_fit = false;
	 OpenMesh::Vec3f first;
	 OpenMesh::Vec3f second;
	 static double dot(OpenMesh::Vec3f v1, OpenMesh::Vec3f v2){
		 return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	 }
	 OpenMesh::Vec3f normal(){
		 return OpenMesh::Vec3f(A, B, C);
	 }
	 PlaneFit(std::vector<OpenMesh::Vec3f> points){
		 is_fit = false;
		 int numPoints = points.size();
		 if (numPoints >= 3){
			 //Compute the mean of the points.
			 OpenMesh::Vec3f mean = OpenMesh::Vec3f(0, 0, 0);
			 for (int i = 0; i < numPoints; ++i){
				 mean = mean + points[i];
			 }
			 mean = mean * (1.0 / numPoints);

			 float xx = 0.0; float xy = 0.0; float xz = 0.0;
			 float yy = 0.0; float yz = 0.0; float zz = 0.0;

			 //Compute the covariance matrix of the points.
			 for (int i = 0; i < numPoints; ++i){
				 OpenMesh::Vec3f diff = points[i] - mean;
				 xx += diff[0] * diff[0];
				 xy += diff[0] * diff[1];
				 xz += diff[0] * diff[2];
				 yy += diff[1] * diff[1];
				 yz += diff[1] * diff[2];
				 zz += diff[2] * diff[2];
			 }

			 xx /= numPoints;
			 xy /= numPoints;
			 xz /= numPoints;
			 yy /= numPoints;
			 yz /= numPoints;
			 zz /= numPoints;

			 OpenMesh::Vec3f weighted_dir(0.0, 0.0, 0.0);
			 OpenMesh::Vec3f axis_dir;
			 double weight;

			 //X-axis
			 double det_x = yy*zz - yz*yz;
			 axis_dir = OpenMesh::Vec3f(det_x, xz*yz - xy*zz, xy*yz - xz*yy);
			 weight = det_x*det_x;
			 if (dot(weighted_dir, axis_dir) < 0.0){
				 weight = -weight;
			 }
			 weighted_dir = axis_dir * weight + weighted_dir;
			 //Y-axis
			 double det_y = xx*zz - xz*xz;
			 axis_dir = OpenMesh::Vec3f(xz*yz - xy*zz, det_y, xy*xz - yz*xx);
			 weight = det_y*det_y;
			 if (dot(weighted_dir, axis_dir) < 0.0){
				 weight = -weight;
			 }
			 weighted_dir = axis_dir * weight + weighted_dir;
			 //Z-axis
			 double det_z = xx*yy - xy*xy;
			 axis_dir = OpenMesh::Vec3f(xy*yz - xz*yy, xy*xz - yz*xx, det_z);
			 weight = det_z*det_z;
			 if (dot(weighted_dir, axis_dir) < 0.0){
				 weight = -weight;
			 }
			 weighted_dir = axis_dir * weight + weighted_dir;
			// cout << "weight_dir:" << weighted_dir;
			 weighted_dir = weighted_dir.normalize();
			 //cout << "after normalize:" << weighted_dir<<endl;
			 if (isfinite(weighted_dir[0]) && isfinite(weighted_dir[1]) && isfinite(weighted_dir[2])){
				 first = mean;
				 A = second[0] = weighted_dir[0];
				 B = second[1] = weighted_dir[1];
				 C = second[2] = weighted_dir[2];
				 D = -dot(mean, weighted_dir);
				 is_fit = true;
			 }
		 }

		 first = OpenMesh::Vec3f(0, 0, 0);
		 second = OpenMesh::Vec3f(0, 0, 0);

	 }

 };
 class polarMap_angleProjection:public polarMap{
 public:
	 polarMap_angleProjection(MyMesh mesh, int poleIndex);
	 virtual vector<vector<int>> calBaseAxis();
 };
 bool contain(vector<int> list, int elem);