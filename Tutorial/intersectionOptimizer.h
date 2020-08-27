#pragma once
#include "polarMap.h";
//嘗試通過貪婪演算法減小變形程度,現在減小變形度的功能已經作廢,因為減小變形度導致了anglelist的不均勻,這個不均勻會逐等高線遞增,最後結果反而更加不好
struct rectGrid
{
	OpenMesh::Vec3f points[4];
	float area = -1;
	float weight = -1;
	float minWeight = -1;
	OpenMesh::Vec3f minDeformPoint[4];
	float minAnglePair[2];
	float minAngleCrossZero[2];
	rectGrid(OpenMesh::Vec3f p11, OpenMesh::Vec3f p12, OpenMesh::Vec3f p21, OpenMesh::Vec3f p22, float w,float min_w){
		points[0] = p11;
		points[1] = p12;
		points[2] = p21;
		points[3] = p22;


		const int x = 0;
		const int y = 1;
		const int z = 2;
		OpenMesh::Vec3f p1 = p12 - p11;
		OpenMesh::Vec3f p2 = p21 - p11;
		OpenMesh::Vec3f faceVec1(p1[y]*p2[z] - p2[y]*p1[z], p1[z]*p2[x] - p2[z]*p1[x], p1[x]*p2[y]-p2[x]*p1[y]);
		
		p1 = p12 - p22;
		p2 = p21 - p22;
		OpenMesh::Vec3f faceVec2(p1[y]*p2[z] - p2[y]*p1[z], p1[z]*p2[x] - p2[z]*p1[x], p1[x]*p2[y] - p2[x]*p1[y]);
		area = faceVec1.length()*0.5f+faceVec2.length()*0.5f;
		weight = w;
		minWeight = min_w;
	}
	rectGrid(){

	}
	void updatePartsPoints(OpenMesh::Vec3f *p11, OpenMesh::Vec3f *p12, OpenMesh::Vec3f *p21, OpenMesh::Vec3f *p22,float w){
		if (p11 != NULL){
			points[0] = *p11;
		}
		if (p12 != NULL){
			points[1] = *p12;
		}
		if (p21 != NULL){
			points[2] = *p21;
		}
		if (p22 != NULL){
			points[3] = *p22;
		}
		const int x = 0;
		const int y = 1;
		const int z = 2;
		OpenMesh::Vec3f p1 = points[1]-points[0];
		OpenMesh::Vec3f p2 = points[2]-points[0];
		OpenMesh::Vec3f faceVec1(p1[y] * p2[z] - p2[y] * p1[z], p1[z] * p2[x] - p2[z] * p1[x], p1[x] * p2[y] - p2[x] * p1[y]);

		p1 = points[1]-points[3];
		p2 = points[2]-points[3];
		OpenMesh::Vec3f faceVec2(p1[y] * p2[z] - p2[y] * p1[z], p1[z] * p2[x] - p2[z] * p1[x], p1[x] * p2[y] - p2[x] * p1[y]);
		area = faceVec1.length()*0.5f + faceVec2.length()*0.5f;
		weight = w;

	}
};
struct floatPair{
	float key;
	float value;
	bool flag = false;
	floatPair(float k, float v){
		key=k;
		value=v;
	}
	floatPair(){
		key = -1;
		value = -1;
	}
};
struct floatGroup
{
	float key;
	vector<floatPair> values;
	bool flag = false;
	floatGroup(){

	}
	floatGroup(float k){
		key = k;
	}
};
struct gridGroups{
	int index;
	vector<floatGroup> groups;
	gridGroups(int idx){
		index = idx;
	}
	gridGroups(){
		index = -1;
	}
};
struct  deformLabel
{
	Vec3f_angle samplePoint;
	float energy;
	deformLabel(){

	}
	deformLabel(Vec3f_angle pt, float e){
		samplePoint = pt;
		energy = e;
	}
};
struct angleArea
{
	float_pair angleBoundary;
	vector<float> fixGridAngles;//不能被改動的格子的左右邊的角度,比如有兩個格子就會有3個值來表示   ...| grid1 | grid2 |...
	float_pair minAngleBoundary(){
		return float_pair(fixGridAngles.front(), fixGridAngles.back());
	}
	angleArea(){

	}
	angleArea(float_pair boundary){
		angleBoundary = boundary;
	}
	angleArea(float_pair boundary,vector<float> angles){
		angleBoundary = boundary;
		fixGridAngles = angles;
	}

};
vector<float> assignGridWithInitArea(vector<angleArea> initArea);
vector<float> newAngleList(int elemNum);
vector<float> newAngleList(float startAngle, int elemNum);
class intersectionOptimizer//嘗試用貪婪演算法減少變形量,最後以失敗告終//已經棄用
{
public:
	float minJumpPercentage = 2.0f;//初始化中取樣能無視對應角度選擇最近點相對於interval的最大範圍
	float minIntervalPercentage = 0.6f;//solveOnce中grid可以被擠壓到的相對interval最小的比例
	float maxIntervalPercentage = 1.5f;
	float coeff_edge=1000;
	vector<gridGroups> deformRecords;
	vector<rectGrid> gridList;
	vector<floatPair> angleMappings;
	void printBestGridInf();
	vector<OpenMesh::Vec3f> showBestGridPoints();
	contour contourList[2];
	int girdNum = -1;
	int sampleRate = 1000;//在計算權重時採樣1000個點;
	intersectionOptimizer();
	intersectionOptimizer(contour last, contour next,int num);
	intersectionOptimizer(contour last, contour next,vector<float> angleList);
	intersectionOptimizer(contour last, contour next, vector<float> angleList, bool negAngel);
	~intersectionOptimizer();
	float calDeform(OpenMesh::Vec3f p11,OpenMesh::Vec3f p12,OpenMesh::Vec3f p21,OpenMesh::Vec3f p22);
	float calLengrhDeform(float avgLength, float angle1, float angle2);
	void solve(int iter_num);
	void solveOnce();
};

struct offsetInf{
	int origenContourIndex = -1;
	vector<float> origenAngleList;//用來記錄初始一環的angle list
	int subLevel = 0;//用來記錄當前等高線是位於其源頭等高線后的第幾層等高線,源頭等高線的subLevel為0
	float totalOffset = 0;//totalOffset為累積的offset_X,源頭等高線后每一級就累積一個offset_X
};
class optimizerGroup{
protected:
	int elemPerRing;
	vector<contour*> contours;
	vector<vector<float>> angleListForContours[2];//每個angleListForContour[0][x]和angleListForContour[1][x]屬於每個ring,這個level由x = contour.id-1作為起始contour,angleListForContour[0][x]表示每個lelem在起始contour上的起點的角度
	//angleListForContour[1][x]為每個elem在contour下一條contour上的起始角度,當然這之適用與僅有一個nextContour的情況,多於1個或者0個的時候實際上不會用到angleListForContour[1]
	
	int OPTIMIZER_SOTP_GRID_NUM =1000;
	float contourInterval = 0;
public:
	optimizerGroup(vector<contour*> contours,int elemPerRing);
	optimizerGroup(vector<contour*> contours, int elemPerRing, float interval);
	//void solve();
	void stage_solve(contour* stageBegin);//舊的stage_solve,已經棄用
	void stage_solve(contour* stageBegin,float elem_startPoint_shift_percent);
	void printAngleListForContours();
	void printAngleListForContours_stage(vector<contour*> contours,int beginIndex);
	vector<float> getAngleListForContour(int index,bool contour1);
	void calAllDeformEnergy();
	int sampleRate = 1000;
	vector<vector<deformLabel>> deformEnergys;
	vector<Vec3f_angle*> sampleRecord;
	vector<int> contourHistory;
	vector<offsetInf> origenInf;
	float getXat(int contour_index,int elem_index,float x);
};
void bubbleSort(vector<floatPair> &marray);