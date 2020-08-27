#pragma once
#include <stdlib.h>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
using namespace std;
//用來處理使用者輸入圖元的類別
struct percentPoint{
	bool inited = false;
	float x=-1;
	float y=-1;
	percentPoint(){

	}
	percentPoint(float px, float py){
		inited = true;
		x = px;
		y = py;
	}
};
struct pixelPoint
{
	int id = -1;
	int x = -1;
	int y = -1;
	float percentage_x;
	float percentage_y;
	bool intersection = false;
	vector<pixelPoint*> lastPoints;
	vector<pixelPoint*> nextPoints;
	pixelPoint(){

	}
	pixelPoint(int posx, int posy,float value_x,float value_y,int newId){
		x = posx;
		y = posy;
		percentage_x = value_x;
		percentage_y = value_y;
		id = newId;
	}
};
struct  col
{
	int height;
	pixelPoint **points;
	col(){
		height = -1;
		points = NULL;
	}
	col(pixelPoint** list, int h){
		height = h;
		points = new pixelPoint*[h];
		for (int i = 0; i < h; i++){
			points[i] = list[i];
		}
		//points = list;
	}
};
struct graphLine{
	vector<pixelPoint> *points;
	vector<graphLine*> leftSide;
	vector<graphLine*> rightSide;
	vector<int> controlPointIndexs;
	vector<pixelPoint*> keyPoints;
	int matWidth = 0;
	int matHeight =0;
	graphLine(){
	}
	graphLine(vector<pixelPoint> *line,int mwidth,int mheight){
		points = line;
		matWidth = mwidth;
		matHeight = mheight;
		for (int i = 0; i < line->size(); i++){
			pixelPoint point = points->operator[](i);
			if (point.percentage_y == 0 || point.percentage_x == 0 || point.percentage_y == 1 || point.percentage_x == 1){
				controlPointIndexs.push_back(i);
				keyPoints.push_back(&points->operator[](i));
				//cout << "keypoint address:" << keyPoints[0] << ", &points[" << i << "]:" << &points[i];
			}
		}
		cout << endl;
	}
	pixelPoint operator[](int index){
		return points->operator[](index);
	}
	pixelPoint leftPoint(){
		return points->front();
	}
	pixelPoint rightPoint(){
		return points->back();
	}
	static bool pixelEqual(pixelPoint p1, pixelPoint p2){
		return p1.x == p2.x&&p1.y == p2.y;
	}
	void tryConnect(graphLine *another){
		if (pixelEqual(leftPoint(), another->leftPoint()) || pixelEqual(leftPoint(), another->rightPoint()))
			leftSide.push_back(another);
		if (pixelEqual(rightPoint(), another->leftPoint()) || pixelEqual(rightPoint(), another->rightPoint()))
			rightSide.push_back(another);
	}
	void finshConnect(){
		if (leftSide.size()>0)
		{
			bool hasLeft = false;
			for each (int  idx in controlPointIndexs)
			{
				if (idx == 0){
					hasLeft = true;
					break;
				}
			}
			if (!hasLeft)
			{
				controlPointIndexs.insert(controlPointIndexs.begin(), 0);
			}
		}
		if (rightSide.size() > 0){
			bool hasRight = false;
			for each (int idx in controlPointIndexs)
			{
				if (idx == points->size() - 1){
					hasRight = true;
					break;
				}
			}
			if (!hasRight){
				controlPointIndexs.push_back(points->size() - 1);
			}
		}
	}
	void updatePixelPoint(pixelPoint &point, float offset_x, float offset_y){
		point.percentage_x += offset_x;
		point.percentage_y += offset_y;
		point.x += offset_x*matWidth;
		point.y += offset_y*matHeight;
	}
	void deform(pixelPoint *controlPoint, percentPoint tragetPoint){
		cout << "controlPoint:(" << controlPoint->percentage_x << "," << controlPoint->percentage_y << ")";
		cout << "tragetPoint:(" << tragetPoint.x << "," << tragetPoint.y << ")";
		cout << "controlPointIndexs:";
		for each (int index in controlPointIndexs)
		{
			cout << index << ",";
		}
		cout << "points:";
		/*for each (pixelPoint point in points)
		{
			cout << "(" << point.percentage_x << "," << point.percentage_y << ");";
		}*/
		for (size_t i = 0; i < points->size(); i++)
		{
			if (&points->operator[](i) == controlPoint){
				int cptIdx = i;
				cout << "cptIdx:" << cptIdx;
				for (int j = 0; j < controlPointIndexs.size(); j++){
					if (controlPointIndexs[j] == cptIdx){
						cout << "j=" << j;
						int last_cpIdx;
						int next_cpIdx;
						if (j==0)
						{
							last_cpIdx = cptIdx;
						}
						else{
							last_cpIdx = controlPointIndexs[j-1];
						}
						
						if (j == controlPointIndexs.size() - 1){
							next_cpIdx = cptIdx;
						}
						else
						{
							next_cpIdx = controlPointIndexs[j + 1];
						}
						cout << "last_cpIdx:" << last_cpIdx << " next_cpIdx:" << next_cpIdx;
						//計算出基準長度,是向左或向右最長的距離
						float length_tolast = 0;
						vector<float> lengthStep_last;
						for (int idx = last_cpIdx; idx < cptIdx; idx++){
							float len = sqrt(powf(points->operator[](idx + 1).percentage_x - points->operator[](idx).percentage_x, 2) + powf(points->operator[](idx + 1).percentage_y - points->operator[](idx).percentage_y, 2));
							length_tolast += len;
							cout << idx << "->" << idx-1 <<" len:"<<len<< ", ";
							lengthStep_last.insert(lengthStep_last.begin(), len);
						}

						cout << endl;
						cout << "length_tolast:" << length_tolast << ";"<<endl;
						float length_tonext = 0;
						vector<float> lengthStep_next;
						for (int idx = cptIdx; idx < next_cpIdx; idx++){
							float len = sqrt(powf(points->operator[](idx + 1).percentage_x - points->operator[](idx).percentage_x, 2) + powf(points->operator[](idx + 1).percentage_y - points->operator[](idx).percentage_y, 2));
							length_tonext += len;
							cout << idx << "->" << idx - 1 << " len:" << len << ", ";
							lengthStep_next.push_back(len);
						}
						cout << "length_tonext:" << length_tonext << ";"<<endl;
						float length_longest = length_tolast;
						if (length_tonext > length_tolast){
							length_longest = length_tonext;
						}
						cout << "length_longest:" << length_longest << ";";
						float deform_x = tragetPoint.x - points->operator[](cptIdx).percentage_x;
						float deform_y = tragetPoint.y - points->operator[](cptIdx).percentage_y;
						//從當前點出發,隨距離遞減地傳遞變形
						updatePixelPoint(points->operator[](cptIdx), deform_x, deform_y);
						float deform_energy = 1;
						int countl = 0;
						for (int l = cptIdx - 1; l >= last_cpIdx; l--){
							//float len = sqrt(powf(points[l].percentage_x - points[l+1].percentage_x, 2) + powf(points[l].percentage_y - points[l+1].percentage_y, 2));
							deform_energy -= lengthStep_last[countl++] / length_longest;
							cout <<"l="<<l<<"deform_energy:"<<deform_energy;
							//cout << l << "->" << l + 1 << " len:" << len << ", ";
							updatePixelPoint(points->operator[](l), deform_x*deform_energy, deform_y*deform_energy);
						}
						cout << endl;
						deform_energy = 1;
						int countn = 0;
						cout << "cptIdx:" << cptIdx << "next_cptidx:" << next_cpIdx << "points.size():" << points->size()<<endl;
						for (int n = cptIdx + 1; n <= next_cpIdx; n++){
							//float len = sqrt(powf(points[n].percentage_x - points[n -1].percentage_x, 2) + powf(points[n].percentage_y - points[n-1].percentage_y, 2));
							deform_energy -= lengthStep_next[countn++] / length_longest;
							cout << "n=" << n;
							cout << "points[n]:(" << (*points)[n].percentage_x << "," << (*points)[n].percentage_y << "),";
							updatePixelPoint((*points)[n], deform_x*deform_energy, deform_y*deform_energy);
						}
						break;
					}
					
				}
				break;
			}
		}
	}
	
};
struct  controlPoint
{
	vector<pixelPoint*>points;
	vector<graphLine*> lines;
	controlPoint* matchPoint = NULL;
	void deform(percentPoint pos){
		for (int i = 0; i < lines.size(); i++){
			lines[i]->deform(points[i], pos);
		}
	}
};
struct controlGraph
{
	vector<graphLine*> mainGraphs;
	vector<controlPoint*> contorlPoints;
	cv::Mat *Image;
};
struct controlPointGroup{
	vector<controlPoint*> group1;
	vector<controlPoint*> group2;
};
class imageConvert
{
public:
	map<int, map<int, int>> connectDirect;
	vector<pair<int, int>> extraStartDirect;
	bool inExtraStartDirect(int id1,int id2);
	bool useSkeleton = false;
	int skeletonThreshold = 50;
	int mergeYThreshold;
	float mergeVerticalThreshold = 0.05f;
	int width;
	col *datas;
	imageConvert(string fileName);
	void calSkeleton();
	void createData();
	~imageConvert();
	void debugDatas();
	vector<vector<pixelPoint>> sample(int sampleNum);
	
	//vector<percentPoint> keyConnectPoints;
	void showGraphicConnect(vector<vector<pixelPoint>> &graph);
	static void drawBaseGraph(cv::Mat *image,vector<graphLine*> graph);
	static void drawBaseGraph(cv::Mat *image, vector<graphLine*> graph, int width, int height);
	static void drawNextGraph(cv::Mat *image, vector<graphLine*> graph, float start_percentX, int width, int height);
	static void drawControlPoints(cv::Mat *image, vector<controlPoint*> ckpts);
	static void drawGraphPoints(cv::Mat *image, vector<graphLine*> graph);
	bool checkConnectOk = false;
	vector<graphLine*> lineGraph;
	float Shift_percent = 0;
	int totalPixelPointNum = 0;
protected:
	vector<pixelPoint*> flatArray;
	vector<pixelPoint*> intersectionPoint;
	vector<pixelPoint*> findLastConnectPoint(int posx,int posy);
	//void findKeyConnectPoints();
	controlPointGroup horizontalMatch(vector<controlPoint*> controls,float margeThreshold);
	bool horizonalIsOk(controlPointGroup match);
	controlPointGroup verticalMatch(vector<controlPoint*> controls, float margeThreshold,bool &success);
	cv::Mat image;
	cv::Mat skeletonImg;
	
};

