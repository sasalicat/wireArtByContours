#pragma once
#include <opencv2/opencv.hpp>
#include <stdio.h>

using namespace std;
//用於在polarMap內將3dmesh畫到圖片里

struct threePoints{
	float value[6];
	threePoints(float* list){
		for (int i = 0; i < 6; i++){
			value[i] = list[i];
		}
	}
	threePoints(float x1, float y1, float x2, float y2, float x3, float y3){
		value[0] = x1;
		value[1] = y1;
		value[2] = x2;
		value[3] = y2;
		value[4] = x3;
		value[5] = y3;
	}
	cv::Point2f operator[](int index){
		if (index >= 3 || index < 0){
			return cv::Point2f(-1, -1);
		}
		else{
			return cv::Point2f(value[2 * index], value[2 * index + 1]);
		}
	}
};
class draw2dConnectMap
{
public:
	cv::Mat map;
	draw2dConnectMap(int map_width, int map_height);
	void showMap();
	void drawTriangle(threePoints triangle);
	cv::Point2i toPixelCoord(cv::Point2f p);
	~draw2dConnectMap();
	cv::Scalar meshColor;
	void drawPoint(float x,float y);
	void drawLine(float x1, float y1, float x2, float y2);
	void writeTxt(float x, float y,string txt);
};

