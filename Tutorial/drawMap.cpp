#include "stdafx.h"
#include "drawMap.h"




draw2dConnectMap::draw2dConnectMap(int width, int height)
{
	map = cv::Mat(height, width, CV_8UC3, cv::Scalar(255, 255, 255));
	cout << "width:" << map.size().width << "height:" << map.size().height << endl;
	meshColor = cv::Scalar(0,0,0);
	//std::cout << "new map:" << map << std::endl;
}

cv::Point2i draw2dConnectMap::toPixelCoord(cv::Point2f textureCoord_p){
	return cv::Point2i(textureCoord_p.x*map.size().width, textureCoord_p.y*map.size().height);
}
draw2dConnectMap::~draw2dConnectMap()
{
	
}

void draw2dConnectMap::showMap(){
	cv::imshow("now map", map);
	cv::waitKey(0);
}
void draw2dConnectMap::drawTriangle(threePoints triangle){
	cv::line(map, toPixelCoord(triangle[0]), toPixelCoord(triangle[1]),meshColor);
	cv::line(map, toPixelCoord(triangle[1]), toPixelCoord(triangle[2]),meshColor);
	cv::line(map, toPixelCoord(triangle[2]), toPixelCoord(triangle[0]),meshColor);
}
void draw2dConnectMap::drawLine(float x1, float y1, float x2, float y2){
	cv::line(map, toPixelCoord(cv::Point2f(x1, y1)), toPixelCoord(cv::Point2f(x2, y2)),meshColor);
}

void draw2dConnectMap::drawPoint(float x, float y){
	//cv::Point(map,toPoint)
	cv::circle(map,toPixelCoord(cv::Point2f(x,y)),3,cv::Scalar(0,0,255),-1);
}
void draw2dConnectMap::writeTxt(float x, float y, string txt){
	cv::putText(map, txt, toPixelCoord(cv::Point2f(x, y)), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255, 0, 0));
}