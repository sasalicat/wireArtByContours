#pragma once
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <iostream>
#include "OpenMesh\Core\IO\MeshIO.hh"
#include "OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh"
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
#include <vector>

using namespace std;
//用於提取離散線段的IS_Match featrue matrix
class IS_FeatureMatrix
{
public:
	IS_FeatureMatrix();
	IS_FeatureMatrix(vector<OpenMesh::Vec3f> c);
	IS_FeatureMatrix(vector<OpenMesh::Vec3f> c, int startIndex, int endIndex);
	~IS_FeatureMatrix();
	const int offset = 3;
	virtual void init(vector<OpenMesh::Vec3f> c);
	virtual void init(vector<OpenMesh::Vec3f> c, int startIndex, int endIndex);
	cv::Mat showAsImage(int multiply);
	cv::Mat matrix;

	
};
vector<cv::Mat> getIntegralImages(cv::Mat ImgN,cv::Mat ImgM);
cv::Mat scaleMat(cv::Mat origenImg, int scale);
void clampIntensityTo01(cv::Mat &image);
cv::Mat ClampIntensityTo01(cv::Mat image);
float calAvgDiff(int s, int t, int l, vector<cv::Mat> IntegImgs);
