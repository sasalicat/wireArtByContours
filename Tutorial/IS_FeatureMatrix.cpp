#include "stdafx.h"
#include "IS_FeatureMatrix.h"
using namespace cv;
float dot(OpenMesh::Vec3f v1, OpenMesh::Vec3f v2){
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
IS_FeatureMatrix::IS_FeatureMatrix(){

}

IS_FeatureMatrix::IS_FeatureMatrix(vector<OpenMesh::Vec3f> c)
{
	init(c);
}
IS_FeatureMatrix::IS_FeatureMatrix(vector<OpenMesh::Vec3f> c, int startIndex, int endIndex){
	init(c, startIndex, endIndex);
}
void IS_FeatureMatrix::init(vector<OpenMesh::Vec3f> c, int startIndex, int endIndex){
	int cLength = endIndex - startIndex + 1;
	matrix = Mat(cLength, cLength, CV_32FC1);
	for (int i = startIndex; i <= endIndex; i++){
		for (int j = startIndex; j <= endIndex; j++){
			bool debug = false;
			if (debug){
				cout << "y=" << i << ",x=" << j <<" " ;
			}
			int j_ = j - offset;
			int fixi = i;
			if (i < 0){
				fixi = c.size() + i;
			}
			int fixj = j;
			if (j < 0){
				fixj = c.size() + j;
			}
			int fixj_ = j_;
			if (j_ < 0){
				fixj_ = c.size() + j_;
			}
			if (i == j){
				matrix.at<float>(i-startIndex, j-startIndex) = 0;
			}
			else if (i == j_){
				matrix.at<float>(i-startIndex, j-startIndex) = 0;
			}
			else{
				OpenMesh::Vec3f Pi = c[fixi];
				OpenMesh::Vec3f Pj = c[fixj];

				OpenMesh::Vec3f Pj_ = c[fixj_];
				if (debug){
					cout << "i=" << i << "j=" << j << "j_ =" << j_;
					cout << "Pi:(" << Pi[0] << "," << Pi[1] << "," << Pi[2] << ") Pj:(" << Pj[0] << "," << Pj[1] << "," << Pj[2] << ") Pj_:(" << Pj_[0] << "," << Pj_[1] << "," << Pj_[2] << ")";
				}
				OpenMesh::Vec3f PjPi = Pi - Pj;
				OpenMesh::Vec3f PjPj_ = Pj_ - Pj;

				float cos = dot(PjPi, PjPj_) / (PjPi.length()*PjPj_.length());
				if (debug)
					cout << "cos:" << cos;
				float angle = acosf(cos) / (M_PI);
				if (cos >= 1.0f){
					angle = 0.0f;
				}
				else if (cos <= -1.0f)
				{
					angle = 1.0f;
				}
				if (debug)
					cout << "angle:" << angle << endl;
				matrix.at<float>(i-startIndex, j-startIndex) = angle;
			}
		}
	}
}
void IS_FeatureMatrix::init(vector<OpenMesh::Vec3f> c){
	int cLength = c.size();
	matrix = Mat(cLength, cLength, CV_32FC1);
	for (int i = 0; i < c.size(); i++){
		for (int j = 0; j < c.size(); j++){
			bool debug = false;//i == 4;//(j == 0 && i == 13);
			int j_ = j - offset;
			if (j_ < 0){
				j_ = c.size() + j_;
			}

			if (i == j){
				matrix.at<float>(i, j) = 0;
			}
			else if(i ==j_){
				matrix.at<float>(i, j) = 0;
			}
			else{
				OpenMesh::Vec3f Pi = c[i];
				OpenMesh::Vec3f Pj = c[j];

				OpenMesh::Vec3f Pj_ = c[j_];
				if (debug){
					cout << "i=" << i << "j=" << j << "j_ =" << j_;
					cout << "Pi:(" << Pi[0] << "," << Pi[1] << "," << Pi[2] << ") Pj:(" << Pj[0] << "," << Pj[1] << "," << Pj[2] << ") Pj_:(" << Pj_[0] << "," << Pj_[1] << "," << Pj_[2] << ")";
				}
				OpenMesh::Vec3f PjPi = Pi - Pj;
				OpenMesh::Vec3f PjPj_ = Pj_ - Pj;
				if (debug){
					cout << "PjPi:(" << PjPi[0] << "," << PjPi[1] << "," << PjPi[2] << ") " << "PjPj_:(" << PjPj_[0] << "," << PjPj_[1] << "," << PjPj_[2] << ")";
				}
				float cos = dot(PjPi, PjPj_)/(PjPi.length()*PjPj_.length());
				if (debug)
					printf("cos:%f", cos);//cout << "cos:" << cos;
				float angle = acosf(cos)/(M_PI);
				if (cos >= 1.0f){
					angle = 0.0f;
				}
				else if (cos<=-1.0f)
				{
					angle = 1.0f;
				}
				if (debug)
					printf("angle:%f\n", angle);//cout << "angle:" << angle<<endl;
				matrix.at<float>(i, j) = angle;
			}
		}
	}

	//showAsImage(10);
}

IS_FeatureMatrix::~IS_FeatureMatrix()
{
}
void clampIntensityTo01(cv::Mat &image){
	float maxIntensity = 0;
	for (int y = 0; y < image.size[1]; y++){
		for (size_t x = 0; x <image.size[0]; x++)
		{
			float nowInt = image.at<float>(y, x);
			if (nowInt > maxIntensity){
				maxIntensity = nowInt;
			}
		}
	}
	for (int y = 0; y < image.size[1]; y++){
		for (size_t x = 0; x < image.size[0]; x++)
		{
			image.at<float>(y, x) = image.at<float>(y, x) / maxIntensity;
		}
	}
}
cv::Mat ClampIntensityTo01(cv::Mat image){
	float maxIntensity = 0;
	for (int y = 0; y < image.size[1]; y++){
		for (size_t x = 0; x <image.size[0]; x++)
		{
			float nowInt = image.at<float>(y, x);
			if (nowInt > maxIntensity){
				maxIntensity = nowInt;
			}
		}
	}
	for (int y = 0; y < image.size[1]; y++){
		for (size_t x = 0; x < image.size[0]; x++)
		{
			image.at<float>(y, x) = image.at<float>(y, x) / maxIntensity;
		}
	}
	return image;
}
cv::Mat scaleMat(cv::Mat origenImg, int scale){
	Mat demoImg(origenImg.size()*scale, origenImg.type());
	cout << "demoImg size:" << demoImg.size[0] << "," << demoImg.size[1] << endl;
	for (int y = 0; y < origenImg.size[1]; y++){
		for (int x = 0; x < origenImg.size[0]; x++)
		{

			for (int offy = 0; offy < scale; offy++){
				for (int offx = 0; offx < scale; offx++){
					demoImg.at<float>(y*scale + offy, x*scale + offx) = origenImg.at<float>(y, x);
				}
			}
		}
	}
	return demoImg;
}
Mat IS_FeatureMatrix::showAsImage(int scale){
	//cout << "matrix size" << matrix.size[0] << "," << matrix.size[1]<<endl;
	Mat demoImg(matrix.size()*scale, matrix.type());
	//cout << "demoImg size:" << demoImg.size[0] << "," << demoImg.size[1]<<endl;
	for (int y = 0; y < matrix.size[1]; y++){
		for (int x = 0; x < matrix.size[0]; x++)
		{

			for (int offy = 0; offy < scale; offy++){
				for (int offx = 0; offx < scale; offx++){
					demoImg.at<float>(y*scale + offy, x*scale + offx) = matrix.at<float>(y, x);
				}
			}
		}
	}
	imshow("IS_FeatureMatrix", demoImg);
	waitKey(0);
	destroyAllWindows();
	return demoImg;
}
vector<cv::Mat> getIntegralImages(cv::Mat ImgN, cv::Mat ImgM){
	int N = ImgN.size[0];
	int M = ImgM.size[0];
	vector<cv::Mat> images;
	for (int i = 0; i < N; i++){
		//cout << "getIntegralImages:i=" << i<<"N="<<N<<"M="<<M;
		Mat Integ_n(M, M, CV_32FC1);
		for (int y = 0; y < M; y++){
			float s_xy = 0;
			for (int x = 0; x < M; x++){
				//cout <<"i="<<i<< " y=" << y << " x=" << x<<"   ";
				bool debug = x == 0 && y == 13;
				int ix = (i + x) %N;
				int iy = (i + y) %N;

				float An = ImgN.at<float>(iy, ix);
				float Am = ImgM.at<float>(y, x);
				s_xy += abs(An - Am);

				float ii_xy_1 = 0;
				if (y>0)
				{
					ii_xy_1 = Integ_n.at<float>(y - 1, x);
				}
				float ii_xy = ii_xy_1 + s_xy;
				Integ_n.at<float>(y, x) = ii_xy;
			}
		}
		
		images.push_back(Integ_n);
	}

	return images;
}
float calAvgDiff(int s, int t, int l, vector<cv::Mat> IntegImgs){
	//s為Am的offset
	//t為An的offset
	//cout << "calAvgDiff s:" << s << " t:" << t << "; ";
	int M = IntegImgs[0].size[0];
	int n = t- s;
	while (n < 0){//修正n,如果n是複數的話通過循環把n轉為正數
		n += IntegImgs.size();
	}
	cv::Mat IntegN = IntegImgs[n];
	//cout << "採用第" << n << "張IntegImg.";


	int b = s - 1;
	int e = b + l;//從s開始長度為l的鏈結束時的索引值會是s+l-1,舉個例子從3開始長度為5的鏈為:3,4,5,6,7 <= (7=3+5-1)
	if (e < M){
		//cout << "e<M,";
		if (b < 0){//說明n==0
			//cout << "b<0. value:" << IntegN.at<float>(e, e) << "total:" << IntegN.at<float>(e, e) << endl;
			return IntegN.at<float>(e,e)/(l*l);
		}
		else{
			//cout << " value-ee:" << IntegN.at<float>(e, e) << "value be:" << IntegN.at<float>(b, e) << "value eb:" << IntegN.at<float>(e, b) << "value bb:" << IntegN.at<float>(b, b) << "total:" << (IntegN.at<float>(e, e) - IntegN.at<float>(b, e) - IntegN.at<float>(e, b) + IntegN.at<float>(b, b)) << endl;
			return (IntegN.at<float>(e, e) - IntegN.at<float>(b, e) - IntegN.at<float>(e, b) + IntegN.at<float>(b, b))/(l*l);
		}
	}
	else
	{
		cout << "錯誤的計算! s=" << s << " t=" << t << " l=" << l << endl;
		//cout << "cycleBoundary:" << b + l - (M - 1)<<" ";
		float area1 = IntegN.at<float>(M - 1, M - 1) - IntegN.at<float>(b, M - 1) - IntegN.at<float>(M - 1, b) + IntegN.at<float>(b, b);
		//cout << "(" << M - 1<<","<<
		//cout << "area1:" << area1 << "value MM:" << IntegN.at<float>(M - 1, M - 1) << "value bM:" << IntegN.at<float>(b, M - 1) << "value Mb:" << IntegN.at<float>(M - 1, b) << " value bb:" << IntegN.at<float>(b, b);
		int cycleBoundary = b + l - (M - 1);
		float area2 = IntegN.at<float>(M - 1, cycleBoundary ) - IntegN.at<float>(b, cycleBoundary);
		//cout << "area2:" << area2 << "value Mc:" << IntegN.at<float>(M - 1, cycleBoundary ) << "value bc:" << IntegN.at<float>(b, cycleBoundary);
		float area3 = IntegN.at<float>(cycleBoundary , M - 1) - IntegN.at<float>(cycleBoundary, b);
		//cout << "area3:" << area3 << "value cM:" << IntegN.at<float>(cycleBoundary , M - 1) << "value cb:" << IntegN.at<float>(cycleBoundary , b);
		float area4 = IntegN.at<float>(cycleBoundary, cycleBoundary);
		//cout << "area4:" << area4 << " total:"<<(area1 + area2 + area3 + area4)  << endl;

		return (area1 + area2 + area3 + area4)/(l*l);
	}
}

