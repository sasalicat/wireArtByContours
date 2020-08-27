#include "stdafx.h"
#include "geodesicAsynFinder.h"
#include <iostream>
//字符串分割函数
std::vector<std::string> split(std::string str, std::string pattern)
{
	std::string::size_type pos;
	std::vector<std::string> result;
	str += pattern;//扩展字符串以方便操作
	int size = str.size();

	for (int i = 0; i<size; i++)
	{

		pos = str.find(pattern, i);

		if (pos<size&&i<size - 1)
		{
			std::string s = str.substr(i, pos - i);

			result.push_back(s);

			i = pos + pattern.size() - 1;

		}
	}
	return result;
}
bool inline isEmpty(string s){
	return s.size()==0||(s.size() == 1 && (int)s[0] == 0);
}
geodesicAsynFinder::geodesicAsynFinder()
{
}


geodesicAsynFinder::~geodesicAsynFinder()
{
}
void geodesicAsynFinder::init(int contourNum,string fname){
	boundarys.resize(contourNum);
	missionList.resize(contourNum);
	fileName = fname;
}
void geodesicAsynFinder::initBoundaryRecord(int index, vector<int> boundary){
	boundarys[index] = boundary;
}
void geodesicAsynFinder::writeBoundaryData(){
	fstream boundaryFile;
	boundaryFile.open(PATH_FILE_BOUNDARY,ios::out);
	for each (vector<int> boundary in boundarys)
	{
		for each (int bidx in boundary)
		{
			boundaryFile << bidx << " ";
		}
		boundaryFile << endl;
	}
	boundaryFile.close();
}
void geodesicAsynFinder::requestPath(int cindex, int faceIdx_s, float* pos_s, int faceIdx_e, float* pos_e){
	missionList[cindex].push_back(gdc_pair(faceIdx_s,pos_s[0],pos_s[1],pos_s[2],faceIdx_e,pos_e[0],pos_e[1],pos_e[2]));
}
bool geodesicAsynFinder::writeToMissionFile(int missionNum){
	//cout << "writeToMissionFile " << missionNum << "筆請求."<<endl;
	fstream missionFile;
	missionFile.open(PATH_FILE_MISSION, ios::out);
	int counter = 0;
	while (missionNum>0||nowPointer.cindex==missionList.size())
	{
		if (missionList.size() == nowPointer.cindex){//到達尾端
			return true;
		}
		else{
			missionFile << nowPointer.cindex<<":";
			if (missionNum < (missionList[nowPointer.cindex].size() - nowPointer.lindex)){//剩餘的missionNum不足以走到line的尾端
				for (int i = 0; i < missionNum; i++){
					gdc_pair pair = missionList[nowPointer.cindex][nowPointer.lindex + i];
					missionFile << pair.to_string();
					counter++;
				}
				missionFile << endl;
				nowPointer.lindex += missionNum;
				missionNum = 0;
			}
			else{
				for (int i = nowPointer.lindex; i < missionList[nowPointer.cindex].size(); i++)//到了尾端missionNum仍然足夠
				{
					gdc_pair pair = missionList[nowPointer.cindex][i];
					missionFile << pair.to_string();
					counter++;
				}
				missionFile << endl;
				missionNum -= missionList[nowPointer.cindex].size() - nowPointer.lindex;
				nowPointer.cindex++;//切到下一contour
				nowPointer.lindex = 0;//line索引切到0
			}
		}
	}
	missionFile.close();
	return nowPointer.cindex >= missionList.size();
}

vector<vector<OpenMesh::Vec3f>> geodesicAsynFinder::readPathResult(){
	//cout << "readPathResult開始.";
	vector<vector<OpenMesh::Vec3f>> results;
	fstream  resultFile(PATH_FILE_RESULT,ios::in);
	char buffer[READ_PATH_BUFFER_SIZE];
	int count = 0;
	while (!resultFile.eof()){
		resultFile.getline(buffer, READ_PATH_BUFFER_SIZE);
		string line(buffer,resultFile.gcount());
		count++;

		vector<OpenMesh::Vec3f> path;
		if (!isEmpty(line)){
			vector<string> tokens = split(line, ";");
			for each (string token in tokens)
			{
				if (!isEmpty(token))
				{
					vector<string> xyz = split(token, " ");
					path.push_back(OpenMesh::Vec3f(stof(xyz[0]), stof(xyz[1]), stof(xyz[2])));
				}
			}
			results.push_back(path);
		}
	}
	return results;
}
void geodesicAsynFinder::doIt(){
	ostringstream oss;
	oss << "geodesicGroupFinder " << fileName;
	system(string(oss.str()).data());
}
void geodesicAsynFinder::printMissionNum(){
	printf("printMissionNum:");
	for (int i = 0; i < missionList.size(); i++){
		printf("i=%d size:%d; ", i, missionList[i].size());
	}
}