#include "stdafx.h"
#include "intersectionOptimizer.h"


intersectionOptimizer::intersectionOptimizer()
{
}
intersectionOptimizer::intersectionOptimizer(contour last, contour next, int num){
	girdNum = num;
	int debugIndex = 0;
	const float interval = 1.0f / (float)num;
	Vec3f_angle* nextSamplePoints = next.sample(sampleRate);
	contourList[0] = last;
	contourList[1] = next;

	for (int i = 0; i < num; i++){
		//if (i ==debugIndex)
		//計算deform
		cout << "第" << i << "個grid: ";

		float now_x = loopClampTo01(i*interval - interval);
		//cout << "now_x:" << now_x;
		float end_x = loopClampTo01(i*interval);
		//cout << "end_x:" << end_x;
		float next_end_x = loopClampTo01(i*interval + interval);

		gridGroups grid(i);
		OpenMesh::Vec3f p11 = last.getPosAtAngle(end_x);            //
		OpenMesh::Vec3f p12 = last.getPosAtAngle(next_end_x);
		int closetPointIndex = 0;//找到在下一條等高線上離p11歐幾里得距離最近的點
		/*for (int j = 1; j < sampleRate; j++){
		if (abs(nextSamplePoints[j].angle - now_x)<abs(nextSamplePoints[closetPointIndex].angle - now_x)){
		closetPointIndex = j;
		}
		}*/
		for (int i = 1; i < sampleRate; i++){
			if ((nextSamplePoints[i].pos - p11).length() < (nextSamplePoints[closetPointIndex].pos - p11).length()){
				closetPointIndex = i;
			}
		}

		if (nextSamplePoints[closetPointIndex].angle < now_x||nextSamplePoints[closetPointIndex].angle>end_x){//最接近的點比now_x跟前面
			//以從clostPointIndex為基準更新 now_x ,end_x,next_end_x這三個值
			now_x = loopClampTo01(nextSamplePoints[closetPointIndex].angle-interval);
			end_x = loopClampTo01(now_x + interval);
			next_end_x = loopClampTo01(end_x + interval);
		}
		int startIndex = now_x*sampleRate;
		if (now_x > end_x){//當前angle值比結束時的大說明過程中經過隊列的結尾 nowx->陣列尾端->endx 例如now_x=0.8 ->1 ->end_x=0.3
			cout << "now_x >end_x:";
			for (int i = startIndex; i<sampleRate; i++)
			{
				floatGroup line(nextSamplePoints[i].angle);//作為以第i個點為基準計算變形量
				if (nextSamplePoints[i + 1].angle>next_end_x){//同理初始點的angle 大於next_end_x則說明跨過了終點
					for (int j = i + 1; j < sampleRate; j++){//先走到終點
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						//cout << "p11:" << p11 << " p12:" << p12 << "p21:" << p21 << "p22:" << p22;
						float weight = calDeform(p11, p12, p21, p22);//+coeff_edge*calLengrhDeform(interval,nextSamplePoints[i].angle,nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "weight:" << weight << endl;
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
					for (int j = 0; nextSamplePoints[j].angle <= next_end_x; j++){//再從起點走到next_end_x
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22); //+ coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
				}
				else{
					for (int j = i+1; nextSamplePoints[j].angle <= next_end_x; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22); //+ coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
				}
				grid.groups.push_back(line);
			}
			for (int i = 0; nextSamplePoints[i].angle <= end_x; i++){
				floatGroup line(nextSamplePoints[i].angle);
				if (nextSamplePoints[i + 1].angle>next_end_x){
					for (int j = i + 1; j < sampleRate; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22); //+ coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
					for (int j = 0; nextSamplePoints[j].angle <= next_end_x; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22); //+ coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
				}
				else{
					for (int j = i + 1; nextSamplePoints[j].angle <= next_end_x; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22); //+ coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
				}
				grid.groups.push_back(line);
			}
		}
		else
		{
			cout << "now_x<=end_x";
			for (int i = startIndex; nextSamplePoints[i].angle <= end_x; i++){
				//cout << "i=" << i << "angle=" << nextSamplePoints[i].angle << ",";
				floatGroup line(nextSamplePoints[i].angle);
				if (nextSamplePoints[i + 1].angle>next_end_x){
					for (int j = i + 1; j < sampleRate; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22); //+ coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
					for (int j = 0; nextSamplePoints[j].angle <= next_end_x; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22); //+ coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
				}
				else{
					for (int j = i + 1; nextSamplePoints[j].angle <= next_end_x; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22); //+ coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
				}
				grid.groups.push_back(line);
			}
		}/*
		cout << "gird size:" << grid.groups.size() << ";"<<endl;
		for each (floatGroup line in grid.groups)
		{
			cout << "line" << line.key <<" size:"<<line.values.size()<<":";
			for each (floatPair pair in line.values)
			{
				cout << "(" << pair.key << "," << pair.value << "),";
			}
			cout << endl;
		}*/
		deformRecords.push_back(grid);
		//初始化交點的angle
		angleMappings.push_back(floatPair(i*interval,i*interval));
		//記錄grid的幾何信息
		OpenMesh::Vec3f p21 = next.getPosAtAngle(end_x);
		OpenMesh::Vec3f p22 = next.getPosAtAngle(next_end_x);
		float weight = calDeform(p11, p12, p21, p22);//+coeff_edge*calLengrhDeform(interval,end_x,next_end_x);
		float min_w=1000000;
		float p21_min_angle;
		float p22_min_angle;

		for each (floatGroup line in grid.groups)
		{

			for each(floatPair pair in line.values){
				if (pair.value < min_w){
					min_w = pair.value;
					p21_min_angle = line.key;
					p22_min_angle = pair.key;
				}
			}
		}
		if (i == 6){
			cout << endl;
		}
		//cout << "計算 grid" << i << ":" << "p11:" << p11 << " p12:" << p12 << " p21:" << p21 << " p22:" << p22 << "weight:" << weight << ";" << endl;
		rectGrid rect(p11, p12, p21, p22, weight, min_w);
		rect.minDeformPoint[0] = p11;
		rect.minDeformPoint[1] = p12;
		rect.minDeformPoint[2] = next.getPosAtAngle(p21_min_angle);
		rect.minDeformPoint[3] = next.getPosAtAngle(p22_min_angle);
		rect.minAnglePair[0] = p21_min_angle;
		rect.minAnglePair[1] = p22_min_angle;
		gridList.push_back(rect);
	}
	//delete nextSamplePoints;
	cout << "初始化之後:" << endl;
	for (int i = 0; i < gridList.size(); i++){
		cout << "grid" << i << ":面積:" << gridList[i].area << "變形量:" << gridList[i].weight << "; ";
		if (i % 2 == 0){
			cout << endl;
		}
	}
	//delete nextSamplePoints;
}
intersectionOptimizer::intersectionOptimizer(contour last, contour next, vector<float> originAngles, bool negAngel){
	cout << "originAngles size:" << originAngles.size() << endl;
	girdNum = originAngles.size();
	int debugIndex = 0;
	const float interval = 1.0f / (float)girdNum;
	Vec3f_angle* nextSamplePoints = next.sample(sampleRate);
	contourList[0] = last;
	contourList[1] = next;
	if (last.id == 13){
		cout << "id13!" << endl;
	}
	cout << "last.id:" << last.id;
	cout << "originAngles:";

	for each (float angle in originAngles)
	{
		cout << angle << ", ";
	}
	for (int g = 0; g < girdNum; g++){
		float now_x = originAngles[g] - interval;
		float end_x = originAngles[g];
		float next_end_x = originAngles[g] + interval;
		cout << "g=" << g << " now_x:" << now_x << " end_x:" << end_x << " next_end_x:" << next_end_x << endl;
		gridGroups grid(g);
		int nowRange_i = floorf(now_x);//算圈數區間,比如當nowRange為-1時now_x處於[-1,0)
		OpenMesh::Vec3f p11 = last.getPosAtAngle(end_x);
		OpenMesh::Vec3f p12 = last.getPosAtAngle(next_end_x);
		int closetPointIndex = 0;//找到在下一條等高線上離p11歐幾里得距離最近的點
		for (int i = 1; i < sampleRate; i++){
			if ((nextSamplePoints[i].pos - p11).length() < (nextSamplePoints[closetPointIndex].pos - p11).length()){
				closetPointIndex = i;
			}
		}
		if (loopTransToLayer(nextSamplePoints[closetPointIndex].angle, nowRange_i)< now_x || loopTransToLayer(nextSamplePoints[closetPointIndex].angle,nowRange_i)>end_x){//最接近的點比now_x跟前面
			//以從clostPointIndex為基準更新 now_x ,end_x,next_end_x這三個值
			//cout << "重新定位";
			/*
			float candiate_dist[3];
			candiate_dist[0] = fabsf(loopTransToLayer(nextSamplePoints[closetPointIndex].angle, nowRange_i - 1) - now_x);
			candiate_dist[1] = fabsf(loopTransToLayer(nextSamplePoints[closetPointIndex].angle, nowRange_i) - now_x);
			candiate_dist[2] = fabsf(loopTransToLayer(nextSamplePoints[closetPointIndex].angle, nowRange_i + 1) - now_x);
			int minIdx = 0;
			for (int i = 1; i < 3; i++){
				if (candiate_dist[i] < candiate_dist[minIdx]){
					minIdx = i;
				}
			}*/
			float closetAngle= loopFindNearestCycle(end_x, nextSamplePoints[closetPointIndex].angle);
			if (abs(closetAngle - end_x) < minJumpPercentage*interval){//最近的點在可以容忍範圍內
				now_x = closetAngle - interval;
				end_x = now_x + interval;
				next_end_x = end_x + interval;
				nowRange_i = floorf(now_x);
				//cout << " now_x:" << now_x << " end_x:" << end_x << " next_end_x:" << next_end_x << endl;
			}
		}
		int startIndex = loopClampTo01(now_x)*sampleRate;

		for (int i = startIndex; loopTransToLayer(nextSamplePoints[i].angle, nowRange_i) <= end_x; i++){
			floatGroup line(loopTransToLayer(nextSamplePoints[i].angle, nowRange_i));
			int nowRange_j = nowRange_i;
			int j = i + 1;
			if (j == sampleRate){
				j = 0;
				nowRange_j++;
			}
			for (; loopTransToLayer(nextSamplePoints[j].angle,nowRange_j) <= next_end_x; j++){
				OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
				OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
				float weight = calDeform(p11, p12, p21, p22);
				floatPair pair(loopTransToLayer(nextSamplePoints[j].angle,nowRange_j), weight);
				line.values.push_back(pair);
				if (j == sampleRate - 1){
					j = 0;
					nowRange_j++;
				}
			}
			if (i == sampleRate - 1){
				i = 0;
				nowRange_i++;
			}
			grid.groups.push_back(line);
		}
		deformRecords.push_back(grid);
		//初始化交點的angle
		angleMappings.push_back(floatPair(end_x, g*interval));
		//記錄grid的幾何信息
		OpenMesh::Vec3f p21 = next.getPosAtAngle(end_x);
		OpenMesh::Vec3f p22 = next.getPosAtAngle(next_end_x);
		float weight = calDeform(p11, p12, p21, p22);
		float min_w = 1000000;
		float p21_min_angle;
		float p22_min_angle;
		int count = 0;
		for each (floatGroup line in grid.groups)
		{
			for each(floatPair pair in line.values){

				if (pair.value < min_w){
					min_w = pair.value;
					p21_min_angle = line.key;
					p22_min_angle = pair.key;
				}
			}
		}

		//cout << "計算 grid" << i << ":" << "p11:" << p11 << " p12:" << p12 << " p21:" << p21 << " p22:" << p22 << "weight:" << weight << ";" << endl;
		rectGrid rect(p11, p12, p21, p22, weight, min_w);
		rect.minDeformPoint[0] = p11;
		rect.minDeformPoint[1] = p12;
		rect.minDeformPoint[2] = next.getPosAtAngle(p21_min_angle);
		rect.minDeformPoint[3] = next.getPosAtAngle(p22_min_angle);
		rect.minAnglePair[0] = p21_min_angle;
		rect.minAnglePair[1] = p22_min_angle;
		gridList.push_back(rect);
	}

}
intersectionOptimizer::intersectionOptimizer(contour last, contour next,vector<float> originAngles){
	cout << "originAngles size:" << originAngles.size()<<endl;
	girdNum = originAngles.size();
	int debugIndex = 0;
	const float interval = 1.0f / (float)girdNum;
	Vec3f_angle* nextSamplePoints = next.sample(sampleRate);
	contourList[0] = last;
	contourList[1] = next;
	cout << "originAngles:";
	for each (float angle in originAngles)
	{
		cout << angle << ", ";
	}
	bool flag_angle1_has_cross_zero = false;//理論上angle1會跨過0度軸1次
	float last_x=-1;//記錄上一次now_x,用來判斷有沒有跨過0度軸
	for (int g = 0; g < girdNum; g++){
		//if (i ==debugIndex)
		//計算deform
		cout << "第" << g << "個grid: ";
		float now_x = loopClampTo01(originAngles[g] - interval);
		if (last_x > now_x){
			flag_angle1_has_cross_zero=true;
		}
		last_x=now_x;
		cout << "now_x:" << now_x<<"originAngles[i]:"<<originAngles[g];
		float end_x = loopClampTo01(originAngles[g]);
		cout << "end_x:" << end_x;
		float next_end_x = loopClampTo01(originAngles[g] + interval);
		cout << "next_end_x:" << next_end_x;
		gridGroups grid(g);
		
		OpenMesh::Vec3f p11 = last.getPosAtAngle(end_x);            //
		OpenMesh::Vec3f p12 = last.getPosAtAngle(next_end_x);
		int closetPointIndex = 0;//找到在下一條等高線上離p11歐幾里得距離最近的點
		/*for (int j = 1; j < sampleRate; j++){
		if (abs(nextSamplePoints[j].angle - now_x)<abs(nextSamplePoints[closetPointIndex].angle - now_x)){
		closetPointIndex = j;
		}
		}*/
		for (int i = 1; i < sampleRate; i++){
			if ((nextSamplePoints[i].pos - p11).length() < (nextSamplePoints[closetPointIndex].pos - p11).length()){
				closetPointIndex = i;
			}
		}
		/*
		if (nextSamplePoints[closetPointIndex].angle < now_x || nextSamplePoints[closetPointIndex].angle>end_x){//最接近的點比now_x跟前面
			//以從clostPointIndex為基準更新 now_x ,end_x,next_end_x這三個值
			now_x = loopClampTo01(nextSamplePoints[closetPointIndex].angle - interval);
			end_x = loopClampTo01(now_x + interval);
			next_end_x = loopClampTo01(end_x + interval);
		}*/
		int startIndex = now_x*sampleRate;

		if (now_x > end_x){//當前angle值比結束時的大說明過程中經過隊列的結尾 nowx->陣列尾端->endx 例如now_x=0.8 ->1 ->end_x=0.3
			//cout << "now_x >end_x:";
			flag_angle1_has_cross_zero = true;//第一次now_x>end_x視為跨過0度軸行為
			for (int i = startIndex; i<sampleRate; i++)
			{
				floatGroup line(nextSamplePoints[i].angle);//作為以第i個點為基準計算變形量
				line.flag = true;
				if (nextSamplePoints[i + 1].angle>next_end_x){//同理初始點的angle 大於next_end_x則說明跨過了終點
					for (int j = i + 1; j < sampleRate; j++){//先走到終點
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						//cout << "p11:" << p11 << " p12:" << p12 << "p21:" << p21 << "p22:" << p22;
						float weight = calDeform(p11, p12, p21, p22) + coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "weight:" << weight << endl;
						floatPair pair(nextSamplePoints[j].angle, weight);
						pair.flag = true;
						line.values.push_back(pair);
					}
					for (int j = 0; nextSamplePoints[j].angle <= next_end_x; j++){//再從起點走到next_end_x
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22) + coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
				}
				else{
					for (int j = i + 1; nextSamplePoints[j].angle <= next_end_x; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22) + coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
				}
				grid.groups.push_back(line);
			}
			for (int i = 0; nextSamplePoints[i].angle <= end_x; i++){
				floatGroup line(nextSamplePoints[i].angle);
				if (nextSamplePoints[i + 1].angle>next_end_x){
					for (int j = i + 1; j < sampleRate; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22) + coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						pair.flag = true;
						line.values.push_back(pair);
					}
					for (int j = 0; nextSamplePoints[j].angle <= next_end_x; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22) + coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
				}
				else{
					for (int j = i + 1; nextSamplePoints[j].angle <= next_end_x; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22) + coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
				}
				grid.groups.push_back(line);
			}
		}
		else
		{
			//cout << "now_x<=end_x"<<" now_x:"<<now_x<<"end_x:"<<end_x;

			for (int i = startIndex; nextSamplePoints[i].angle <= end_x; i++){
				//cout << "i=" << i << "angle=" << nextSamplePoints[i].angle <<" ";
				//cout << "flag:" << flag_angle1_has_cross_zero<<", ";
				floatGroup line(nextSamplePoints[i].angle);
				if (!flag_angle1_has_cross_zero){//如果還沒有跨過0度軸,即使是順位也被標記成0度之前
					line.flag = true;
				}
				if (nextSamplePoints[i + 1].angle>next_end_x){
					for (int j = i + 1; j < sampleRate; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22) + coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						pair.flag = true;
						line.values.push_back(pair);
					}
					for (int j = 0; nextSamplePoints[j].angle <= next_end_x; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						float weight = calDeform(p11, p12, p21, p22) + coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						line.values.push_back(pair);
					}
				}
				else{

					for (int j = i + 1; nextSamplePoints[j].angle <= next_end_x; j++){
						OpenMesh::Vec3f p21 = nextSamplePoints[i].pos;
						OpenMesh::Vec3f p22 = nextSamplePoints[j].pos;
						if (last.id == 9 && g == 7){
							//cout << "angle1:" << nextSamplePoints[i].angle << " angle2:" << nextSamplePoints[j].angle << " rect deform:" << calDeform(p11, p12, p21, p22) << "length deform:" << coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle)<<endl;
						}
						float weight = calDeform(p11, p12, p21, p22) + coeff_edge*calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						//cout << "calLengrhDeform(" << interval << "," << nextSamplePoints[i].angle << "," << nextSamplePoints[j].angle << "):" << calLengrhDeform(interval, nextSamplePoints[i].angle, nextSamplePoints[j].angle);
						floatPair pair(nextSamplePoints[j].angle, weight);
						pair.flag = line.flag;//順位這繼承angle1的狀態,因為順位不可能跨過0度軸
						line.values.push_back(pair);
					}
				}
				grid.groups.push_back(line);
			}
		}/*
		 cout << "gird size:" << grid.groups.size() << ";"<<endl;
		 for each (floatGroup line in grid.groups)
		 {
		 cout << "line" << line.key <<" size:"<<line.values.size()<<":";
		 for each (floatPair pair in line.values)
		 {
		 cout << "(" << pair.key << "," << pair.value << "),";
		 }
		 cout << endl;
		 }*/
		deformRecords.push_back(grid);
		//初始化交點的angle
		angleMappings.push_back(floatPair(end_x, g*interval));
		//記錄grid的幾何信息
		OpenMesh::Vec3f p21 = next.getPosAtAngle(end_x);
		OpenMesh::Vec3f p22 = next.getPosAtAngle(next_end_x);
		float weight = calDeform(p11, p12, p21, p22) + coeff_edge*calLengrhDeform(interval, end_x, next_end_x);
		float min_w = 1000000;
		float p21_min_angle;
		bool p21_cross_zero;
		float p22_min_angle;
		bool p22_cross_zero;
		int count = 0;
		for each (floatGroup line in grid.groups)
		{
			for each(floatPair pair in line.values){

				if (pair.value < min_w){
					min_w = pair.value;
					p21_min_angle = line.key;
					p21_cross_zero = line.flag;
					p22_min_angle = pair.key;
					p22_cross_zero = pair.flag;
				}
			}
		}
		//cout << "計算 grid" << i << ":" << "p11:" << p11 << " p12:" << p12 << " p21:" << p21 << " p22:" << p22 << "weight:" << weight << ";" << endl;
		rectGrid rect(p11, p12, p21, p22, weight, min_w);
		rect.minDeformPoint[0] = p11;
		rect.minDeformPoint[1] = p12;
		rect.minDeformPoint[2] = next.getPosAtAngle(p21_min_angle);
		rect.minDeformPoint[3] = next.getPosAtAngle(p22_min_angle);
		rect.minAnglePair[0] = p21_min_angle;
		rect.minAnglePair[1] = p22_min_angle;
		rect.minAngleCrossZero[0] = p21_cross_zero;
		rect.minAngleCrossZero[1] = p22_cross_zero;
		gridList.push_back(rect);
	}
	//delete nextSamplePoints;
	cout << "初始化之後:" << endl;
	for (int i = 0; i < gridList.size(); i++){
		cout << "grid" << i << ":面積:" << gridList[i].area << "變形量:" << gridList[i].weight << "; ";
		if (i % 2 == 0){
			cout << endl;
		}
	}
	//delete nextSamplePoints;
}

intersectionOptimizer::~intersectionOptimizer()
{
}

float intersectionOptimizer::calDeform(OpenMesh::Vec3f p11, OpenMesh::Vec3f p12, OpenMesh::Vec3f p21, OpenMesh::Vec3f p22){
	OpenMesh::Vec3f e1_ = p12 - p11;//contour1的邊
	OpenMesh::Vec3f e2_ = p21 - p22;//contour2的邊
	OpenMesh::Vec3f e_1 = p21 - p11;//第一組intersection point的向量
	OpenMesh::Vec3f e_2 = p22 - p12;//第二組intersection point的向量
	
	//return dot(e1_.normalize(), e2_.normalize())+dot(e_1.normalize(),e_2.normalize());
	return abs(dot(e1_.normalize(), e_1.normalize()))+abs(dot(e_1.normalize(),e2_.normalize())) +abs(dot(e2_.normalize(), e_2.normalize()))+abs(dot(e1_.normalize(),e_2.normalize()));
}
float intersectionOptimizer::calLengrhDeform(float avgLength, float angle1, float angle2){
	if (angle1 > angle2){//跨過了0度軸
		angle1 -= 1;//把angle1變成負數
	}
	return powf((angle2 - angle1) - avgLength, 2);
}
void intersectionOptimizer::solve(int iter_num){
	/*floatPair *largestWeightAnglePair=new floatPair[girdNum];
	for (int i = 0; i = deformRecords.size(); i++){
		gridGroups now_grid = deformRecords[i];
		floatPair maxAnglePair;
		float maxWeight = -1;
		for each (floatGroup group in now_grid.groups)
		{
			for each (floatPair calResult in group.values)
			{
				if (calResult.value > maxWeight){
					maxAnglePair = calResult;
				}
			}
		}
		largestWeightAnglePair[i] = maxAnglePair;
	}*/
	const float ignoreThreshold = 1.2f;
	for (int t = 0; t < iter_num; t++){
		cout <<endl<< "第" << t << "次迭代";
		int maxIndex = 0;
		for (int i = 1; i < gridList.size(); i++){
			if (gridList[i].weight > gridList[i].minWeight*ignoreThreshold){//接近最小變形量一定程度后就不再改變
				if (gridList[i].weight*gridList[i].area > gridList[maxIndex].weight*gridList[maxIndex].area){//找到最需要更新的grid
					maxIndex = i;
				}
			}
		}
		cout << "處理" << maxIndex << ":";
		rectGrid last;
		int lastIndex;
		rectGrid next;
		int nextIndex;
		if (maxIndex <= 0){
			last = gridList.back();
			lastIndex = gridList.size() - 1;
		}
		else
		{
			last = gridList[maxIndex - 1];
			lastIndex = maxIndex - 1;
		}
		if (maxIndex >= gridList.size() - 1){
			next = gridList[0];
			nextIndex = 0;
		}
		else
		{
			next = gridList[maxIndex + 1];
			nextIndex = maxIndex + 1;
		}
		//變形程度越小且面積越大的格子會分到越多的變形量
		//面積越大體現為那個格子的面積除以當前格子的面積,即相對面積
		//相對面積除以變形量,這樣變形量較小且相對面積較大的格子會得到一個較大的權重
		float weight_last = (last.area / gridList[maxIndex].area) / last.weight;
		float weight_next = (next.area / gridList[maxIndex].area)/next.weight;
		float weight_current = 1/gridList[maxIndex].weight;

		//找到所有組合中變形量最少的邊組合
		gridGroups nowGroups= deformRecords[maxIndex];
		float angle1 = -1;
		float angle2 = -1;
		float minDeform = 1000000;
		for each (floatGroup group in nowGroups.groups)
		{
			for each(floatPair pair in group.values){
				if (pair.value<minDeform)
				{
					angle1 = group.key;
					angle2 = pair.key;
					minDeform = pair.value;
				}
			}
		}
		cout << "grid" << maxIndex << "最小變形量為:" << minDeform;
		cout << "起始點1:" << gridList[maxIndex].points[0] << "起始點2:" << gridList[maxIndex].points[1]<<endl;
		cout << "最小變形量 angle1:" << angle1 << " 位置:" << contourList[1].getPosAtAngle(angle1)<< endl;
		cout << "最小變形量 angle2:" << angle2 << " 位置" << contourList[1].getPosAtAngle(angle2) << endl;
		//按權重移動角度
		float old_angle1 = angleMappings[maxIndex].value;
		float new_angle1 = (angle1 -old_angle1)*(weight_last/(weight_last+weight_current))+old_angle1;
		float old_angle2 = angleMappings[(maxIndex + 1) % angleMappings.size()].value;
		float new_angle2 = (angle2 - old_angle2)*(weight_next / (weight_next + weight_current)) + old_angle2;
		cout << "更新left  angle:" << old_angle1 << "=>" << new_angle1<<endl;
		cout << "更新right angle:" << old_angle2 << "=>" << new_angle2 << endl;
		angleMappings[maxIndex].value = new_angle1;//更新angleMapping
		angleMappings[(maxIndex + 1) % angleMappings.size()].value = new_angle2;
		OpenMesh::Vec3f pos1 = contourList[1].getPosAtAngle(new_angle1);
		OpenMesh::Vec3f pos2 = contourList[1].getPosAtAngle(new_angle2);
		//更新當前格子
		OpenMesh::Vec3f p11 = gridList[maxIndex].points[0];
		OpenMesh::Vec3f p12 = gridList[maxIndex].points[1];
		OpenMesh::Vec3f p21 = pos1;
		OpenMesh::Vec3f p22 = pos2;
		float new_w = calDeform(p11, p12, p21, p22);
		float old_w = gridList[maxIndex].weight;
		gridList[maxIndex].updatePartsPoints(NULL,NULL,&p21,&p22,new_w);
		cout << "變形量更新:" << old_w << " =>" << new_w<< endl;
		//更新前個格子
		p11 = last.points[0];
		p12 = last.points[1];
		p21 = last.points[2];
		p22 = pos1;
		new_w = calDeform(p11, p12, p21, p22);
		old_w = gridList[lastIndex].weight;
		gridList[lastIndex].updatePartsPoints(NULL, NULL, NULL, &p22, new_w);
		cout << "左變形量更新:" << old_w << "=>" << new_w << endl;
		//更新后個格子
		p11 = next.points[0];
		p12 = next.points[1];
		p21 = pos2;
		p22 = next.points[3];
		new_w = calDeform(p11, p12, p21, p22);
		old_w = gridList[nextIndex].weight;
		gridList[nextIndex].updatePartsPoints(NULL, NULL, &p21, NULL, new_w);
		cout << "右變形量更新:" << old_w << "=>" << new_w << endl;

	}
}
void bubbleSort(vector<floatPair> &marray){
	for (int end_i = marray.size() - 1; end_i >= 0; end_i--)
	{
		//cout << "end_i:" << end_i;
		for (size_t i = 0; i < end_i; i++)
		{
			//cout << "i:" << i;
			if (marray[i].value > marray[i + 1].value){
				float temp = marray[i].value;
				marray[i].value = marray[i + 1].value;
				marray[i + 1].value = temp;
			}
			//cout << "交換后marray[" << i << "]:" << marray[i].value << " marray[" << i + 1 << "]:" << marray[i + 1].value;
		}
		//cout << "排序end_i:"<<end_i<<"后的結果:";
		/*for each (floatPair pair in marray)
		{
			cout << "(" << pair.key << "," << pair.value << ")";
		}*/
		cout << endl;
		cout << endl;
	}
}
void intersectionOptimizer::solveOnce(){
	//根據energy function分佈交點
	//cout << "solve once for index:" << contourList[0].id - 1 << "angleMapping size:"<<angleMappings.size()<<endl;
	for (int i = 0; i < gridList.size(); i++){
		int lastIdx = i - 1;
		if (lastIdx < 0){
			lastIdx += gridList.size();
		}

		//變形量越大面積越小的圖源會獲得更大的權重
		//cout << "grid i=" << i  << "minAngle1:"<<gridList[i].minAnglePair[0]<<"minAngle2:"<<gridList[i].minAnglePair[1]<<";"<<endl;
		float weight_last = gridList[lastIdx].minWeight / gridList[lastIdx].area;
		float weight_current = gridList[i].minWeight / gridList[i].area;
		float t = weight_last / (weight_last + weight_current);
		float angle = 0;
		//if (i==0)
		//{
		//if (gridList[lastIdx].minAngleCrossZero[1] != gridList[i].minAngleCrossZero[0]){
		/*if (gridList[lastIdx].minAnglePair[1]>gridList[i].minAnglePair[0]&&gridList[lastIdx].minAnglePair[1]-gridList[i].minAnglePair[0]>0.5){//說明上一個點到下一個點之間跨過了0度軸
				float lastAngle = gridList[lastIdx].minAnglePair[1] - 1;
				float mixAngle = (lastAngle - gridList[i].minAnglePair[0])*t + gridList[i].minAnglePair[0];
				if (mixAngle < 0){
					angle = mixAngle + 1.0f;
				}
				else{
					angle = mixAngle;
				}
			}
		else{
				angle = (gridList[lastIdx].minAnglePair[1] - gridList[i].minAnglePair[0])*t + gridList[i].minAnglePair[0];
			}*/
		//cout << "grid"<<i<<"gridList[lastIdx].minAnglePair[1]:" << gridList[lastIdx].minAnglePair[1] << " gridList[i].minAnglePair[0]:" << gridList[i].minAnglePair[0] << " t:" << t << endl;
		float lastValue = gridList[lastIdx].minAnglePair[1];
		if (i == 0 && lastValue>gridList[i].minAnglePair[0]){//特例
			lastValue = loopFindNearestCycle(gridList[i].minAnglePair[0],lastValue);
		}
		//cout << "lastValue:" << lastValue<<endl;
		angle = (lastValue - gridList[i].minAnglePair[0])*t + gridList[i].minAnglePair[0];
		//cout << "算出的angle為:" << angle << endl;
		//}
		/*else
		{
			if ((contourList[0].id == 20)&&(i==1)){
				cout << "第19contour 第1grid cross zero1:" << gridList[lastIdx].minAngleCrossZero[1] << "grid cross zero2:" << gridList[i].minAngleCrossZero[0] << endl;
				cout << "第19contour 第1grid gridList[lastIdx].minAnglePair[1]:" << gridList[lastIdx].minAnglePair[1] << " gridList[i].minAnglePair[0]:" << gridList[i].minAnglePair[0] << " t:" << t<<endl;
			}
			angle = (gridList[lastIdx].minAnglePair[1] - gridList[i].minAnglePair[0])*t + gridList[i].minAnglePair[0];
		}*/

		angleMappings[i].value = angle;
	}
	//根據最小寬度限制調整結果

	/*for (size_t i = 0; i < angleMappings.size()-1; i++)
	{
		if (angleMappings[i].value > angleMappings[i + 1].value){
			float temp = angleMappings[i].value;
			angleMappings[i].value = angleMappings[i + 1].value;
			angleMappings[i + 1].value = temp;
		}
	}*/
	bubbleSort(angleMappings);
	if (contourList[0].id == 10){
		cout << "排序后angleMappings:";
		for (int i = 0; i < angleMappings.size(); i++){
			cout << "i=" << i << ":" << angleMappings[i].value;
		}
		cout << endl;
	}
	
	if (((int)floor(angleMappings.front().value)) == ((int)floor(angleMappings.back().value)) && angleMappings.back().value>angleMappings.front().value){

	}
	else if (loopFindNearestCycle(angleMappings.front().value,angleMappings.back().value)>angleMappings.front().value){
		float front_aft_trans = loopTransToLayer(angleMappings.front().value, floor(angleMappings.back().value));
		float back_aft_trans = loopTransToLayer(angleMappings.back().value, floor(angleMappings.front().value));
		angleMappings[0].value = back_aft_trans;
		angleMappings[angleMappings.size() - 1].value = front_aft_trans;
	}
	//cout << "交換后angleMappings:";
	//for (int i = 0; i < angleMappings.size(); i++){
	//	cout << "i=" << i << ":" << angleMappings[i].value;
	//}
	//
	const float SMALL_NUM = 0.000001f;
	int counter = 0;
	while (true)
	{
		//cout << "搜尋contour id" << contourList[0].id << "到 id" << contourList[1].id << "之間";
		float interval = 1.0f / (float)angleMappings.size();
		int minGridIndex = -1;
		for (int i = 0; i < angleMappings.size(); i++){
			float startAnglePoint = angleMappings[i].value;
			float endAnglePoint = angleMappings[(i + 1) % angleMappings.size()].value;
			//cout << "i:" << i << " startAnglePoint:" << startAnglePoint << " endAnglePoint:" << endAnglePoint << "len" << endAnglePoint - startAnglePoint << endl;
			if (i == angleMappings.size() - 1){
				//int layer = floorf(endAnglePoint);
				//startAnglePoint = loopTransToLayer(startAnglePoint,layer);
				endAnglePoint = loopFindNearestCycle(startAnglePoint, endAnglePoint);
				//cout << "startAnglePoint:" << startAnglePoint << " endAnglePoint:" << endAnglePoint << "endAnglePoint-startAnglePoint:" << endAnglePoint - startAnglePoint<<endl;
			}
			//if (contourList[0].id == 48)
				//cout << "i=" << i << "e-s=" << endAnglePoint - startAnglePoint << ",";
			if (endAnglePoint - startAnglePoint +SMALL_NUM < interval*minIntervalPercentage ){
				//cout << "i=" << i << "符合擴張條件e-s:" << endAnglePoint - startAnglePoint << "threshold:" << interval*minIntervalPercentage<<endl;
				if (minGridIndex < 0)
				{
					minGridIndex = i;
				}
				else
				{
					float angleMin = angleMappings[(minGridIndex + 1) % angleMappings.size()].value - angleMappings[minGridIndex].value;
					if ((endAnglePoint - startAnglePoint) < angleMin)
					{
						minGridIndex = i;
					}
				}
			}
		}

		if (minGridIndex<0)
		{
			break;//如果所有格子都滿足了最小間隔限制,則停止迴圈
		}
		else//擴張最小格子
		{
			if (contourList[0].id == 10){
				cout << "開始擴張" << minGridIndex << ":" << endl;
				float startAnglePoint = angleMappings[minGridIndex].value;
				float endAnglePoint = angleMappings[(minGridIndex + 1) % angleMappings.size()].value;
				cout << "length:" << endAnglePoint - startAnglePoint;
			}
			
			float expansionAngle_forward = 0;
			float expansionAngle_backward = 0;
			//int expansionCounter_forward = 0;
			//int expansionCounter_backward = 0;
			bool oddExpansion = angleMappings.size() % 2 == 0;//奇數擴張的意思是當角度從兩邊走的時候會相交于最後一格,而偶數擴張相交于最後一條邊
			float bestAngle1 = gridList[minGridIndex].minAnglePair[0];
			float bestAngle2 = gridList[minGridIndex].minAnglePair[1];
			float nowAngle1 = angleMappings[minGridIndex].value;
			bestAngle1 = loopFindNearestCycle(nowAngle1, bestAngle1);
			float nowAngle2 = angleMappings[(minGridIndex + 1) % angleMappings.size()].value;
			nowAngle2 = loopFindNearestCycle(nowAngle1,nowAngle2);
			bestAngle2 = loopFindNearestCycle(nowAngle2, bestAngle2);
			/*if (nowAngle2 - nowAngle1 < -0.5f){//判斷nowAngle1->nowAngle2之間有沒有跨0度軸
				nowAngle1 = 1 - nowAngle1;
				if (abs(bestAngle1 - nowAngle1) < 0.5)//在跨0度軸前提下,判斷bestAngle1是否在nowAngle1同一側
					bestAngle1 = 1 - bestAngle1;
				if (abs(bestAngle2 - nowAngle1) < 0.5)//判斷bestAngle2...
					bestAngle2 = 1 - bestAngle2;
			}*/
			cout << "interval*minIntervalPercentage:" << interval*minIntervalPercentage;
			if (bestAngle1 < nowAngle1 && nowAngle2 < bestAngle2){
				float offset1 = nowAngle1 - bestAngle1;
				float offset2 = bestAngle2 - nowAngle2;
				cout << "offset1:" << offset1 << "offset2:" << offset2 << "(offset2 / (offset1 + offset2)):" << (offset2 / (offset1 + offset2)) << "(offset1 / (offset1 + offset2)):" << offset1 << " abs:"<<abs(nowAngle2 - nowAngle1);
				expansionAngle_forward = (interval*minIntervalPercentage - abs(nowAngle2 - nowAngle1))*(offset2 / (offset1 + offset2));
				expansionAngle_backward = -(interval*minIntervalPercentage - abs(nowAngle2 - nowAngle1))*(offset1 / (offset1 + offset2));

			}
			else if (bestAngle1>nowAngle1 &&nowAngle2>bestAngle2)//這種情況無法處理
			{
				expansionAngle_forward = (interval*minIntervalPercentage - abs(nowAngle2 - nowAngle1)) / 2;
				expansionAngle_backward = -(interval*minIntervalPercentage - abs(nowAngle2 - nowAngle1)) / 2;
			}
			else if (bestAngle1>=nowAngle1)
			{
				expansionAngle_backward = bestAngle1 - nowAngle1;
				expansionAngle_forward = interval*minIntervalPercentage - (nowAngle2 - bestAngle1);

			} 
			else if (bestAngle2<=nowAngle2)
			{
				expansionAngle_forward = bestAngle2 - nowAngle2;
				expansionAngle_backward = - (interval*minIntervalPercentage-(bestAngle2-nowAngle1));
			}
			else
			{
				cout << "intersectionOptimizer::solveOnce() 無法處理的case nowAngle1:" << nowAngle1 << ", baseAngle1:" << bestAngle1 << ", nowAngle2:" << nowAngle2 << ", bestAngle2:" << bestAngle2 << ";" << endl;
				system("pause");
			}
			if (contourList[0].id == 10){
				cout << "nowAngle1:" << nowAngle1 << ", bestAngle1:" << bestAngle1 << ", nowAngle2:" << nowAngle2 << ", bestAngle2:" << bestAngle2 << ";" << endl;
				cout << "expansionAngle_forward:" << expansionAngle_forward << ",expansionAngle_backward:" << expansionAngle_backward << endl;
			}
			vector<float> absorb_value(angleMappings.size() - 1);//每個格子吸收的變形量比例
			float total_len = 0;//1 - interval*minIntervalPercentage;
			for (int i = 1; i < angleMappings.size(); i++){
				float angle1 = angleMappings[(minGridIndex + i) % angleMappings.size()].value;
				float angle2 = angleMappings[(minGridIndex + i + 1) % angleMappings.size()].value;
				if (minGridIndex + i == angleMappings.size() - 1){
					angle1 = loopTransToLayer(angle1, floorf(angle2));
				}
				if (angle2 - angle1 > interval*minIntervalPercentage+SMALL_NUM){//如果區域小於最小間隔的話,則這區區域是不會分攤變形量的
					total_len += angle2 - angle1;
				}
			}
			for (int i = 1; i < angleMappings.size(); i++){
				float angle1 = angleMappings[(minGridIndex + i) % angleMappings.size()].value;
				float angle2 = angleMappings[(minGridIndex + i + 1) % angleMappings.size()].value;
				if (minGridIndex + i==angleMappings.size()-1){
					angle1 = loopTransToLayer(angle1, floorf(angle2));
				}
				//cout << "i=" << i << "angle2-angle1:" << angle2 - angle1<<",";
				if (angle2 - angle1 <= interval*minIntervalPercentage+SMALL_NUM)
					absorb_value[i - 1] = 0;//吸收量為0意味著不吸收任何變形
				else
					absorb_value[i - 1] = (angle2 - angle1) / total_len;
			}
			//cout <<endl<< "absorb_value:";
			float total = 0;
			for (int i = 0; i < absorb_value.size(); i++){
				//cout << " " << i << ":" << absorb_value[i] << ",";
				total += absorb_value[i];
			}
			//cout << "總吸收量百分比之和:" << total<<endl;
			//cout << "expansionAngle_forward:" << expansionAngle_forward << "expansionAngle_backward:" << expansionAngle_backward<<endl;
			float expansionAngle_total = expansionAngle_forward - expansionAngle_backward;
			float now_expansionAngle = expansionAngle_forward;
			angleMappings[(minGridIndex + 1) % angleMappings.size()].value += now_expansionAngle;
			for (int i = 1; i < angleMappings.size(); i++){//擴張
				//cout << "i=" << i << "angle1:" << angleMappings[(minGridIndex + i) % angleMappings.size()].value << " +=" << now_expansionAngle;
				//angleMappings[(minGridIndex + i) % angleMappings.size()].value += now_expansionAngle;
				//now_expansionAngle -= expansionAngle_total*absorb_value[i - 1];//吸收所以減小了擴張角度
				//angleMappings[(minGridIndex + i+1) % angleMappings.size()].value += now_expansionAngle;
				//cout << "angle2:" << angleMappings[(minGridIndex + i + 1) % angleMappings.size()].value << "+=" << now_expansionAngle << ";    ";
				now_expansionAngle -= expansionAngle_total*absorb_value[i - 1];
				angleMappings[(minGridIndex + i + 1) % angleMappings.size()].value += now_expansionAngle;
			}
		}
		if (contourList[0].id == 13){
			cout << "id:13 第" << counter++ << "次擴張結果:";
			for (int i = 0; i < angleMappings.size(); i++){
				cout << i << ":" << angleMappings[i].value << " ";
			}
			cout << endl;
		}
	}
	cout << "solveOnce 結束.";
	if (contourList[0].id == 10){
		cout << "contour10最終結果:";
		for (int i = 0; i < angleMappings.size(); i++){
			cout << "i=" << i << ":" << angleMappings[i].value << "; ";
		}
		cout << endl;
	}
}
void intersectionOptimizer::printBestGridInf(){
	for (int i = 0; i < gridList.size(); i++){
		//cout << "第i個grid:" << gridList[i].
	}
}
vector<OpenMesh::Vec3f> intersectionOptimizer::showBestGridPoints(){
	vector<OpenMesh::Vec3f> result;
	for each (rectGrid grid in gridList)
	{
		for (size_t i = 0; i < 4; i++)
		{
			result.push_back(grid.minDeformPoint[i]);
		}
	}
	return result;
}
optimizerGroup::optimizerGroup(vector<contour*> groups,int elemNum){
	contours = groups;
	elemPerRing = elemNum;
	float interval = 1.0f / (float)elemPerRing;
	angleListForContours[0].resize(contours.size());
	angleListForContours[1].resize(contours.size());
	for (int a = 0; a < angleListForContours[0].size();a++){
		angleListForContours[0][a].resize(elemNum);
		angleListForContours[1][a].resize(elemNum);
		for (int i = 0; i < elemPerRing; i++){
			angleListForContours[0][a][i] = interval*i;
			angleListForContours[1][a][i] = interval*i;
		}
	}
	//printAngleListForContours();
}
/*void optimizerGroup::solve(){
	contourHistory.clear();
	vector<contour*> headers;
	for each (contour* c in contours)
	{
		cout << "contour" << c->id << " lastContours:";
		for each (contour* l in c->lastContours)
		{
			cout << l->id << ", ";
		}
		cout << endl;
		if (c->lastContours.size() == 0){
			headers.push_back(c);
		}
		if (c->nextContours.size()>1)
		{
			for each (contour* next in c->nextContours)
			{
				headers.push_back(next);
			}
		}
	}
	cout << "header size" << headers.size()<<":";
	for each (contour* c in headers)
	{
		cout << c->id << ",";
	}
	bool shutdown = false;
	for each(contour* current in headers)
	{
		while (current->nextContours.size()==1)
		{
			int index = current->id;
			contourHistory.push_back(index);
			cout << "current index = " << index;
			contour* next = current->nextContours[0];
			intersectionOptimizer optim(*current,*next,angleListForContours[index][0],false);
			optim.solveOnce();
			int next_index = next->id;
			cout << "nextindex = " << next_index;
			cout << "angleListForContours size:" << angleListForContours[0].size()<<"angleMappings size:"<<optim.angleMappings.size();
			for (int i = 0; i < elemPerRing; i++)
			{
				angleListForContours[0][next_index][i] = optim.angleMappings[i].value;
			}
			current = next;
		}
		if (shutdown){
			break;
		}
	}
	cout << "optimizerGroup solve over"<<endl;
	cout << "log contour index history:";
	for each (int index in contourHistory)
	{
		cout << index << ",";
	}
	cout << endl;
}*/

void optimizerGroup::stage_solve(contour* stageBegin){
	int nowIndex = stageBegin->id;
	contour* current = stageBegin;
	while (current->nextContours.size() == 1){//對於沒有分歧的contour
		//cout << "處理contour:" << nowIndex << ", ";
		contour* next = current->nextContours[0];
		int next_index = next->id;
		if (angleListForContours[0][nowIndex].size() <= OPTIMIZER_SOTP_GRID_NUM){//直接不做了
			int stageElemNum = angleListForContours[0][nowIndex].size();
			angleListForContours[1][next_index].resize(stageElemNum);
			angleListForContours[0][next_index].resize(stageElemNum);
			for (int i = 0; i < stageElemNum; i++)
			{
				angleListForContours[1][nowIndex][i] = angleListForContours[0][nowIndex][i];
				angleListForContours[0][next_index][i] = angleListForContours[1][nowIndex][i];
			}
		}
		else
		{
			intersectionOptimizer optim(*current, *next, angleListForContours[0][nowIndex], false);
			optim.solveOnce();
			int stageElemNum = angleListForContours[0][nowIndex].size();
			angleListForContours[1][nowIndex].resize(stageElemNum);
			angleListForContours[0][next_index].resize(stageElemNum);
			for (int i = 0; i < stageElemNum; i++)
			{
				angleListForContours[1][nowIndex][i] = optim.angleMappings[i].value;
				angleListForContours[0][next_index][i] = optim.angleMappings[i].value;
			}
		}
		current = next;
		nowIndex = next_index;
	}
	if (current->nextContours.size()>1){//最後停在分歧contour

		vector<vector<angleArea>> totalAreas(current->divergenceSegm.size());
		//cout << "current->divergenceSegm size():" << current->divergenceSegm.size();
		for (size_t i = 0; i < current->divergenceSegm.size(); i++)
		{
			//cout << "處理i:" << i << " contour index:" << current->nextContours[i]->id - 1 << endl;
			vector<angleArea> Areas;
			vector<fp_mapping> segm = current->divergenceSegm[i];
			for (size_t j = 0; j < segm.size(); j++)
			{
				fp_mapping mapping = segm[j];
				//cout << "處理 fp_mapping:(" << mapping[0][0] << "," << mapping[0][1] << ")>>(" << mapping[1][0] << "," << mapping[1][1] << ")";
				if (mapping[1][1] - mapping[1][0] < 0.01){

				}
				else
				{
					angleArea nowArea(mapping.group[1]);
					//cout << "angleListForContours[nowIndex].size():" << angleListForContours[nowIndex].size()<<endl;
					float boundary1 = mapping[0][0];
					float boundary2 = mapping[0][1];//mapping.group[0]的角度是直接從pfts的anglebyCal中複製的,所以只會是[0-1]
					
					bool crossZero = boundary1 > boundary2;
					//cout << "boundary1:" << boundary1 << " boundary2:" << boundary2<<"crossZero:"<<crossZero;
					vector<float> angleInside;
					for each (float angle in angleListForContours[0][nowIndex])
					{
						float angle01 = loopClampTo01(angle);//angleListForContours的值並不是位於[0-1]之內的,為了比較先拉到[0-1]之間
						if (crossZero){
							if (angle01>=boundary1 || angle01<=boundary2){
								if (angle01 > boundary2){//0點之前的點
									angle01 -=1;//將angle01折疊到[-1,1]之間,這樣方便角度之間比大小
								}
								//插入排序法,保證角度從小到大排列
								if (angleInside.size() == 0){
									angleInside.push_back(angle01);
								}
								else
								{
									
									int insertIndex;
									for (insertIndex = 0; insertIndex < angleInside.size(); insertIndex++){
										if (angleInside[insertIndex] > angle01){
											angleInside.insert(angleInside.begin() + insertIndex, angle01);
											break;
										}
									}
									if (insertIndex == angleInside.size()){//到了尾端還沒有找到插入點
										angleInside.push_back(angle01);
									}
								}
							}
						}
						else
						{
							if (angle01 >= boundary1&&angle01 <= boundary2){
								if (angleInside.size() == 0){
									angleInside.push_back(angle01);
								}
								else
								{

									int insertIndex;
									for (insertIndex = 0; insertIndex < angleInside.size(); insertIndex++){
										if (angleInside[insertIndex] > angle01){
											angleInside.insert(angleInside.begin() + insertIndex, angle01);
											break;
										}
									}
									if (insertIndex == angleInside.size()){//到了尾端還沒有找到插入點
										angleInside.push_back(angle01);
									}
								}
							}
						}
					}
					if (crossZero)
					{
						boundary1 -= 1;
					}
					if (angleInside.size() > 1)
					{
						float group0Length = boundary2 - boundary1;
						float group1Length = mapping[1][1] - mapping[1][0];//mapping.group[1]是由取樣索引點除以點數量產生的,所以位於[0,2]之間
						nowArea.fixGridAngles.resize(angleInside.size());
						//current->divergenceSegm[i][j].intersectionMapping.resize(angleInside.size());
						for (size_t p = 0; p < angleInside.size(); p++)
						{
							
							float angle_inNext = ((angleInside[p] - boundary1) / group0Length)*group1Length + mapping[1][0];
							float angle_ori = angleInside[p];
							if (angle_ori < 0){
								angle_ori += 1;
							}
							//cout<<"angle:"<<angle_ori << " angle01:" << angleInside[i]<<"angleInNext:"<<angle_inNext<<endl;
							nowArea.fixGridAngles[p] = angle_inNext;
							//current->divergenceSegm[i][j].intersectionMapping[p] = float_pair(angleInside[p], angle_inNext);
						}
						Areas.push_back(nowArea);
					}
					/*cout << "nowArea:";
					for each (float angle in nowArea.fixGridAngles)
					{
						cout << angle << ",";
					}
					cout << endl;*/
				}
			}//for (size_t j = 0; j < segm.size(); j++)
			totalAreas[i]= Areas;
		}
		for (size_t i = 0; i < totalAreas.size(); i++)
		{
			contour *next = current->nextContours[i];
			int next_index = next ->id;
			//cout << "處理分歧等高線index:" << current->id - 1 << " next_index:" << next_index << endl;
			if (totalAreas[i].size() != 0){
			   angleArea  area0= totalAreas[i][0];
				angleListForContours[0][next_index] = assignGridWithInitArea(totalAreas[i]);

			}
			else
				angleListForContours[0][next_index] = newAngleList(elemPerRing);
			stage_solve(next);
		}
	}//if (current->nextContours.size()>1)
	else//最後停在末端
	{
		//cout << "contour end" << endl;
		return;
	}
	

}
int indexOf(vector<float> list, float elem){
	for (int i = 0; i < list.size(); i++){
		if (list[i] == elem){
			return i;
		}
	}
	return -1;
}
void insertSort(vector<float> &list,float value){
	int index = 0;
	for (; index < list.size(); index++){
		if (list[index]>value){
			break;
		}
	}
	if (index==list.size())
	{
		list.push_back(value);
	}
	else
	{
		list.insert(list.begin() + index, value);
	}
}
void optimizerGroup::stage_solve(contour* stageBegin,float shift_percent){
	origenInf.resize(contours.size());
	int nowIndex = stageBegin->id;
	contour* current = stageBegin;
	int level = 0;
	vector<float> origenAngleList = angleListForContours[0][nowIndex];
	int origenIndex = nowIndex;
	while (current->nextContours.size() == 1){//對於沒有分歧的contour
		//cout << "處理contour:" << nowIndex << ", ";
		contour* next = current->nextContours[0];
		int next_index = next->id;
		if (angleListForContours[0][nowIndex].size() <= OPTIMIZER_SOTP_GRID_NUM){//直接不做了
			int stageElemNum = angleListForContours[0][nowIndex].size();//該圈的第二條等高線直接繼承第一條的交點角度
			angleListForContours[1][nowIndex].resize(stageElemNum);
			angleListForContours[0][next_index].resize(stageElemNum);
			for (int i = 0; i < stageElemNum; i++)
			{
				angleListForContours[1][nowIndex][i] = angleListForContours[0][nowIndex][i];
				float len = loopFindNearestBiggerValue(angleListForContours[0][nowIndex][i], angleListForContours[0][nowIndex][(i + 1)% angleListForContours[0][nowIndex].size()])-angleListForContours[0][nowIndex][i];
				//cout << "v2=" << loopFindNearestBiggerValue(angleListForContours[0][nowIndex][i], angleListForContours[0][nowIndex][(i + 1) % angleListForContours[0][nowIndex].size()]) << " v1=" << angleListForContours[0][nowIndex][i] << " len=" << len << "; ";
				angleListForContours[0][next_index][i] = angleListForContours[0][nowIndex][i]+len*shift_percent;//下一圈的角度增加一個shift
			}
			origenInf[nowIndex].origenContourIndex = origenIndex;
			origenInf[nowIndex].origenAngleList = origenAngleList;
			origenInf[nowIndex].subLevel = level;
			origenInf[nowIndex].totalOffset = level*shift_percent;
		}
		/*else
		{
			intersectionOptimizer optim(*current, *next, angleListForContours[0][nowIndex], false);
			optim.solveOnce();
			int stageElemNum = angleListForContours[0][nowIndex].size();
			angleListForContours[0][next_index].resize(stageElemNum);
			angleListForContours[1][nowIndex].resize(stageElemNum);
			for (int i = 0; i < stageElemNum; i++)
			{
				angleListForContours[1][nowIndex][i] = optim.angleMappings[i].value;
			}
			for (int i = 0; i < stageElemNum; i++)
			{
				float len = loopFindNearestBiggerValue(angleListForContours[1][nowIndex][i], angleListForContours[1][nowIndex][(i + 1) % angleListForContours[1][nowIndex].size()]) - angleListForContours[1][nowIndex][i];
				angleListForContours[0][next_index][i] = angleListForContours[1][nowIndex][i] + len*shift_percent;
			}
		}*/
		current = next;
		nowIndex = next_index;
		level++;
	}
	origenInf[nowIndex].origenContourIndex = origenIndex;
	origenInf[nowIndex].origenAngleList = origenAngleList;
	origenInf[nowIndex].subLevel = level;
	origenInf[nowIndex].totalOffset = level*shift_percent;
	if (current->nextContours.size()>1){//最後停在分歧contour
		vector<vector<angleArea>> totalAreas(current->divergenceSegm.size());
		//cout << "current->divergenceSegm size():" << current->divergenceSegm.size();
		//cout << "divergence cindex:" << current->id - 1;
		if (nowIndex == 29){
			cout << "index 29 contour";
		}
		for (size_t i = 0; i < current->divergenceSegm.size(); i++)
		{
			//cout << "處理i:" << i << " contour index:" << current->nextContours[i]->id - 1 << endl;
			vector<angleArea> Areas;
			vector<fp_mapping> segm = current->divergenceSegm[i];
			for (size_t j = 0; j < segm.size(); j++)
			{
				fp_mapping mapping = segm[j];
				if (nowIndex == 29)
					cout << "處理 fp_mapping:(" << mapping[0][0] << "," << mapping[0][1] << ")>>(" << mapping[1][0] << "," << mapping[1][1] << ")";
				if (mapping[1][1] - mapping[1][0] < 0.01){

				}
				else if (loopFindNearestBiggerValue(mapping[0][0],mapping[0][1]) - mapping[0][0] < 0.01){
					//cout << "mapping[0][1]:" << mapping[0][1] << " mapping[0][0]:" << mapping[0][0];
					//cout << endl;
				}
				else
				{
					angleArea nowArea(mapping.group[1]);
					//cout << "angleListForContours[nowIndex].size():" << angleListForContours[nowIndex].size()<<endl;
					
					float boundary1 = mapping[0][0];
					float boundary2 = mapping[0][1];//mapping.group[0]的角度是直接從pfts的anglebyCal中複製的,所以只會是[0-1]
					
					bool crossZero = boundary1 > boundary2;
					//cout << "boundary1:" << boundary1 << " boundary2:" << boundary2<<"crossZero:"<<crossZero;
					//for (int i = 0; i < angleListForContours[nowIndex].size(); i++){
					vector<int> gridInside;
					//for each (float angle in angleListForContours[0][nowIndex])
					/*for (int i = 0; i < angleListForContours[0][nowIndex].size();i++)
					{
						//int listSize = angleListForContours[0][nowIndex].size();
						float angle1 = angleListForContours[0][nowIndex][i];
						
						float angle01 = loopClampTo01(angle1);//angleListForContours的值並不是位於[0-1]之內的,為了比較先拉到[0-1]之間
						if (crossZero){
							if (angle01>=boundary1 || angle01<=boundary2){
								if (angle01 > boundary2){//0點之前的點
									angle01 -=1;//將angle01折疊到[-1,1]之間,這樣方便角度之間比大小
								}
								//插入排序法,保證角度從小到大排列
								if (angleInside.size() == 0){
									angleInside.push_back(angle01);
								}
								else
								{
									
									int insertIndex;
									for (insertIndex = 0; insertIndex < angleInside.size(); insertIndex++){
										if (angleInside[insertIndex] > angle01){
											angleInside.insert(angleInside.begin() + insertIndex, angle01);
											break;
										}
									}
									if (insertIndex == angleInside.size()){//到了尾端還沒有找到插入點
										angleInside.push_back(angle01);
									}
								}
							}
						}
						else
						{
							if (angle01 >= boundary1&&angle01 <= boundary2){
								if (angleInside.size() == 0){
									angleInside.push_back(angle01);
								}
								else
								{

									int insertIndex;
									for (insertIndex = 0; insertIndex < angleInside.size(); insertIndex++){
										if (angleInside[insertIndex] > angle01){
											angleInside.insert(angleInside.begin() + insertIndex, angle01);
											break;
										}
									}
									if (insertIndex == angleInside.size()){//到了尾端還沒有找到插入點
										angleInside.push_back(angle01);
									}
								}
							}
						}
					}*/
					for (int i = 0; i < origenInf[nowIndex].origenAngleList.size(); i++){
						float angle01_1 =loopClampTo01(getXat(nowIndex, i, 0));
						bool angle1_inside = false;
						if (crossZero){
							if (angle01_1 >= boundary1 || angle01_1 <= boundary2){
								angle1_inside = true;
								
							}
						}
						else
						{
							if (angle01_1 >= boundary1&&angle01_1 <= boundary2){
								angle1_inside = true;
								
							}
						}
						float angle01_2 = loopClampTo01(getXat(nowIndex, (i+1)%origenInf[nowIndex].origenAngleList.size(), 0));
						bool angle2_inside = false;
						if (crossZero){
							if (angle01_2 >= boundary1 || angle01_2 <= boundary2){
								angle2_inside = true;
							}
						}
						else
						{
							if (angle01_2 >= boundary1&&angle01_2 <= boundary2){
								angle2_inside = true;
							}
						}
						if (angle1_inside&&angle2_inside){
							gridInside.push_back(i);
						}
					}
					if (nowIndex == 29){
						cout << "gridInside:" << endl;
						for each (int gidx in gridInside)
						{
							cout << gidx << "," << endl;
						}
					}
					if (crossZero)//將原來比boundary2大的boundary1降1級方便後面和角度比大小
					{
						boundary1 -= 1;
					}
					if(true) //(angleInside.size() > 1)
					{
						vector<float> angleInside;
						current->divergenceSegm[i][j].gridMarkers.resize(gridInside.size());
						float group0Length = boundary2 - boundary1;
						float group1Length = mapping[1][1] - mapping[1][0];//mapping.group[1]是由取樣索引點除以點數量產生的,所以位於[0,2]之間,並不會有cross zero的問題
						for (int g = 0; g < gridInside.size();g++)
						{
							//cout << "gridInside:" << gridInside[g];
							current->divergenceSegm[i][j].gridMarkers[g].gridIndex = gridInside[g];
							float angle1 = getXat(nowIndex, gridInside[g], 0);
							if (crossZero&&angle1 > boundary2)
							{
								angle1 -= 1;
							}
							
							float angle1_inNext = ((loopFindNearestBiggerValue(boundary1,angle1) -boundary1) / group0Length)*group1Length + mapping[1][0];
							//cout << "angle1:" << angle1 << " inNext:" << angle1_inNext <<"loopFindNearestBiggerValue(angle1,boundary1):"<< loopFindNearestBiggerValue(angle1, boundary1);
							current->divergenceSegm[i][j].gridMarkers[g].angleCorresponds[0] = angle1_inNext;
							if (indexOf(angleInside, angle1_inNext)<0)//angleInside內沒有angle1_inNext
							{
								insertSort(angleInside, angle1_inNext);
							}
							float angle2 = getXat(nowIndex,gridInside[g],1);
							if (crossZero&&angle2 > boundary2)
							{
								angle2 -= 1;
							}
							
							float angle2_inNext = ((loopFindNearestBiggerValue(boundary1, angle2) - boundary1) / group0Length)*group1Length + mapping[1][0];
							//cout << "angle2:" << angle2<<" inNext:"<<angle2_inNext;
							current->divergenceSegm[i][j].gridMarkers[g].angleCorresponds[1] = angle2_inNext;
							if (indexOf(angleInside, angle2_inNext)<0)//angleInside內沒有angle1_inNext
							{
								insertSort(angleInside, angle2_inNext);
							}
							cout << endl;
						}
						nowArea.fixGridAngles = angleInside;
						//nowArea.fixGridAngles.resize(angleInside.size());

						//current->divergenceSegm[i][j].intersectionMapping.resize(angleInside.size());
						//current->divergenceSegm[i][j].gridMarkers = gridInside;
						/*for (size_t p = 0; p < angleInside.size(); p++)
						{
							
							float angle_inNext = ((angleInside[p] - boundary1) / group0Length)*group1Length + mapping[1][0];
							float angle_ori = angleInside[p];
							if (angle_ori < 0){
								angle_ori += 1;
							}
							//cout<<"angle:"<<angle_ori << " angle01:" << angleInside[i]<<"angleInNext:"<<angle_inNext<<endl;
							nowArea.fixGridAngles[p] = angle_inNext;
							current->divergenceSegm[i][j].intersectionMapping[p] = float_pair(angleInside[p], angle_inNext);
						}*/
						if (nowArea.fixGridAngles.size() < 2){//特殊情況,如果angleInside小於2個,說明是不正常連接的切割
							if (Areas.size()==0)//只有Area是空的的時候才會放進去當作旗標
							{
								Areas.push_back(nowArea);
							}
						}
						else{
							Areas.push_back(nowArea);
						}
					}
				}
			}//for (size_t j = 0; j < segm.size(); j++)
			totalAreas[i]= Areas;
			if (nowIndex == 29){
				cout << endl;
				for each (vector<angleArea> areas in totalAreas)
				{

					for each (angleArea area in areas)
					{
						cout << "(";
						for each (float angle in area.fixGridAngles)
						{
							cout << angle << ",";
						}
						cout << ")+";
					}
					cout << endl;
				}
				cout << "---------------------------------------------------------";
			}
			
		}

		for (size_t i = 0; i < totalAreas.size(); i++)
		{
			contour *next = current->nextContours[i];
			int next_index = next ->id;
			if (nowIndex == 29){
				cout << "處理分歧等高線index:" << current->id<< " next_index:" << next_index << endl;
			}
			if (totalAreas[i].size() != 0){//對應第i個次級等高線的區域存在
				if (totalAreas[i].size()==1&&totalAreas[i][0].fixGridAngles.size()==0){//只有一個等高線之間的區域且區域內沒有劃分表面的連線
					float length1 = current->length;//直接重新生成angleList
					float length2 = next->length;
					int pointNum = angleListForContours[0][nowIndex].size();
					int newNum = ceilf((length2 / length1)*pointNum);
					angleListForContours[0][next_index] = newAngleList(newNum);
				}
				else if (totalAreas[i].size() == 1 && totalAreas[i][0].fixGridAngles.size() == 1){//只有一個等高線之間的區域且區域內只有一條劃分表面的連線
					float length1 = current->length;//以那個劃分連線在次級等高線上的點為其中一個angle創建新的angleList
					float length2 = next->length;
					int pointNum = angleListForContours[0][nowIndex].size();
					int newNum = ceilf((length2 / length1)*pointNum);
					angleListForContours[0][next_index] = newAngleList(totalAreas[i][0].fixGridAngles[0], newNum);
				}
				else{//以多個區域的劃分連線在下級等高線內的端點為一定要存在的angle,其餘沒有端點的部分以劃分之間的最短間隔內插出新的angle
					if (totalAreas[i][0].fixGridAngles.size()<2)//???
					{
						totalAreas[i].erase(totalAreas[i].begin());
					}
					vector<float> newRing = assignGridWithInitArea(totalAreas[i]);//使用已有區域的端點內插出整條angleList
					vector<float> newRing_shift(newRing.size());
					for (size_t i = 0; i < newRing.size(); i++)
					{
						float angle1 = newRing[i];
						float angle2 = newRing[(i+1)%newRing.size()];
						if (i == newRing.size() - 1){
							angle2 = angle2 + 1;
						}
						newRing_shift[i] = angle1 + (angle2 - angle1)*shift_percent;
					}
					angleListForContours[0][next_index] = newRing_shift;
					if (nowIndex == 29){
						for each (angleArea area in totalAreas[i])
						{
							cout << "(";
							for each (float angle in area.fixGridAngles)
							{
								cout << angle << ",";
							}
							cout << ")+";
						}
						cout << endl;
						cout << "newRing:";
						for each (float angle in newRing)
						{
							cout << angle << ",";
						}
						cout << endl;
					}

					/*angleListForContours[0][next_index].resize(newRing.size());
					for (size_t j = 0; j < newRing.size(); j++)
					{
						float angle1 = newRing[j];
						float angle2 = loopFindNearestBiggerValue(newRing[j], newRing[(j + 1) % newRing.size()]);
						angleListForContours[0][next_index][j] = angle1 + (angle2 - angle1)*shift_percent;
					}*/
				}
			}
			else{
				float length1 = current->length;
				float length2 = next->length;
				int pointNum = angleListForContours[0][nowIndex].size();
				int newNum = ceilf((length2 / length1)*pointNum);
				angleListForContours[0][next_index] = newAngleList(newNum);
			}
			stage_solve(next, shift_percent);
		}
	}//if (current->nextContours.size()>1)
	else//最後停在末端
	{
		//cout << "contour end" << endl;
		return;
	}
	

}
vector<float> optimizerGroup::getAngleListForContour(int index, bool contour1){
	if (contour1)
		return angleListForContours[0][index];
	else
		return angleListForContours[1][index];
}
void optimizerGroup::printAngleListForContours(){
	for (size_t i = 0; i < angleListForContours[0].size(); i++)
	{
		vector<float> ring1 = angleListForContours[0][i];
		vector<float> ring2 = angleListForContours[1][i];

		cout << "ring1 for contour" << contours[i]->id << ":";
		for each (float angle in ring1)
		{
			cout << angle << ",";
		}
		cout << endl;
		if (contours[i]->nextContours.size() == 1){
			cout << "ring2 for contour" << contours[i]->nextContours[0] << ":";
			for each (float angle in ring2)
			{
				cout << angle << ",";
			}
			cout << endl;
		}
		else{
			cout << "ring2 沒有意義.";
		}
	}
}
void optimizerGroup::printAngleListForContours_stage(vector<contour*> contours,int beginIndex){
	int nowIndex = beginIndex;
	cout << "print stage root index:" << nowIndex << "---------------------------";
	while (contours[nowIndex]->nextContours.size()==1)
	{
		cout << "angleList " << nowIndex << "(size"<<angleListForContours[nowIndex].size()<<"):";
		
		for (size_t i = 0; i < angleListForContours[0][nowIndex].size(); i++)
		{
			cout << i << ": " << angleListForContours[0][nowIndex][i] << "=>" << angleListForContours[1][nowIndex][i] << ";";
		}
		nowIndex = contours[nowIndex]->nextContours[0]->id;
	}
	if (contours[nowIndex]->nextContours.size()>0){
		for each (contour* next in contours[nowIndex]->nextContours)
		{
			cout << ">>next index" << next->id<< ":" << endl;
			printAngleListForContours_stage(contours, next->id);

		}
	}
}
float calDeformEnergy(contour current, float at_angle,float interval){
	if (current.nextContours.size() != 1){
		return -1;
	}
	contour next = *current.nextContours[0];
	float distance = (current.getPosAtAngle(at_angle) - next.getPosAtAngle(at_angle)).length();
	return distance / interval;
}
optimizerGroup::optimizerGroup(vector<contour*> contours, int elemPerRing, float interval){
	this->contours = contours;
	this->contourInterval = interval;
}
void optimizerGroup::calAllDeformEnergy(){
	cout << "contour size:" << contours.size();
	for each (contour *c in contours)
	{
		cout << "contour" << c->id<<":"<<endl;
		if (c->nextContours.size() != 1){
			vector<deformLabel> energyList(sampleRate);
			Vec3f_angle *sampleList = c->sample(sampleRate);
			for (int i = 0; i < sampleRate; i++){
				energyList[i] = deformLabel(sampleList[i], 1);
			}
			deformEnergys.push_back(energyList);
		}
		else
		{
			vector<deformLabel> energyList(sampleRate);
			Vec3f_angle *sampleList = c->sample(sampleRate);
			for (int i = 0; i < sampleRate; i++){
				float energy= calDeformEnergy(*c, sampleList[i].angle, contourInterval);
				energyList[i]=deformLabel(sampleList[i],energy);
				//cout << i << ":" << energyList[i].energy << ",";
				if (energy == 0){
					//cout << "找到energy為0的點";
				}
			}
			//cout << "energyList size:" << energyList.size();
			deformEnergys.push_back(energyList);
			
		}
	}
}
float optimizerGroup::getXat(int contour_index, int elem_index, float x){//輸入第elem_index為指定第幾個elem,x是位於此elem中的百分比
	offsetInf inf = origenInf[contour_index];
	/*cout << "印出origenAngleList:";
	for (size_t i = 0; i < inf.origenAngleList.size(); i++)
	{
		cout << inf.origenAngleList[i] << ",";
	}
	cout << endl;*/
	float total_offset = inf.totalOffset;
	int index_offset = floorf(total_offset + x);
	elem_index += index_offset;//加上offset
	elem_index = elem_index%inf.origenAngleList.size();
	//cout << "elem_index:" << elem_index;
	x = (total_offset + x) - index_offset;

	float angle1 = inf.origenAngleList[elem_index];
	float angle2 = 0; //inf.origenAngleList[elem_index + 1];
	if (elem_index==origenInf[contour_index].origenAngleList.size()-1)//最後一個grid
	{
		angle2 = inf.origenAngleList[0] + 1;
	}
	else
	{
		angle2 = inf.origenAngleList[elem_index + 1];
	}
	//cout << "angle1:" << angle1 << " angle2:" << angle2 << "angle1 + (angle2 - angle1)*x:" << angle1 + (angle2 - angle1)*x;
	return angle1 + (angle2 - angle1)*x;
}
vector<float> assignGridWithInitArea(vector<angleArea> initArea){//用於初始話分歧等高線的初始交點分佈.已知固有區域initArea,其交點角度組為initArea[i].fixGridAngles,回傳由此內插出的等高線的所有交點的角度
	vector<float_pair> maxEmptyArea;//最大的空白區域意味著減掉了最小的有東西的區域(minAngleBoundary)
	//vector<float_pair> minEmptyArea;
	vector<angleArea> sortedInitArea;
	//cout << "initArea size:" << initArea.size()<<" ";
	while (initArea.size()>0)//插入排序發,讓
	{
		int minStartAngleIndex = 0;
		for (size_t j = 1; j < initArea.size(); j++)
		{
			if (initArea[j].minAngleBoundary()[0] < initArea[minStartAngleIndex].minAngleBoundary()[0]){
				minStartAngleIndex = j;
			}
		}
		sortedInitArea.push_back(initArea[minStartAngleIndex]);
		initArea.erase(initArea.begin()+minStartAngleIndex);
	}
	/*cout << "sort 后結果:";
	for each(angleArea area in sortedInitArea){
		cout << "[";
		for each (float angle in area.fixGridAngles)
		{
			cout << angle << ",";
		}
		cout << "];";
	}*/
	angleArea firstarea = sortedInitArea[0];
	//環切割開始
	float minGridAngle = 1000000;//計算一個最小的交點間間隔所佔的角度,作為內插出空白交點的基準
	for each (angleArea area in sortedInitArea)
	{
		for (size_t i = 0; i < area.fixGridAngles.size()-1; i++)
		{
			float interval = area.fixGridAngles[i + 1] - area.fixGridAngles[i];
			if (interval < minGridAngle){
				minGridAngle = interval;
			}
		}
	}
	//cout << "minGridAngle:" << minGridAngle;

	vector<float> finalAngleList;
	for (int i = 0; i < sortedInitArea.size(); i++){
		//cout << "i=" << i << ":";
		angleArea now_fixArea = sortedInitArea[i];
		angleArea next_fixArea = sortedInitArea[(i + 1) % sortedInitArea.size()];
		finalAngleList.insert(finalAngleList.end(), now_fixArea.fixGridAngles.begin(), now_fixArea.fixGridAngles.end());
		float next_Angle = loopFindNearestCycle(now_fixArea.minAngleBoundary()[1], next_fixArea.minAngleBoundary()[0]);//loopTransToLayer(next_fixArea.minAngleBoundary()[0], floor(now_fixArea.minAngleBoundary()[1]));
		//cout << "next_Angle:" << next_Angle << "now_fixArea.minAngleBoundary()[1]:" << now_fixArea.minAngleBoundary()[1];
		float emptyLength = next_Angle - now_fixArea.minAngleBoundary()[1];
		//cout << "emptyLength:" << emptyLength;
		int newGridNum = emptyLength / minGridAngle;//求出在emptyLength中最多能放置幾個長度為minGridAngle的格子
		//cout << " newGridNum:" << newGridNum;
		float gridLen = emptyLength / newGridNum;
		//cout << "gridLen:" << gridLen;
		for (int j = 1; j < newGridNum; j++)//0個格子或者1個格子的話就不用插入新的交點,2個格子插入一個交點把區域分為2段,以此類推
		{
			finalAngleList.push_back(now_fixArea.minAngleBoundary()[1] + gridLen*j);
		}
		//cout << endl;
	}
	return finalAngleList;
}
vector<float> newAngleList(int elemNum){
	float interval = 1.0f /(float)elemNum;
	vector<float> result(elemNum);
	for (int i = 0; i < elemNum; i++){
		result[i] = i*interval;
	}
	return result;
}
vector<float> newAngleList(float startAngle, int elemNum){
	float interval = 1.0f / (float)elemNum;
	vector<float> result(elemNum);
	for (int i = 0; i < elemNum; i++){
		result[i] = startAngle+ i*interval;
	}
	return result;
}
