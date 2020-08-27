#include "stdafx.h"
#include "polarMap.h"
#include <sstream>
#include <time.h>

const int debug_num =1;
int debug_pointIdx[debug_num] = {4809};
int show_vertexList[96]{477,462,0,10,25,40,55,70, 85,100,115,130,145,160,175,190, 206,221,263,251,266,281,297,312, 327,342,357,372,387,402,417,432
,328,343,358,373,388,403,418,433, 448,463,1,11,26,41,56,71, 87,101,116,131,146,161,176,191, 207,222,237,252,267,282,298,313
,208,223,238,253,268,283,299,314, 329,344,259,374,389,404,419,434, 449,464,2,12,27,42,57,72, 87,102,117,132,147,162,177,192};
float calEdgeTime(OpenMesh::Vec3f startPoint, OpenMesh::Vec3f edgePoint1, OpenMesh::Vec3f edgePoint2){
	float x_s = startPoint[0]; float x_1 = edgePoint1[0]; float x_2 = edgePoint2[0];
	float y_s = startPoint[1]; float y_1 = edgePoint1[1]; float y_2 = edgePoint2[1];
	float z_s = startPoint[2]; float z_1 = edgePoint1[2]; float z_2 = edgePoint2[2];
	float denominator = -(x_2*x_2 - 2 * x_2*x_1 + x_1*x_1
		+ y_2*y_2 - 2 * y_2*y_1 + y_1*y_1
		+ z_2*z_2 - 2 * z_2*z_1 + z_1*z_1);
	float numerator = x_1*x_2 - x_s*x_2 - x_1*x_1 + x_1*x_s
		+ y_1*y_2 - y_s*y_2 - y_1*y_1 + y_1*y_s
		+ z_1*z_2 - z_s*z_2 - z_1*z_1 + z_1*z_s;
	return numerator / denominator;
}

OpenMesh::Vec3f calCloestPoint(OpenMesh::Vec3f startPoint, OpenMesh::Vec3f edgePoint1, OpenMesh::Vec3f edgePoint2){
	/*float x_s = startPoint[0]; float x_1 = edgePoint1[0]; float x_2 = edgePoint2[0];
	float y_s = startPoint[1]; float y_1 = edgePoint1[1]; float y_2 = edgePoint2[1];
	float z_s = startPoint[2]; float z_1 = edgePoint1[2]; float z_2 = edgePoint2[2];
	float denominator = -(x_2*x_2 - 2 * x_2*x_1 + x_1*x_1
		+ y_2*y_2 - 2 * y_2*y_1 + y_1*y_1
		+ z_2*z_2 - 2 * z_2*z_1 + z_1*z_1);
	float numerator = x_1*x_2 - x_s*x_2 - x_1*x_1 + x_1*x_s
		+ y_1*y_2 - y_s*y_2 - y_1*y_1 + y_1*y_s
		+ z_1*z_2 - z_s*z_2 - z_1*z_1 + z_1*z_s;*/
	float t = calEdgeTime(startPoint, edgePoint1, edgePoint2);//numerator / denominator;
	return edgePoint1 + (edgePoint2 - edgePoint1)*t;
}
float angleOf(OpenMesh::Vec3f v1, OpenMesh::Vec3f v2){
	 if (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]){
		 return 0;
	 }
	float ans= OpenMesh::dot(v1, v2);
	//cout << "dot result:" << ans << ";v1 length:" << v1.length() << ";v2 length:" << v2.length()<<endl;
	ans = ans / (v1.length()*v2.length());

	return acosf(ans);
}
 anglepair cal_phki_phij(OpenMesh::Vec3f ek, OpenMesh::Vec3f ej, float Uk, float Uj){
	 //cout << "ek:" << ek << "length" << ek.length() << " ej:" << ej <<"length:"<<ej.length()<< " dot:" << OpenMesh::dot(ek, ej)<<endl;
	float Ukj = sqrtf(powf(ek.length(), 2) + powf(ej.length(), 2) - 2 * OpenMesh::dot(ek, ej));
	//cout << "Ukj 為:" << Ukj<<endl;
	float angle_ikj = acosf((powf(ej.length(), 2) - powf(ek.length(), 2) - powf(Ukj, 2)) / (-2 * ek.length()*Ukj));
	float angle_jks = acosf((powf(Uj, 2) - powf(Uk, 2) - powf(Ukj, 2)) / (-2*Ukj*Uk));
	float angle_k = angle_ikj + angle_jks;
	float Uijk = sqrtf(powf(ek.length(),2) + powf(Uk,2) - 2*ek.length()*Uk*cosf(angle_k));//distant for vi to virtual source s
	float ph_ki = acosf((powf(ek.length(), 2) - powf(Uk, 2) - powf(Uijk, 2)) / (-2*Uk*Uijk));
	float ph_kj = acosf((powf(Ukj,2)-powf(Uk,2)-powf(Uj,2))/(-2*Uk*Uj));
	float ph_ij = ph_kj - ph_ki;
	return anglepair(ph_ki,ph_ij);
}
 anglepair cal_phki_phij(OpenMesh::Vec3f ek, OpenMesh::Vec3f ej, float Uk, float Uj,bool debug){
	 if (debug){
		 //cout << "ek:" << ek << "length" << ek.length() << " ej:" << ej <<"length:"<<ej.length()<< " dot:" << OpenMesh::dot(ek, ej)<<endl;
		 cout << "cal_phki_phij debug:";
		 float Ukj = sqrtf(powf(ek.length(), 2) + powf(ej.length(), 2) - 2 * OpenMesh::dot(ek, ej));
		 cout << "Ukj:" << Ukj;
		 //cout << "Ukj 為:" << Ukj<<endl;
		 float angle_ikj = acosf((powf(ej.length(), 2) - powf(ek.length(), 2) - powf(Ukj, 2)) / (-2 * ek.length()*Ukj));
		 cout << "before acos:" << (powf(ej.length(), 2) - powf(ek.length(), 2) - powf(Ukj, 2)) / (-2 * ek.length()*Ukj);
		 cout << "angle_ikj:" << angle_ikj;
		 float angle_jks = acosf((powf(Uj, 2) - powf(Uk, 2) - powf(Ukj, 2)) / (-2 * Ukj*Uk));
		 cout << "angle_jks:" << angle_jks;
		 float angle_k = angle_ikj + angle_jks;
		 cout << "angle_k:" << angle_k;
		 float Uijk = sqrtf(powf(ek.length(), 2) + powf(Uk, 2) - 2 * ek.length()*Uk*cosf(angle_k));//distant for vi to virtual source s
		 cout << "Uijk:" << Uijk;
		 float ph_ki = acosf((powf(ek.length(), 2) - powf(Uk, 2) - powf(Uijk, 2)) / (-2 * Uk*Uijk));
		 cout << "ph_ki:" << ph_ki;
		 float ph_kj = acosf((powf(Ukj, 2) - powf(Uk, 2) - powf(Uj, 2)) / (-2 * Uk*Uj));
		 cout << "ph_kj:" << ph_kj;
		 float ph_ij = ph_kj - ph_ki;
		 cout << "ph_ij:" << ph_ij<<endl;
		 return anglepair(ph_ki, ph_ij);
	 }
	 else{
		 return cal_phki_phij(ek, ej, Uk, Uj);
	 }
 }
 anglepair cal_phki_phij(OpenMesh::Vec3f ek, OpenMesh::Vec3f ej, float Uk, float Uj, bool debug,ostringstream &ostm){
	 if (debug){
		 //cout << "ek:" << ek << "length" << ek.length() << " ej:" << ej <<"length:"<<ej.length()<< " dot:" << OpenMesh::dot(ek, ej)<<endl;
		 ostm << "cal_phki_phij debug:";
		 float Ukj = sqrtf(powf(ek.length(), 2) + powf(ej.length(), 2) - 2 * OpenMesh::dot(ek, ej));
		 ostm << "Ukj:" << Ukj;
		 //cout << "Ukj 為:" << Ukj<<endl;
		 float angle_ikj = acosf((powf(ej.length(), 2) - powf(ek.length(), 2) - powf(Ukj, 2)) / (-2 * ek.length()*Ukj));
		 ostm << "before acos:" << (powf(ej.length(), 2) - powf(ek.length(), 2) - powf(Ukj, 2)) / (-2 * ek.length()*Ukj);
		 ostm << "angle_ikj:" << angle_ikj;
		 float angle_jks = acosf((powf(Uj, 2) - powf(Uk, 2) - powf(Ukj, 2)) / (-2 * Ukj*Uk));
		 ostm << "angle_jks:" << angle_jks;
		 float angle_k = angle_ikj + angle_jks;
		 ostm << "angle_k:" << angle_k;
		 float Uijk = sqrtf(powf(ek.length(), 2) + powf(Uk, 2) - 2 * ek.length()*Uk*cosf(angle_k));//distant for vi to virtual source s
		 ostm << "Uijk:" << Uijk;
		 float ph_ki = acosf((powf(ek.length(), 2) - powf(Uk, 2) - powf(Uijk, 2)) / (-2 * Uk*Uijk));
		 ostm << "ph_ki:" << ph_ki;
		 float ph_kj = acosf((powf(Ukj, 2) - powf(Uk, 2) - powf(Uj, 2)) / (-2 * Uk*Uj));
		 ostm << "ph_kj:" << ph_kj;
		 float ph_ij = ph_kj - ph_ki;
		 ostm << "ph_ij:" << ph_ij << endl;
		 return anglepair(ph_ki, ph_ij);
	 }
	 else{
		 return cal_phki_phij(ek, ej, Uk, Uj);
	 }
 }
polarMap::polarMap(MyMesh mesh, int poleIndex)
{
	this->mesh = mesh;
	poleIdx = poleIndex;
	
}

OpenMesh::ArrayKernel polarMap::getMesh(){
	return mesh;
}
polarMap::~polarMap()
{
}
int polarMap::pointNum(){
	return mesh.n_vertices();
}
polarPoint polarMap:: getPointFrom(int idx){
	return polarPoint(idx,mesh.point(MyMesh::VertexHandle(idx)),U[idx],theta[idx]);
}
bool contain(vector<MyMesh::VertexHandle> list, MyMesh::VertexHandle elem){
	for (int i = 0; i < list.size(); i++){
		if (list[i].idx() == elem.idx()){
			return true;
		}
	}
	return false;
}
bool contain(vector<int> list, int elem){
	for (int i = 0; i < list.size(); i++){
		if (list[i] == elem){
			return true;
		}
	}
	return false;
}
bool contain(int *list, int elem,int size){
	for (int i = 0; i < size; i++){
		if (list[i] == elem){
			return true;
		}
	}
	return false;
}
void polarMap::debugDist(int ring){
	MyMesh::VertexHandle vh = OpenMesh::VertexHandle(poleIdx);
	//MyMesh::VertexFaceIter vf_it = mesh.vf_begin(vh);
	//MyMesh::FaceHandle fh = (*vf_it);
	vector<MyMesh::VertexHandle> outerPoint=vector<MyMesh::VertexHandle>();
	vector<MyMesh::VertexHandle> before = vector<MyMesh::VertexHandle>();
	outerPoint.push_back(vh);
	for (int i = 0; i < ring; i++){
		vector<MyMesh::VertexHandle> newPoints = vector<MyMesh::VertexHandle>();
		cout << "第" << i + 1 << "圈:";
		for (int p = 0; p < outerPoint.size(); p++){
			vh = outerPoint[p];
			for (MyMesh::VertexVertexCWIter vv_it = mesh.vv_cwbegin(vh); vv_it != mesh.vv_cwend(vh); vv_it++){
				if (!contain(outerPoint, *vv_it) && !contain(newPoints, *vv_it)&&!contain(before,*vv_it)){
					int idx = (*vv_it).idx();
					cout << "v" << idx << " U:" << U[idx] << " theta:" << theta[idx] << "; ";
					newPoints.push_back(*vv_it);
					before.push_back(*vv_it);
				}
			}
		}
		outerPoint = newPoints;
		cout << endl;
	}
}
struct vertexPath{
	float totalDistance = 0;
	vector<int> vidx_list;
	vertexPath(int origen_vidx, int next_vidx, float firstDist){
		vidx_list.push_back(origen_vidx);
		vidx_list.push_back(next_vidx);
		totalDistance = firstDist;
	}
	vertexPath(const vertexPath &before, int new_vidx, float new_dist){
		vidx_list.assign(before.vidx_list.begin(), before.vidx_list.end());
		vidx_list.push_back(new_vidx);
		totalDistance = before.totalDistance + new_dist;
	}
};

void polarMap::calRingVertex( vector<OpenMesh::VertexHandle> &candiates, float* &distArray, int* &lastVertexArray){
	OpenMesh::VertexHandle Vertex = candiates[0];
	//cout << "candiate size:" << candiates.size();
	candiates.erase(candiates.begin());
	//cout << "after erase:" << candiates.size();
	for (MyMesh::VertexVertexCWIter vvc_it = mesh.cvv_cwbegin(Vertex); vvc_it != mesh.vv_cwend(Vertex); vvc_it++)
	{
		//cout << "Vertex.idx():" << Vertex.idx();
		//cout << "(*vvc_it).idx():" << (*vvc_it).idx();
		
		float distance = (mesh.point(*vvc_it) - mesh.point(Vertex)).length();
		float new_dist = distArray[Vertex.idx()]+ distance;
		if (new_dist < distArray[(*vvc_it).idx()]){//更新
			distArray[(*vvc_it).idx()] = new_dist;
			lastVertexArray[(*vvc_it).idx()] = Vertex.idx();
			bool inserted= false;
			for (int i = 0; i < candiates.size(); i++){
				if (distArray[candiates[i].idx()]>new_dist){
					candiates.insert(candiates.begin() + i, *vvc_it);
					inserted = true;
					break;
				}
			}
			if (!inserted){
				candiates.push_back(*vvc_it);
			}
			//candiates.push_back(*vvc_it);
		}
	}
	//cout << "after insert:" << candiates.size();
}
void polarMap::calRingVertex(vector<OpenMesh::VertexHandle> &candiates, float* &distArray, int* &lastVertexArray, vector<int> vidx_domain){
	//cout << "calRingVertex 處理" << candiates[0].idx();
	OpenMesh::VertexHandle Vertex = candiates[0];
	//cout << "candiate size:" << candiates.size();
	candiates.erase(candiates.begin());
	//cout << "after erase:" << candiates.size();
	//cout << "domain size:" << vidx_domain.size();
	for (MyMesh::VertexVertexCWIter vvc_it = mesh.cvv_cwbegin(Vertex); vvc_it != mesh.vv_cwend(Vertex); vvc_it++)
	{
		//cout << "搜尋" << (*vvc_it).idx();
		if (contain(vidx_domain, (*vvc_it).idx())){
			//cout << "搜尋到!  ";
			float distance = (mesh.point(*vvc_it) - mesh.point(Vertex)).length();
			float new_dist = distArray[Vertex.idx()] + distance;
			if (new_dist < distArray[(*vvc_it).idx()]){//更新
				distArray[(*vvc_it).idx()] = new_dist;
				lastVertexArray[(*vvc_it).idx()] = Vertex.idx();
				bool inserted = false;
				for (int i = 0; i < candiates.size(); i++){
					if (distArray[candiates[i].idx()]>new_dist){
						candiates.insert(candiates.begin() + i, *vvc_it);
						inserted = true;
						break;
					}
				}
				if (!inserted){
					candiates.push_back(*vvc_it);
				}
				//candiates.push_back(*vvc_it);
			}
		}
	}
}

float calAvgPointDiff(const vector<OpenMesh::Vec3f> &ringM,int idxBeg,int idxEnd,const vector<OpenMesh::Vec3f> &ringN, int offset,float maxThreshold){
	float energy=0;
	//cout << "maxThreshold = " << maxThreshold << endl;
	for (int i = idxBeg, j=0; i <= idxEnd; i++,j++)//i遍歷ringM,j遍歷ringN
	{
		int fixi = i;
		if (i < 0){
			fixi = i + ringM.size();
		}
		float dist = (ringM[fixi] - ringN[(j + offset) % ringN.size()]).length();
		//cout << "i=" << i << ",j=" << j << "Pi:(" << ringM[fixi][0] << "," << ringM[fixi][1] << "," << ringM[fixi][2] << ") Pj:(" << ringN[(j + offset) % ringN.size()]<<")"<< " dist:" << dist;
		if (dist > maxThreshold){
			energy += (dist - maxThreshold);
		}
	}
	return energy / ringM.size();
}
void polarMap::assignDivergence(){
	clock_t time_start, time_end;
	time_start = clock();
	const int SAMPLE_POINT_NUM = 200;
	for (int i = 0; i < contours->size(); i++){
		if (contours->operator[](i)->nextContours.size()>1)//分歧的contour
		{
			//cout << "處理contour index:" << i << endl;
			map<int,vector<int_pair>> asignGroups;//索引值是對應的下一條等高線的id,值是一個陣列記錄所有劃分給這個等高線的分斷的起點pfts的索引值和終點pfts的索引值
			contour *current = contours->operator[](i);
			int currentBelongIdx =-1;
			int_pair record;
			int firstIdx = current->pointsBelong[0];
			for (int v = 0; v < current->pfts.size(); v++){
				if (current->pointsBelong[v] != currentBelongIdx){
					record = int_pair(v);
					currentBelongIdx = current->pointsBelong[v];
					//cout << "v=" << v << "創建新record;   ";
				}
				if (v == current->pfts.size()-1){
					if (currentBelongIdx == firstIdx){//說明首尾是聯通的
						asignGroups[firstIdx][0].values[0] = record.values[0]-current->pfts.size();
						//cout << "v=" << v << "設置初始值;   ";
					}
				}
				else{
					if (current->pointsBelong[v + 1] != currentBelongIdx){
						record.values[1] = v;
						asignGroups[currentBelongIdx].push_back(record);
						//cout << "v=" << v << "記錄于"<< currentBelongIdx << "[" << record.values[0] << "," << record.values[1] << "];   ";
					}
				}
			}
			/*cout << "印出asignGroups(size:" << current->pfts.size()<<"):"<<endl;
			for each (pair<int,vector<int_pair>> pair in asignGroups)
			{
				cout << "idx" << pair.first << ":";
				for each (int_pair value in pair.second)
				{
					cout << "[" << value.values[0] << "," << value.values[1] << "]";
				}
				cout << endl;
			}
			cout << "contour" << i << " assign divergence finsh" << endl;*/
			/*vector<vector<OpenMesh::Vec3f>> resampleUnion= combineSample_withoutAngle(current->nextContours, SAMPLE_POINT_NUM);
			cout << "印出resampleUnion數量:";
			for (int i = 0; i < resampleUnion.size(); i++){
				cout << " " << i << ":" << resampleUnion[i].size();
			}*/
			current->divergenceSegm.resize(current->nextContours.size());
			for each (pair<int, vector<int_pair>> group in asignGroups)
			{
				vector<OpenMesh::Vec3f> resample1 = current->sample_withoutAngle(SAMPLE_POINT_NUM);
				int contour_inNext_index = 0;
				vector<float_pair> partSegments;//(group.second.size());
				//cout << "處理group:" << group.first;
				for (size_t i = 0; i < current->nextContours.size(); i++)
				{
					//cout <<i<< ">id-1:" << current->nextContours[i]->id - 1<<" ";
					if (group.first == current->nextContours[i]->id){
						//cout << "get!";
						contour_inNext_index = i;
						break;
					}
				}
				//cout << "對應下一條contour:" << group.first;
				for (int p = 0; p < group.second.size();p++)
				//for each(int_pair pair in group.second)
				{
					clock_t start_t, end_t;
					start_t = clock();
					int_pair pair = group.second[p];
					float_pair segmemt;
					if(true){//if (i == 23 && pair.values[0] == -1){
						//cout << "subset:[" << pair.values[0] << "," << pair.values[1] << "]印出對應的最近點:";
						int idx1 = pair.values[0];
						if (idx1 < 0){
							idx1 += current->pfts.size();
						}
						float angleBeg = current->pfts[idx1].angleByCal;
						float angleEnd = current->pfts[pair.values[1]].angleByCal;
						if (loopFindNearestBiggerValue(angleBeg,angleEnd)-angleBeg<AssignIgnorePartAngleThreshold){
							break;//如果起始段太過短則結果不具有參考意義
						}
						int indexBeg = angleBeg*resample1.size();
						int indexEnd = angleEnd*resample1.size();
						if (indexBeg > indexEnd){
							indexBeg -= resample1.size();
						}
						IS_FeatureMatrix IS_M(resample1, indexBeg, indexEnd);
						//IS_FeatureMatrix IS_M(current->toVec3fList(), pair.values[0], pair.values[1]);
						int M = IS_M.matrix.size[0];

						//cout << "contour_inNext_index:" << contour_inNext_index << "contour index:" << group.first;
						float mindiff = 1000000;
						float inteval = (assign_maxPercent - assign_minPercent) / assign_segmentNum; 
						int snum = assign_segmentNum;
						//cout << "interval:" << inteval << endl;
						for (int s = 0; s <= snum; s++){
							//cout << "s:" << s;
							//IS_FeatureMatrix IS_N(resampleUnion[contour_inNext_index]);
							//IS_FeatureMatrix IS_N(contours->operator[](group.first - 1)->toVec3fList());
							float proportion = assign_minPercent + s*inteval;
							//cout << "比例:" << proportion << "; ";
							int N = M / proportion;
							vector<OpenMesh::Vec3f> resample2 = current->nextContours[contour_inNext_index]->sample_withoutAngle(N);
							
							IS_FeatureMatrix IS_N(resample2);

							//int N = IS_N.matrix.size[0];

							if (M > N){
								cout << "wdnmd M>N! where M:" << M << "N:" << N << endl;
							}
							else{
								//cout << "處理:[" << pair.values[0] << "," << pair.values[1] << "]" << endl;
								vector<cv::Mat> Integs = getIntegralImages(IS_N.matrix, IS_M.matrix);

								int possibleT = N;
								int mint = -1;
								float minDist = 1000000;
								for (int t = 0; t < possibleT; t++){
									float distEnergy = calAvgPointDiff(resample1, indexBeg, indexEnd, resample2, t, intervalBtwLevels);
									if (distEnergy < minDist)
										minDist = distEnergy;
									float diff = calAvgDiff(0, t, M, Integs) + assigh_distEnergy_weight*distEnergy;
									//if (pair.values[0]==120)
									//cout << "t:" << t << " diff:" << diff << " dist:"<<distEnergy<<"  ";
									if (diff < mindiff){
										mindiff = diff;
										mint = t;
										//cout << "mint=" << mint<<"diff:"<<diff << " dist:" << distEnergy;
										//cout << endl;
									}
								}
								//cout << "minDist:" << minDist<<endl;
								//cout << "mint:" << mint << "mindiff:" << mindiff << endl;

								if (mint >= 0){//在這次嘗試中找到了一個比之間相差量更小的值,其第一個點為mint

									for (int i = pair.values[0]; i <= pair.values[1]; i++)
									{
										int fixi = i;
										if (i<0)
										{
											fixi += current->pfts.size();
										}
										
										float angle = current->pfts[fixi].angleByCal;
										int offset = (angle - angleBeg)*resample1.size();//適用於未跨過0度軸的case

										if (angleBeg>angleEnd){//跨過0度軸
											if (angle >= angleBeg){
												//當前點從angleBeg出發還未跨0度軸
											}
											else{//跨過0度軸了
												offset = (1 - angleBeg + angle)*resample1.size();
											}
										}

										current->closetNextPoint[fixi] = resample2[(mint + offset) % resample2.size()];
	
									}
									segmemt.values[0] = (float)mint / (float)N;
									segmemt.values[1] = (float)(mint + M - 1) / (float)N;


								}

							}
						}


					}
					end_t = clock();
					double time = ((double)(end_t - start_t) / CLOCKS_PER_SEC);
					cout << "asignGroups耗時:" << time << "; ";
					partSegments.push_back(segmemt);//[p] = segmemt;
				}//for (int p = 0; p < group.second.size();p++) end here
				/*if (partSegments.size() > 1){
					cout << "調整之前parSegments:";
					for each (float_pair pair in partSegments)
					{
						cout << "(" << pair[0] << "," << pair[1] << "); ";
					}
					cout << endl;
				}
				cout << "group:" << group.first << "part segment:";*/
				//糾正部分分配的區域重複的問題
				for (int n = 0; n < partSegments.size();n++)//對所有
				{
					for (int m = 0; m < partSegments.size(); m++){
						if (n != m){
							//cout << "n=" << n << " m=" << m << ":" ;
							if ( partSegments[n][1]<1){//未跨過0度軸
								if (partSegments[m][0] < partSegments[n][1] && partSegments[m][0]>partSegments[n][0]){//第m個分段的第一端點位於n分段內
									float avgAngle = (partSegments[n][1] + partSegments[m][0]) / 2;
									partSegments[m].values[0] = avgAngle;
									partSegments[n].values[1] = avgAngle;
									//cout << "case1:" << partSegments[m].values[0]<<endl;
								}
							}
							else
							{
								float fix_n1 = partSegments[n][1] - 1;
								//cout << "fix_n1:" << fix_n1;
								if (partSegments[m][0]<fix_n1 || partSegments[m][0]>partSegments[n][0]){
									float avgAngle = (fix_n1 + partSegments[m][0]) / 2;
									partSegments[m].values[0] = avgAngle;
									partSegments[n].values[1] = avgAngle+1;
									//cout << "case2:" << partSegments[m].values[0]<<endl;
								}
							}
						}
					}
				}
				/*if (partSegments.size() > 1){
					cout << "之後parSegments:";
					for each (float_pair pair in partSegments)
					{
						cout << "(" << pair[0] << "," << pair[1] << "); ";
					}
					cout << endl;
				}*/
				current->divergenceSegm[contour_inNext_index].resize(partSegments.size());
				for (int n = 0; n < partSegments.size(); n++){
					int_pair idx_pair1 = group.second[n];
					int idx1 = idx_pair1[0];
					if (idx1 < 0){
						idx1 += current->pfts.size();
					}
					int idx2 = idx_pair1[1];
					float_pair angle_pair1(current->pfts[idx1].angleByCal, current->pfts[idx2].angleByCal);
					float_pair angle_pair2 = partSegments[n];
					current->divergenceSegm[contour_inNext_index][n] = fp_mapping(angle_pair1, angle_pair2);
				}

				//cout << endl;
			}//for each (pair<int, vector<int_pair>> group in asignGroups) end here
			/*
			for each (vector<fp_mapping> list in current->divergenceSegm)
			{
				for each (fp_mapping mapping in list)
				{
					cout << "(" << mapping[0][0] << "," << mapping[0][1] << ")>>(" << mapping[1][0] << "," << mapping[1][1] << "); ";
				}
				cout << endl;
			}
			cout << endl;*/
		}//if (contours->operator[](i)->nextContours.size()>1)

	}
	time_end = clock();
	double cost = ((double)(time_end - time_start)) / CLOCKS_PER_SEC;
	cout << "assignDivergence結束,耗時:"<<cost <<"秒"<< endl;
}
struct  candiateAxis
{
	float weight = 0;
	vector<int> vidx;
	candiateAxis(vector<int> vlist, float len){
		vidx = vlist;
		weight = len;
	}
};
struct floatInt_pair{
	float weight;
	int vidx;
	floatInt_pair(float w, int v){
		weight = w;
		vidx = v;
	}
};
vector<vector<int>> polarMap::calBaseAxis(){
	/*
	vector<vertexPath> paths4CalAxis;
	OpenMesh::VertexHandle origenVertex(poleIdx);
	for (MyMesh::VertexVertexCWIter vvc_it = mesh.cvv_cwbegin(origenVertex); vvc_it != mesh.vv_cwend(origenVertex); vvc_it++){
		int vidx = (*vvc_it).idx();
		float dist= (mesh.point(*vvc_it) - mesh.point(origenVertex)).length;
		paths4CalAxis.push_back(vertexPath(poleIdx, vidx, dist));
	}*/
	distFormVertex = new float[mesh.n_vertices()]{BIG_NUMBER};//初始化距離矩陣
	//cout << "初始化calBaseAxis:";
	for (int i = 0; i < mesh.n_vertices(); i++){

		distFormVertex[i] = BIG_NUMBER;
		
	}
	distFormVertex[poleIdx] = 0;
	int* lastVertexIndex = new int[mesh.n_vertices()];//初始化上一個點索引
	for (int idx = 0; idx < mesh.n_vertices(); idx++){
		lastVertexIndex[idx] = idx;
	}
	vector<OpenMesh::VertexHandle> candiates;
	candiates.push_back(OpenMesh::VertexHandle(poleIdx));
	while (candiates.size() > 0){//計算所有點到原點(poleIdx)的距離和路徑,這裡不是測地距離而是通過頂點之間的直接連接來就算
		calRingVertex(candiates, distFormVertex, lastVertexIndex);
	}
	for (int i = 0; i < mesh.n_vertices(); i++){
		if (distFormVertex[i] > farestFromVertex){
			farestFromVertex = distFormVertex[i];
		}
	}
	vector<vector<int>> axises;
	vector<candiateAxis> candiate_axis;

	/*cout << "localExtremum_idx:";
	for (int i = 0; i < localExtremum_idx.size(); i++){
		cout << localExtremum_idx[i] << ",";
	}*/

	for each (int traget_vidx in localExtremum_idx)//對於每一個局部極值點
	{//組裝一條從其實點到局部極值點的路徑
		//cout << "局部極值點來自:";
		vector<int> now_axis;
		now_axis.push_back(traget_vidx);
		int next_vidx = traget_vidx;
		while (next_vidx != poleIdx)
		{
			//cout << next_vidx << "<=";
			next_vidx = lastVertexIndex[next_vidx];
			//now_axis.push_back(next_vidx);
			now_axis.insert(now_axis.begin(),next_vidx);
		}
		now_axis.insert(now_axis.begin(),poleIdx);
		float total_length=0;
		for (int i = 0; i < now_axis.size()-1; i++){
			float len = (mesh.point(OpenMesh::VertexHandle(now_axis[i + 1])) - mesh.point(OpenMesh::VertexHandle(now_axis[i]))).length();
			total_length += len;
		}
		//cout << endl;
		candiate_axis.push_back(candiateAxis(now_axis,total_length));
		axises.push_back(now_axis);
	}
	
	for each(contour *c in *contours){//每一條等高線
		bool debug = false;
		vector<floatInt_pair> candiate_startPoint;
		if (debug)
			cout << "contour" << c->id <<"長度為:"<<c->pfts.size()<< " :";
		for each(candiateAxis axis in candiate_axis){

			for (int i = 0; i < axis.vidx.size()-1; i++){
				int index = c->indexOfTwoVidx(axis.vidx[i], axis.vidx[i + 1]);
				if (index>=0)
				{
					if (debug)
						cout << "找到一個可能起始點index為:" << index;
					candiate_startPoint.push_back(floatInt_pair(axis.weight,index));
				}
			}
		}
		if (debug)
			cout << endl;
		floatInt_pair max_pair = candiate_startPoint[0];
		for (int i = 1; i < candiate_startPoint.size();i++)
		{
			if (candiate_startPoint[i].weight > max_pair.weight){
				max_pair = candiate_startPoint[i];
			}
		}

		if (debug)
			cout << "contour 長度" << c->length << endl;
		int startIdx = max_pair.vidx;
		float now_angle = 0;
		c->pfts[startIdx].angleByCal = now_angle;
		//計算contour每一個點的角度
		for (int i = startIdx; i < c->pfts.size()-1; i++){
			now_angle += ((c->pfts[i + 1].pointPos - c->pfts[i].pointPos).length() / c->length);
			if (debug)
				cout << i+1 << ":" <<"nowangle:"<<now_angle<<"| ";
			c->pfts[i + 1].angleByCal = now_angle;
		}
		if (debug)
			cout << endl << "到達contour尾端" << endl;
		now_angle +=((c->pfts[c->pfts.size()-1].pointPos - c->pfts[0].pointPos).length() / c->length);
		if (debug)
			cout << 0 << ":"  << "nowangle:" << now_angle << "| ";
		c->pfts[0].angleByCal = now_angle;
		for (int i = 0; i < startIdx; i++){
			now_angle += ((c->pfts[i + 1].pointPos - c->pfts[i].pointPos).length() / c->length);
			if (debug)
				cout << i + 1 << ":" << "dist:" << (c->pfts[i + 1].pointPos - c->pfts[i].pointPos).length() << "_nowangle:" << now_angle << "| ";
			c->pfts[i + 1].angleByCal = now_angle;
		}

		c->pfts[startIdx].angleByCal = 1;
		if (c->id == 2){
			cout << "contour index:" << c->id << "startIdx:" << startIdx;
			printf("pfts[startIdx]:%e/n", c->pfts[startIdx]);
			system("pause");
		}
	}
	cout << "calBaseAxis 結束"<<endl;
	return axises;

}
float polarMap::getPercentageFromVertex(int vidx){
	return distFormVertex[vidx] / farestFromVertex;
}
void polarMap::initStartRing(){
	U[poleIdx] = 0;
	MyMesh::VertexHandle startPoint = MyMesh::VertexHandle(poleIdx);
	float angleBefore=0;
	OpenMesh::Vec3f vecBefore;
	for (MyMesh::VertexVertexCWIter vvc_it = mesh.cvv_cwbegin(startPoint); vvc_it != mesh.vv_cwend(startPoint); vvc_it++){
		MyMesh::VertexHandle vertex = (*vvc_it);
		U[vertex.idx()] = (mesh.point(vertex) - mesh.point(startPoint)).length();//設置距離
		OpenMesh::Vec3f outer = mesh.point(vertex) - mesh.point(startPoint);//計算本邊的向量

		if (!inited){
			baseDirection = outer;
			vecBefore = baseDirection;
			inited = true;
			zero_anix_vertexIdx= vertex.idx();
		}
		//cout << "vec before:" << vecBefore;
		//cout << " outer:" << outer << endl;
		angleBefore += angleOf(vecBefore, outer);//將自己的角度累加上去
		theta[vertex.idx()] = angleBefore;//設置角度
		vecBefore = outer;//設置的邊給下一個頂點算角度
		//cout << "初始化vertex:" << vertex.idx() << " distance:" << U[vertex.idx()] << " angle:" << theta[vertex.idx()]<<endl;
	}


}
void polarMap::findExtremumVertexs(){

	localExtremum_idx.clear();
	for (int i = 0;i<mesh.n_vertices(); i++){
		//cout << "i=" << i<<":";
		MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
		bool isExtremum = true;//極值點的旗標
		for (MyMesh::VertexVertexCWIter vv_cwit = mesh.vv_cwbegin(vh); vv_cwit != mesh.vv_cwend(vh); vv_cwit++){
			int vidx = (*vv_cwit).idx();
			if (U[vidx] > U[i]){//如果周圍的任意一點比當前點離原點更遠,則當前點不是起始點
				//cout << " > ";
				isExtremum = false;
				break;
			}
		}
		//cout << ";"<< endl;
		if (isExtremum){//找到極值點
			localExtremum_idx.push_back(i);
			//cout << "U[" << i << "]:" << U[i] << "=>";
			for (MyMesh::VertexVertexCWIter vv_cwit = mesh.vv_cwbegin(vh); vv_cwit != mesh.vv_cwend(vh); vv_cwit++){
				int vidx = (*vv_cwit).idx();
			//	cout << "idx " << vidx << ":" << U[vidx];
			}
			//cout << endl;
		}
	}
	cout << "找到的所有極值點:";
	for each (int idx in localExtremum_idx)
	{
		cout << idx << ",";
	}
	cout << endl;

}
template <typename T>
bool remove(vector<T> &list, T elem){
	for (size_t i = 0; i < list.size(); i++)
	{
		if (list[i] == elem)
		{
			list.erase(list.begin() + i);
			return true;
		}
	}
}
vector<contour*> polarMap::getPossibleBoundary(float U){
	vector<contour*> result;
	if (U < contours->operator[](0)->distance)//在第一圈等高線之前
	{
		result.push_back(contours->operator[](0));
		int count = 1;
		while (contours->operator[](count)->distance == result[0]->distance)
		{
			result.push_back(contours->operator[](count++));
		}
		return result;
	}
	else if (U>contours->back()->distance){//在最後一圈等高線之後
		result.push_back(contours->back());
		int count = contours->size() - 2;
		while (contours->operator[](count)->distance == result[0]->distance){
			result.push_back(contours->operator[](count--));
			cout << "count=" << count << "result[0]->distance:" << result[0]->distance;
		}
		return result;
	}
	for (size_t i = 0; i < contours->size()-1; i++)
	{
		if (contours->operator[](i)->distance <= U &&contours->operator[](i+1)->distance >= U)
		{
			//把所有同一距離的前邊界加入result
			float lowerDist = contours->operator[](i)->distance;
			result.push_back(contours->operator[](i));
			int count = i-1;
			while (count >= 0 && contours->operator[](count)->distance==lowerDist)
			{
				result.insert(result.begin(), contours->operator[](count));
				count--;
			}
			//把所有同一距離的后邊界加入result
			float upperDist = contours->operator[](i + 1)->distance;
			result.push_back(contours->operator[](i + 1));
			count = i + 2;
			while (count < contours->size() && contours->operator[](count)->distance == upperDist){
				result.push_back(contours->operator[](count));
				count++;
			}
		}
	}
	return result;
}
vector<contour*> polarMap::getPossibleBoundary(float U, int extraLayer){
	vector<contour*> results;
	vector<float> tragetDist;
	//int firstUpBondary_i = 0;;

	//cout << "U:" << U << "intervalBtwLevels:" << intervalBtwLevels<<endl;
	/*for (int i = 1; i <= totalLevel; i++){
		if (U < i*intervalBtwLevels){
			firstUpBondary_i = i;
			break;
		}
	}*/
	int firstUpBondary_i = U / intervalBtwLevels + 1;
	//cout << "firstUpBondary_i:" << firstUpBondary_i;
	for (size_t i = firstUpBondary_i; i <= firstUpBondary_i + extraLayer; i++)
	{
		if (i <= totalLevel&&i >= 1){
			tragetDist.push_back(i*intervalBtwLevels);
			//cout << "push back i" << i << "i*intervalBtwLevels:" << i*intervalBtwLevels;
		}
	}
	for (size_t i = firstUpBondary_i - 1; i >= firstUpBondary_i-1 - extraLayer; i--){
		if (i <= totalLevel&&i >= 1){
			tragetDist.push_back(i*intervalBtwLevels);
			//cout << "push back i" << i << "i*intervalBtwLevels:" << i*intervalBtwLevels;
		}
	}
	/*cout << "tragetDist:";
	for each (float dist in tragetDist)
	{
		cout << dist << ",";
	}*/
	for each (contour* c in *contours)
	{
		//cout << "c id:" << c->id << "dist:" << c->distance;
		for each (float v in tragetDist)
		{
			if (c->distance == v){
				results.push_back(c);
			}
		}
	}
	return results;
}
void polarMap::exploreVertex(int nowVIdx, vector<int> &regionVertexIdx, set<int> &boundaryContourIdx,bool debug){
	//debug = false;
	if (debug){
		cout << "處理vertex" << nowVIdx<<":距離"<<U[nowVIdx]<<" ";
	}
	OpenMesh::VertexHandle now_vh(nowVIdx);
	regionVertexIdx.push_back(nowVIdx);
	float vdist = U[nowVIdx];
	if (debug)
		cout << "dist:" << vdist<<" possibleBoundary:";
	//ector<contour*> possibleBoundary = getPossibleBoundary(vdist);
	if (debug){
		for each (contour* boundary in explorePossibleBoundary)
		{
			cout << boundary->id << ", ";

		}
	}
	for (MyMesh::VertexVertexCWIter vv_cwit = mesh.vv_cwbegin(now_vh); vv_cwit != mesh.vv_cwend(now_vh); vv_cwit++){
		if (contain(regionVertexIdx, (*vv_cwit).idx()))//如果是之前的點則啥都不做
		{
			continue;
		}
		if (debug){
			cout << "周圍點" << (*vv_cwit).idx() << " 可能邊界" << explorePossibleBoundary.size() << "條";
		}
		bool hasMeetBoundary = false;
		/*
		for each(contour* boundary in possibleBoundary){
			if (debug)
				boundary->debug_contact((*vv_cwit).idx());
			if (boundary->contact((*vv_cwit).idx())){
				if (debug)
					cout << " 周圍點:" << (*vv_cwit).idx() << "meetBoundary" << boundary->id<<"!";
				boundaryContourIdx.insert(boundary->id);
				meetBoundary = true;
				//break;
			}
		}*/
		/*
		if (explorePossibleBoundary.size() == 2){//兩條等高線
			meetBoundary = explorePossibleBoundary[0]->contact_asSmall((*vv_cwit).idx()) || possibleBoundary[1]->contact_asLarge((*vv_cwit).idx());//作為小的等高線的小源頭或者作為大的等高線的大源頭,說明檢索到了超出本區域的點
		}
		else if (explorePossibleBoundary.size() == 1){//一條等高線,說明是極值存在區域
			meetBoundary = explorePossibleBoundary[0]->contact_asSmall((*vv_cwit).idx());//作為唯一邊界等高線的小源頭,說明超出了本區域
		}
		else{
			cout << "!!!無法處理的情況,可能邊界數量超過2條!!!";
		}*/
		for each (contour *boundary in explorePossibleBoundary)
		{
			if (boundary->distance < U[exploreOrigenVIdx]){//距離較小的等高線
				if (debug)
					cout << "較小等高線" << boundary->id<<"距離:"<<boundary->distance;
				bool meetBoundary = boundary->contact_asSmall((*vv_cwit).idx());//是較小等高線的小源頭之一,說明該點應該位於上一個區域
				if (meetBoundary){
					boundaryContourIdx.insert(boundary->id);
					if (debug){
						cout << "meetboundary!!!";
						pointFromTwoSource pfts = boundary->get_contact_asSmall_pfts((*vv_cwit).idx());
						cout << "contact_asSmall found vidx:" << pfts.sourceIdx_Small << "U[sourceIdx_Small]=" << U[pfts.sourceIdx_Small] << "sourceIdx_Large:" << pfts.sourceIdx_Large << "U[sourceIdx_Large]=" << U[pfts.sourceIdx_Large] << ";";
					}
					//break;
					if (!hasMeetBoundary){
						hasMeetBoundary = true;
					}
				}
				
			}
			else{
				bool meetBoundary = boundary->contact_asLarge((*vv_cwit).idx());//是較大登高線的大源頭之一,說明該點應該位於下一個區域
				if (debug)
					cout << "較大等高線" << boundary->id << "距離:" << boundary->distance;
				if (meetBoundary){
					boundaryContourIdx.insert(boundary->id);
					if (debug)
						cout << "meetboundary!!!";
					//break;
					if (!hasMeetBoundary){
						hasMeetBoundary = true;
					}
				}
			}
		}
		if (debug)
			cout << endl;
		if (!hasMeetBoundary)
			exploreVertex((*vv_cwit).idx(), regionVertexIdx, boundaryContourIdx,debug);

	}
}
void polarMap::exploreFace(int contourIdx, int nowFaceIdx, MyMesh::HalfedgeHandle *fromEdge, vector<int> &allFaceIdx){
	if (!contain(allFaceIdx, nowFaceIdx)){
		allFaceIdx.push_back(nowFaceIdx);
		MyMesh::FaceHandle nowFace(nowFaceIdx);
		vector<contour*> boundaries;
		contour *c1 = (contours->operator[](contourIdx));
		boundaries.assign(c1->nextContours.begin(), c1->nextContours.end());
		boundaries.push_back(c1);
		
		for (MyMesh::FaceHalfedgeCWIter fh_cwit = mesh.fh_cwbegin(nowFace); fh_cwit != mesh.fh_cwend(nowFace); fh_cwit++){
			if (fromEdge != NULL){
				if (*fh_cwit == *fromEdge){//如果half edge 是
					continue;
				}
			}
			bool meetboundary=false;
			for each (contour *c in boundaries)
			{
				if(c->contact(mesh.from_vertex_handle(*fh_cwit).idx(), mesh.to_vertex_handle(*fh_cwit).idx())){
					meetboundary = true;
				}
			}
			if (!meetboundary){
				MyMesh::FaceHandle nextFace = mesh.opposite_face_handle(*fh_cwit);
				if (!contain(allFaceIdx, nextFace.idx())){
					exploreFace(contourIdx,nextFace.idx(),&mesh.opposite_halfedge_handle(*fh_cwit),allFaceIdx);
				}
			}
		}
	}
}
void polarMap::removeIncorrectExtremumVertexs(){//用來刪掉findExtremumVertexs找出來的並不是真正的
	//cout << "removeIncorrectExtremumVertexs開始:" << endl;
	/*cout << "localExtremum_idx:";
	for each (int idx in localExtremum_idx)
	{
		cout << idx << ",";
	}
	cout << endl;*/
	vector<int> tempEVs;
	tempEVs.assign(localExtremum_idx.begin(), localExtremum_idx.end());//複製初始localExtremum_idx
	while (tempEVs.size() > 0){

		int nowVIdx = tempEVs[0];
		//cout << "點" << nowVIdx << ":";
		tempEVs.erase(tempEVs.begin() + 0);

		vector<int> candidates;
		candidates.push_back(nowVIdx);

		vector<int> region_vidx;
		set<int> boundary_index;
		exploreOrigenVIdx = nowVIdx;
		explorePossibleBoundary = getPossibleBoundary(U[nowVIdx],1);
		/*cout << "explorePossibleBoundary:";
		for each (contour* b in explorePossibleBoundary)
		{
			cout << b->id << ",";
		}
		cout << "     ";*/
		exploreVertex(nowVIdx, region_vidx, boundary_index,nowVIdx == 1830);
		if (nowVIdx==1830)
		{
			for each (int v in region_vidx)
			{
				debugVertexIndex.insert(v);
			}
		}
		//移除同一區域重複的頂點并加入候選人

		for (int i = 0; i < tempEVs.size();i++)
		{
			int vidx = tempEVs[i];
			if (contain(region_vidx, vidx)){
				tempEVs.erase(tempEVs.begin() + i--);
				//這裡是插入排序法,結束時candidates會是一個U值從大到小的隊列
				bool inserted = false;
				for (int i = 0; i < candidates.size(); i++){
					if (U[candidates[i]] < U[vidx]){//第一次找到比自己U值小的候選人
						candidates.insert(candidates.begin() + i, vidx);
						inserted = true;
						break;
					}
				}
				if (!inserted)//這說明在candidates里找不到比vidx的U值更小的elem,換言之vidx的U值是目前最小
					candidates.push_back(vidx);
			}
		}
		bool ExtremumArea = (boundary_index.size() == 1);
		if (boundary_index.size() > 1){//即使多於1也可能是因為getPossibleBoundary的額外extraLayer邊界被不小心碰觸到
			//這裡判斷是不是所有接觸到的邊界都比極值點遠或近,如果該極值點比所有邊界都遠或都近的話那依舊是極值點
			float max_dist = -1;
			float min_dist = 1000000;
			for each (int bid in boundary_index)
			{
				if (contours->operator[](bid-1)->distance>max_dist)
				{
					max_dist = contours->operator[](bid - 1)->distance;
				}
				if (contours->operator[](bid - 1)->distance<min_dist)
				{
					min_dist = contours->operator[](bid - 1)->distance;
				}
			}
			//極值點比大小的方式就只要和所有邊界的最大值和最小值比就知道極值點是否夾在中間了
			if (U[nowVIdx] > max_dist || U[nowVIdx] < min_dist){
				ExtremumArea = true;
			}
		}
		//從候選人中找到離原點最遠的頂點作為真正的局部極值點
		//if (boundary_index.size() == 1){//只被一條等高線圍繞正確的case
		if (ExtremumArea){//如果從上面的判斷得知這是個極值區域的話,判斷就簡單了
			if (candidates.size() > 1){//這說明同一區域存在複數'局部極值點',而我們只需要一個局部極值點
				for (int i = 1; i < candidates.size(); i++){
					int fakeVertex = candidates[i];
					remove(localExtremum_idx, fakeVertex);
					//cout << "boundary_index size1 remove:" << fakeVertex;
				}
			}
		}
		else{//不正確的case
			for (int i = 0; i < candidates.size(); i++){
				int fakeVertex = candidates[i];
				remove(localExtremum_idx, fakeVertex);
				//cout << "size0 remove:" << fakeVertex;
			}
		}
	}
	cout << "removeIncorrectExtremumVertexs結束:" << endl;
	cout << "localExtremum_idx:";
	for each (int idx in localExtremum_idx)
	{
		cout << idx << ",";
	}
	cout << endl;
}
void polarMap::calU(int vidx, vector<MyMesh::VertexHandle> &candidates,bool debug){
	MyMesh::VertexHandle center = MyMesh::VertexHandle(vidx);
	MyMesh::VertexVertexCWIter vvc_it = mesh.vv_cwbegin(center);
	MyMesh::VertexHandle lastPoint =*(mesh.vv_ccwbegin(center));
	MyMesh::VertexHandle nowPoint;
	ostringstream ostm;
	if (debug)
		ostm <<endl<< "處理點:" << vidx<<"theta:"<<theta[vidx]<<endl;
	bool update = false;
	for (vvc_it++; vvc_it != mesh.vv_cwend(center); vvc_it++){
		nowPoint = *vvc_it;
		int k = lastPoint.idx();
		int j = nowPoint.idx();
		if (debug)
			ostm << ">k為:" << k << "j為:" << j;
		OpenMesh::Vec3f ek = mesh.point(lastPoint) - mesh.point(center);
		OpenMesh::Vec3f ej = mesh.point(nowPoint) - mesh.point(center);
		OpenMesh::Vec3f ejk = mesh.point(lastPoint) - mesh.point(nowPoint);
		lastPoint = nowPoint;
		float A = OpenMesh::cross(ej, ek).length();
		float ejk_2 = OpenMesh::dot(ejk, ejk);
		float H = sqrtf((ejk_2-powf(U[j]-U[k],2))*(powf(U[j]+U[k],2)-ejk_2));
		float xj = (A*(ejk_2 + powf(U[k], 2) - powf(U[j], 2)) + OpenMesh::dot(ek, ejk)*H) / (2 * A*ejk_2);
		float xk = (A*(ejk_2 + powf(U[j], 2) - powf(U[k], 2)) - OpenMesh::dot(ej, ejk)*H) / (2 * A*ejk_2);
		float newUi;
		float newTheta;
		int i = vidx;
		//cout << "U[" << i << "]:" << U[i]<<" ";
		/*if (vidx==1){
			cout << "j:" << j << "theta:" << theta[j] << " k:" << k << "theta:" << theta[k]<<endl;
		}*/
		if (xj > 0 && xk > 0){
			if (debug)
				ostm << "_xjxk:";
			newUi=(xj*ej + xk*ek).length();
			if (debug)
				ostm << "newUi:" << newUi;
			anglepair apair = cal_phki_phij(ek, ej, U[k], U[j]);
			if (debug){
				cal_phki_phij(ek, ej, U[k], U[j], debug,ostm);
				ostm << "~~~ek:(" << ek[0] << "," << ek[1] << "," << ek[2] << ") ej:(" << ej[0] << "," << ej[1] << "," << ej[2] << ") Uk:" << U[k] << " Uj:" << U[j] << endl;
			}
			float a = apair.ph_ij / (apair.ph_ij + apair.ph_ki);
			//cout << "calU"<<" a:" << a<<" ";
			if (debug){
				ostm << "calTheta:a:" << a << " theta[" << j << "]:" << theta[j] << " theta[" << k << "]:" << theta[k]<<" phij:"<<apair.ph_ij<<" phki"<<apair.ph_ki;
			}
			if (abs(theta[j] - theta[k]) > M_PI){
				if (theta[j] > theta[k]){
					newTheta = (1 - a)*theta[j] + a*(theta[k] + 2 * M_PI);

				}
				else{
					newTheta = (1 - a)*(theta[j] + 2 * M_PI) + a*theta[k];

				}
				if (newTheta >= 2 * M_PI){
					newTheta -= 2 * M_PI;
				}
			}
			else{
				newTheta = (1 - a)*theta[j] + a*theta[k];
			}
			/*if (abs(theta[j] - theta[k]) > M_PI){
				newTheta  += M_PI;
				while (newTheta > 2 * M_PI){
					newTheta -= (2 * M_PI);
				}
			}*/
			if (debug)
				ostm << " newTheta:" << newTheta;
		}
		else{
			float minLk = U[k] + ek.length();
			float minLj = U[j] + ej.length();
			if (minLk < minLj){
				newUi = minLk;
				newTheta = theta[k];
			}
			else{
				newUi = minLj;
				newTheta = theta[j];
			}
			if (debug){
				ostm << "else:newUi:" << newUi;
				ostm << " newTheta:" << newTheta;
			}
		}
		
		if (newUi < U[i]){
			if (debug)
				ostm << "set!" << endl;
			U[i] = newUi;
			theta[i] = newTheta;
			update = true;
			/*if (contain(debug_pointIdx, 8, i)){
				
			}*/
		}
		if (debug)
			ostm << endl;
		//OpenMesh::Vec3f xk = / 2 * A* (OpenMesh::dot(ejk, ejk));
	}//結束后nowPoint應該是iter_end 


	lastPoint = nowPoint;
	nowPoint = (*mesh.vv_cwbegin(center));
	int k = lastPoint.idx();
	int j = nowPoint.idx();
	if (debug)
		ostm << "->k為:" << k << "j為:" << j;
	OpenMesh::Vec3f ek = mesh.point(lastPoint) - mesh.point(center);
	OpenMesh::Vec3f ej = mesh.point(nowPoint) - mesh.point(center);
	OpenMesh::Vec3f ejk = mesh.point(lastPoint) - mesh.point(nowPoint);
	lastPoint = nowPoint;
	float A = OpenMesh::cross(ej, ek).length();
	float ejk_2 = OpenMesh::dot(ejk, ejk);
	float H = sqrtf((ejk_2 - powf(U[j] - U[k], 2))*(powf(U[j] + U[k], 2) - ejk_2));
	float xj = (A*(ejk_2 + powf(U[k], 2) - powf(U[j], 2)) + OpenMesh::dot(ek, ejk)*H) / (2 * A*ejk_2);
	float xk = (A*(ejk_2 + powf(U[j], 2) - powf(U[k], 2)) - OpenMesh::dot(ej, ejk)*H) / (2 * A*ejk_2);
	float newUi;
	float newTheta;
	int i = vidx;//(*vvc_it).idx();
	if (xj > 0 && xk > 0){
		if (debug)
			ostm << "_xjxk:";
		newUi = (xj*ej + xk*ek).length();
		if (debug)
			ostm << "newUi:" << newUi;
		anglepair apair = cal_phki_phij(ek, ej, U[k], U[j]);
		if (debug){
			cal_phki_phij(ek, ej, U[k], U[j], debug, ostm);
			ostm << "~~~ek:(" << ek[0] << "," << ek[1] << "," << ek[2] << ") ej:(" << ej[0] << "," << ej[1] << "," << ej[2] << ") Uk:" << U[k] << " Uj:" << U[j] << endl;
		}
		float a = apair.ph_ij / (apair.ph_ij + apair.ph_ki);
		if (abs(theta[j] - theta[k]) > M_PI){
			if (theta[j] > theta[k]){
				newTheta = (1 - a)*theta[j] + a*(theta[k]+2*M_PI);

			}
			else{
				newTheta = (1 - a)*(theta[j] + 2 * M_PI) + a*theta[k];
				
			}
			if (newTheta >= 2 * M_PI){
				newTheta -= 2 * M_PI;
			}
		}
		else{
			newTheta = (1 - a)*theta[j] + a*theta[k];
		}
		/*
		if (abs(theta[j] - theta[k]) > M_PI){
			newTheta += M_PI;
			while (newTheta > 2 * M_PI){
				newTheta -= (2 * M_PI);
			}
		}*/
		//if()
		if (debug)
			ostm << " newTheta:" << newTheta;
	}
	else{
		float minLk = U[k] + ek.length();
		float minLj = U[j] + ej.length();
		if (minLk < minLj){
			newUi = minLk;
			newTheta = theta[k];
		}
		else{
			newUi = minLj;
			newTheta = theta[j];
		}
		if (debug){
			ostm << "else:newUi:" << newUi;
			ostm << " newTheta:" << newTheta;
		}
	}
	if (newUi < U[i]){
		if (debug)
			ostm << "set!" << endl;
		U[i] = newUi;
		theta[i] = newTheta;
		update = true;
	}
	if (debug)
		ostm << "最後theta為:" << theta[i] << endl;
	if (update){
		cout << ostm.str();
		candidates.push_back(center);
	}
}
void polarMap::calU(int vidx, vector<MyMesh::VertexHandle> &candidates, int debugIdx){
	MyMesh::VertexHandle center = MyMesh::VertexHandle(vidx);
	MyMesh::VertexVertexCWIter vvc_it = mesh.vv_cwbegin(center);
	MyMesh::VertexHandle lastPoint = *(mesh.vv_ccwbegin(center));
	MyMesh::VertexHandle nowPoint;
	ostringstream ostm;
	bool debug = vidx == debugIdx;
	if (debug)
		ostm << endl << "處理點:" << vidx << "theta:" << theta[vidx] << endl;
	bool update = false;
	for (vvc_it++; vvc_it != mesh.vv_cwend(center); vvc_it++){
		nowPoint = *vvc_it;
		int k = lastPoint.idx();
		int j = nowPoint.idx();

		
		OpenMesh::Vec3f ek = mesh.point(lastPoint) - mesh.point(center);
		OpenMesh::Vec3f ej = mesh.point(nowPoint) - mesh.point(center);
		OpenMesh::Vec3f ejk = mesh.point(lastPoint) - mesh.point(nowPoint);
		lastPoint = nowPoint;
		float A = OpenMesh::cross(ej, ek).length();
		float ejk_2 = OpenMesh::dot(ejk, ejk);
		float H = sqrtf((ejk_2 - powf(U[j] - U[k], 2))*(powf(U[j] + U[k], 2) - ejk_2));
		float xj = (A*(ejk_2 + powf(U[k], 2) - powf(U[j], 2)) + OpenMesh::dot(ek, ejk)*H) / (2 * A*ejk_2);
		float xk = (A*(ejk_2 + powf(U[j], 2) - powf(U[k], 2)) - OpenMesh::dot(ej, ejk)*H) / (2 * A*ejk_2);
		float newUi;
		float newTheta;
		int i = vidx;
		//cout << "U[" << i << "]:" << U[i]<<" ";
		/*if (vidx==1){
		cout << "j:" << j << "theta:" << theta[j] << " k:" << k << "theta:" << theta[k]<<endl;
		}*/
		if (xj > 0 && xk > 0){
			//if (debug)
			//	ostm << "_xjxk:";
			newUi = (xj*ej + xk*ek).length();
			//if (debug)
			//	ostm << "newUi:" << newUi;
			anglepair apair = cal_phki_phij(ek, ej, U[k], U[j]);
			//if (debug){
			//	cal_phki_phij(ek, ej, U[k], U[j], debug, ostm);
			//	ostm << "~~~ek:(" << ek[0] << "," << ek[1] << "," << ek[2] << ") ej:(" << ej[0] << "," << ej[1] << "," << ej[2] << ") Uk:" << U[k] << " Uj:" << U[j] << endl;
			//}
			float a = apair.ph_ij / (apair.ph_ij + apair.ph_ki);
			//cout << "calU"<<" a:" << a<<" ";
			//if (debug){
			//	ostm << "calTheta:a:" << a << " theta[" << j << "]:" << theta[j] << " theta[" << k << "]:" << theta[k] << " phij:" << apair.ph_ij << " phki" << apair.ph_ki;
			//}
			if (abs(theta[j] - theta[k]) > M_PI){
				if (theta[j] > theta[k]){
					newTheta = (1 - a)*theta[j] + a*(theta[k] + 2 * M_PI);

				}
				else{
					newTheta = (1 - a)*(theta[j] + 2 * M_PI) + a*theta[k];

				}
				if (newTheta >= 2 * M_PI){
					newTheta -= 2 * M_PI;
				}
			}
			else{
				newTheta = (1 - a)*theta[j] + a*theta[k];
			}
			/*if (abs(theta[j] - theta[k]) > M_PI){
			newTheta  += M_PI;
			while (newTheta > 2 * M_PI){
			newTheta -= (2 * M_PI);
			}
			}*/
			//if (debug)
			//	ostm << " newTheta:" << newTheta;
		}
		else{
			float minLk = U[k] + ek.length();
			float minLj = U[j] + ej.length();
			if (minLk < minLj){
				newUi = minLk;
				newTheta = theta[k];
			}
			else{
				newUi = minLj;
				newTheta = theta[j];
			}
			//if (debug){
			//	ostm << "else:newUi:" << newUi;
			//	ostm << " newTheta:" << newTheta;
			//}
		}

		if (newUi < U[i]){
			if (k == debugIdx || j == debugIdx){
				ostm << "作為來源 k=" << k << "j=" << j << " U[" << k << "]=" << U[k] << " U[" << j << "]=" << U[j] << " 目標vidx=" << vidx;
				ostm << "old:" << U[i] << ">new:" << newUi << endl;
			}
			if (debug){
				ostm << ">k為:" << k << "j為:" << j << " U[" << k << "]=" << U[k] << " U[" << j << "]=" << U[j];
				ostm << "old:" << U[i] << ">new:" << newUi << endl;
			}
			U[i] = newUi;
			theta[i] = newTheta;
			update = true;
			/*if (contain(debug_pointIdx, 8, i)){

			}*/
		}
		//OpenMesh::Vec3f xk = / 2 * A* (OpenMesh::dot(ejk, ejk));
	}//結束后nowPoint應該是iter_end 


	lastPoint = nowPoint;
	nowPoint = (*mesh.vv_cwbegin(center));
	int k = lastPoint.idx();
	int j = nowPoint.idx();


	OpenMesh::Vec3f ek = mesh.point(lastPoint) - mesh.point(center);
	OpenMesh::Vec3f ej = mesh.point(nowPoint) - mesh.point(center);
	OpenMesh::Vec3f ejk = mesh.point(lastPoint) - mesh.point(nowPoint);
	lastPoint = nowPoint;
	float A = OpenMesh::cross(ej, ek).length();
	float ejk_2 = OpenMesh::dot(ejk, ejk);
	float H = sqrtf((ejk_2 - powf(U[j] - U[k], 2))*(powf(U[j] + U[k], 2) - ejk_2));
	float xj = (A*(ejk_2 + powf(U[k], 2) - powf(U[j], 2)) + OpenMesh::dot(ek, ejk)*H) / (2 * A*ejk_2);
	float xk = (A*(ejk_2 + powf(U[j], 2) - powf(U[k], 2)) - OpenMesh::dot(ej, ejk)*H) / (2 * A*ejk_2);
	float newUi;
	float newTheta;
	int i = vidx;//(*vvc_it).idx();
	if (xj > 0 && xk > 0){
		//if (debug)
		//	ostm << "_xjxk:";
		newUi = (xj*ej + xk*ek).length();
		//if (debug)
		//	ostm << "newUi:" << newUi;
		anglepair apair = cal_phki_phij(ek, ej, U[k], U[j]);
		//if (debug){
		//	cal_phki_phij(ek, ej, U[k], U[j], debug, ostm);
		//	ostm << "~~~ek:(" << ek[0] << "," << ek[1] << "," << ek[2] << ") ej:(" << ej[0] << "," << ej[1] << "," << ej[2] << ") Uk:" << U[k] << " Uj:" << U[j] << endl;
		//}
		float a = apair.ph_ij / (apair.ph_ij + apair.ph_ki);
		if (abs(theta[j] - theta[k]) > M_PI){
			if (theta[j] > theta[k]){
				newTheta = (1 - a)*theta[j] + a*(theta[k] + 2 * M_PI);

			}
			else{
				newTheta = (1 - a)*(theta[j] + 2 * M_PI) + a*theta[k];

			}
			if (newTheta >= 2 * M_PI){
				newTheta -= 2 * M_PI;
			}
		}
		else{
			newTheta = (1 - a)*theta[j] + a*theta[k];
		}
		/*
		if (abs(theta[j] - theta[k]) > M_PI){
		newTheta += M_PI;
		while (newTheta > 2 * M_PI){
		newTheta -= (2 * M_PI);
		}
		}*/
		//if()
		//if (debug)
		//	ostm << " newTheta:" << newTheta;
	}
	else{
		float minLk = U[k] + ek.length();
		float minLj = U[j] + ej.length();
		if (minLk < minLj){
			newUi = minLk;
			newTheta = theta[k];
		}
		else{
			newUi = minLj;
			newTheta = theta[j];
		}
		//if (debug){
		//	ostm << "else:newUi:" << newUi;
		//	ostm << " newTheta:" << newTheta;
		//}
	}
	if (newUi < U[i]){
		if (k==debugIdx||j==debugIdx){
			ostm << "作為來源 k=" << k << "j=" << j << " U[" << k << "]=" << U[k] << " U[" << j << "]=" << U[j] << " 目標vidx=" << vidx;
			ostm << "old:" << U[i] << ">new:" << newUi << endl;
		}
		if (debug){
			ostm << "->k為:" << k << "j為:" << j << " U[" << k << "]=" << U[k] << " U[" << j << "]=" << U[j];
			ostm << "old:" << U[i] << ">new:" << newUi << endl;
		}
		U[i] = newUi;
		theta[i] = newTheta;
		update = true;
	}
	//if (debug)
	//	ostm << "最後theta為:" << theta[i] << endl;
	if (update){
		cout << ostm.str();
		candidates.push_back(center);
	}
}
void polarMap::showPoints(int range){
	if (range > mesh.n_vertices()){
		range = mesh.n_vertices();
	}
	for (int i = 0; i < range; i++){
		cout << "頂點" << i << ": U:" << U[i] << " theta/pi:" << theta[i]/M_PI << "|  "<<endl;
	}
}
void polarMap::showPoints(vector<int> vList){
	for (int i=0; i < vList.size(); i++){
		int vIdx = vList[i];
		OpenMesh::VertexHandle vh = OpenMesh::VertexHandle(vIdx);
		OpenMesh::Vec3f pos = mesh.point(vh);
		float len = sqrtf(powf(pos[0], 2) + powf(pos[2], 2));
		float angle = acosf(pos[2]/len);
		cout << "vertex" << vIdx << ":" << pos[0] << "," << pos[1] << "," << pos[2] <<"U:"<<U[vIdx]<<"theta:"<<theta[vIdx]/M_PI<< " cal angle:"<<angle/M_PI<<endl;
	}
}
/*
float polarMap::computeDistance(int vidx){
	MyMesh::VertexHandle center = MyMesh::VertexHandle(vidx);
	MyMesh::VertexHandle lastPoint = *(mesh.vv_cwend(center));
	float minDist = BIG_NUMBER;
	for (MyMesh::VertexVertexCWIter vvc_it = mesh.cvv_cwbegin(center); vvc_it != mesh.vv_cwend(center); vvc_it++){
		MyMesh::VertexHandle nowPoint = *vvc_it;
		int k = lastPoint.idx;
		int j = nowPoint.idx;
		OpenMesh::Vec3f ek = mesh.point(lastPoint) - mesh.point(center);
		OpenMesh::Vec3f ej = mesh.point(nowPoint) - mesh.point(center);
		OpenMesh::Vec3f ejk = mesh.point(lastPoint) - mesh.point(nowPoint);
		lastPoint = nowPoint;
		float A = OpenMesh::cross(ej, ek).length();
		float ejk_2 = OpenMesh::dot(ejk, ejk);
		float H = sqrtf((ejk_2 - powf(U[j] - U[k], 2))*(powf(U[j] + U[k], 2) - ejk_2));
		float xj = (A*(ejk_2 + powf(U[k], 2) - powf(U[j], 2)) + OpenMesh::dot(ek, ejk)*H) / (2 * A*ejk_2);
		float xk = (A*(ejk_2 + powf(U[j], 2) - powf(U[k], 2)) - OpenMesh::dot(ej, ejk)*H) / (2 * A*ejk_2);
		float newUi;
		if (xj > 0 && xk > 0){
			newUi = (xj*ej + xk*ek).length;
		}
		else{
			float minLk = U[k] + ek.length();
			float minLj = U[j] + ej.length();
			if (minLk < minLj){
				newUi = minLk;
			}
			else{
				newUi = minLj;
			}
		}
		if (newUi < minDist){
			minDist = newUi;
		}
	}
	return minDist;
}*/
MyMesh::VertexHandle polarMap::popSmallestNode(vector<MyMesh::VertexHandle> &list){
	int sindex = 0;
	MyMesh::VertexHandle smallest = list[0];
	for (int i = 1; i < list.size(); i++){
		int index = list[i].idx();
		if (U[index] < U[smallest.idx()]){
			smallest = list[i];
			sindex = i;
		}
	}
	list.erase(list.begin()+sindex);
	return smallest;
	
}
vector<MyMesh::VertexHandle> polarMap::neighbour(MyMesh::VertexHandle vertex){
	vector<MyMesh::VertexHandle> result;
	for (MyMesh::VertexVertexCWIter vvc_it = mesh.cvv_cwbegin(vertex); vvc_it != mesh.vv_cwend(vertex); vvc_it++){
		result.push_back(*vvc_it);
	}
	return result;
}
void polarMap::createMap(){
	cout << "開始create map"<<endl;
	U = new float[mesh.n_vertices()];
	cout << "mesh.n_vertices():" << mesh.n_vertices();
	for (int i = 0; i < mesh.n_vertices(); i++){
		U[i] = BIG_NUMBER;
	}
	theta = new float[mesh.n_vertices()];
	for (int i = 0; i < mesh.n_vertices(); i++){
		theta[i] = 0;
	}
	initStartRing();
	cout << "用來作為基本軸的頂點idx為:" << zero_anix_vertexIdx<<endl;
	vector<MyMesh::VertexHandle> candidates;
	MyMesh::VertexHandle startPoint = MyMesh::VertexHandle(poleIdx);
	//for (MyMesh::VertexVertexCWIter vvc_it = mesh.cvv_cwbegin(startPoint); vvc_it != mesh.vv_cwend(startPoint); vvc_it++){
	//	candidates.push_back(*vvc_it);
	//}
	candidates = neighbour(startPoint);
	int count = 0;
	while (candidates.size() > 0){
		MyMesh::VertexHandle node= popSmallestNode(candidates);
		
		vector<MyMesh::VertexHandle> ring = neighbour(node);
		for (int i = 0; i < ring.size(); i++){
			MyMesh::VertexHandle vertex = ring[i];
			calU(vertex.idx(), candidates, -1);

		}
		count++;
		//cout << "迭代后candidates的數量:" + candidates.size();
	}
	//cout << "迭代次數:" + count;
	//printf("count %d",count);
	//for (int i = 0; i < mesh.n_vertices(); i++){
	//	cout << "vertice" << i <<"U:"<<U[i] << " theta:" << theta[i] << endl;
	//}
	vector<int> vl;
	for (int i = 0; i < 96; i++){
		vl.push_back(show_vertexList[i]);
	}
	//showPoints(vl);
	farthestIdx = 0;
	for (int v = 0; v < mesh.n_vertices(); v++){
		if (U[farthestIdx] < U[v]){
			farthestIdx = v;
		}
	}

	findExtremumVertexs();
}
OpenMesh::Vec2f normalizePos(float* border, OpenMesh::Vec2f ori){
	float width = border[2] - border[0];
	float height = border[3] - border[1];
	return OpenMesh::Vec2f((ori[0] - border[0]) / width, (ori[1] - border[1]) / height);
}
OpenMesh::Vec3f polarMap::calPercentagetPosAtPath(vector<OpenMesh::Vec3f> path, float percentage){
	float totalLength = 0;
	for (int i = 0; i < path.size() - 1; i++){
		totalLength += (path[i + 1] - path[i]).length();
	}
	float nowLength = 0;
	for (int i = 0; i < path.size() - 1; i++){
		float len = (path[i + 1] - path[i]).length();
		if ((nowLength + len) / totalLength >= percentage){//點就在這個範圍內
			float sub = ((percentage*totalLength) - nowLength) / len;
			return path[i] + (path[i + 1] - path[i])*sub;
		}
		else
		{
			nowLength += len;
		}
	}
	cout << "percent:" << percentage << "搜尋不到點!" << endl;
	return path.back();
}
string int2str(int i) {
	string s;
	stringstream ss(s);
	ss << i;
	return ss.str();
}
void polarMap::drawMeshToImage(){
	painter = new draw2dConnectMap(801,801);
	//計算border
	OpenMesh::Vec2f *euclidList = new OpenMesh::Vec2f[mesh.n_vertices()];
	float border[4] = {0,0,0,0};
	for (int v = 0; v < mesh.n_vertices(); v++){
		euclidList[v] = getPointFrom(v).toEuclid2d();
		if (euclidList[v][0] < border[0]){
			border[0] = euclidList[v][0];
		}
		if (euclidList[v][1] < border[1]){
			border[1] = euclidList[v][1];
		}
		if (euclidList[v][0] > border[2]){
			border[2] = euclidList[v][0];
		}
		if (euclidList[v][1] > border[3]){
			border[3] = euclidList[v][1];
		}
	}
	//添加空白
	float center[2] = {(border[0]+border[2])/2,(border[1]+border[3])/2};
	float radius[2] = { border[2] - center[0], border[3] - center[1] };
	radius[0] *= 1.1;//寬製造0.1倍的空白
	radius[1] *= 1.1;//高製造0.1倍的空白
	border[0] = center[0] - radius[0];//修改border
	border[1] = center[1] - radius[1];
	border[2] = center[0] + radius[0];
	border[3] = center[1] + radius[1];
	for (int f = 0; f < mesh.n_faces();f++){
		OpenMesh::FaceHandle fh(f);
		vector<float> list;
		bool containEnd = false;
		for (MyMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it != mesh.fv_end(fh);fv_it++){
			//OpenMesh::Vec2f pos= getPointFrom((*fv_it).idx()).toEuclid2d();
			OpenMesh::Vec2f pos = euclidList[(*fv_it).idx()];
			pos = normalizePos(border, pos);//將點的x,y值歸一化到0到1之間
			list.push_back(pos[0]);
			list.push_back(pos[1]);
			if ((*fv_it).idx() == farthestIdx){
				containEnd = true;
			}
		}
		float* ldata= list.data();
		if(!containEnd)
			painter->drawTriangle(threePoints(ldata));
	}
	//OpenMesh::Vec2f traget_point = normalizePos(border, euclidList[zero_anix_vertexIdx]);
	//painter->drawPoint(traget_point[0], traget_point[1]);
	for (int i = 0; i < debug_num; i++){
		OpenMesh::Vec2f traget_point = normalizePos(border, euclidList[debug_pointIdx[i]]);
		painter->drawPoint(traget_point[0], traget_point[1]);
	}
	if (Image_show_vertex_idx){
		for (int v = 0; v < mesh.n_vertices(); v++){
			string no;
			no = int2str(v);
			OpenMesh::Vec2f traget_point = normalizePos(border, euclidList[v]);
			//painter->drawPoint(traget_point[0], traget_point[1]);
			painter->writeTxt(traget_point[0], traget_point[1], no);
		}
	}
	painter->showMap();
	
	delete []euclidList;
}
void polarMap::drawContourToImage(vector<PointCurve> contours){
	draw2dConnectMap painter = draw2dConnectMap(1001, 1001);
	int contourNum = contours.size();
	float radius = 0.45;
	//OpenMesh::Vec2f center = OpenMesh::Vec2f(radius, radius);
	for (int i = 0; i < contourNum; i++){
		float now_radius = 0.45 *((float)(i + 1) /(float)contourNum);
		if (i==14)
			cout << "contour"<<i<<" 由" << contours[i].points.size() << "点组成   ";
		for (int p = 0; p < contours[i].points.size()-1; p++){
			//cout << "center:(" << center[0] << "," << center[1] << ")" << endl;
			if (i == 14)
				cout <<"angle:"<<contours[i].points[p].angle << "x:" << cosf(2 * M_PI*contours[i].points[p].angle) << "y:" << sin(2 * M_PI*contours[i].points[p].angle) << " "<<endl;
			OpenMesh::Vec2f point1(cosf(2 * M_PI*contours[i].points[p].angle)*now_radius, sinf(2 * M_PI*contours[i].points[p].angle)*now_radius);
			OpenMesh::Vec2f point2(cosf(2 * M_PI*contours[i].points[p + 1].angle)*now_radius, sinf(2 * M_PI*contours[i].points[p + 1].angle)*now_radius);
			//point1 = point1 + center;
			//point2 = point2 + center;
			//cout << "畫線 (" << point1[0]+radius << "," << point1[1]+radius << ")" << "to" << "(" << point2[0]+radius << "," << point2[1]+radius << ") "<<endl;
			
			painter.drawLine(point1[0]+radius,point1[1]+radius,point2[0]+radius,point2[1]+radius);
			painter.drawPoint(point1[0] + radius, point1[1] + radius);
		}
		OpenMesh::Vec2f point_end(cosf(2 * M_PI*contours[i].points[contours[i].points.size() - 1].angle)*now_radius, sinf(2 * M_PI*contours[i].points[contours[i].points.size() - 1].angle)*now_radius);
		//OpenMesh::Vec2f point2(cosf(2 * M_PI*contours[i].points[p + 1].angle)*now_radius, sinf(2 * M_PI*contours[i].points[p + 1].angle)*now_radius);
		painter.drawPoint(point_end[0], point_end[1]);
		if (i==14)
			cout << endl;
	}
	painter.showMap();
}
int polarMap::getFasthestIdx(){
	return farthestIdx;
}
struct indexCurve{
	vector<int> indexs;
};
bool connect(MyMesh mesh, int vidx1, int vidx2){
	OpenMesh::VertexHandle vh = OpenMesh::VertexHandle(vidx1);
	for (MyMesh::VertexVertexCWIter vc_it = mesh.vv_cwbegin(vh); vc_it != mesh.vv_cwend(vh);vc_it++){
		if (vidx2==(*vc_it).idx()){
			return true;
		}
	}
	return false;
}
int indexOf(vector<pointFromTwoSource> list, int vidx, bool small){
	if (small){
		for(int i = 0; i<list.size();i++){
			pointFromTwoSource pfts = list[i];
			if (pfts.sourceIdx_Small == vidx){
				return i;
			}
		}
		return -1;
	}
	else{
		for (int i = 0; i<list.size(); i++){
			pointFromTwoSource pfts = list[i];
			if (pfts.sourceIdx_Large == vidx){
				return i;
			}
		}
		return -1;
	}
}
int indexOf(vector<pointFromTwoSource> list, int vidx1, int vidx2){
	for(int i = 0; i < list.size();i++){
		pointFromTwoSource pfts = list[i];
		if (pfts.sourceIdx_Small == vidx1 &&pfts.sourceIdx_Large == vidx2){
			return i;
		}
	}
	return -1;
}
vector<pointFromTwoSource> getPFTSfromIdx(vector<pointFromTwoSource> list, int pidx, bool small){
	vector<pointFromTwoSource> result;
	for each(pointFromTwoSource pfts in list){
		if (small&&pfts.sourceIdx_Small == pidx){
			result.push_back(pfts);
		}
		else if (!small&&pfts.sourceIdx_Large == pidx){
			result.push_back(pfts);
		}
	}
	return result;
}
void removePFTS(vector<pointFromTwoSource> &list, pointFromTwoSource traget){
	for (int i = 0; i < list.size(); i++){
		if (list[i].sourceIdx_Small == traget.sourceIdx_Small && list[i].sourceIdx_Large == traget.sourceIdx_Large){
			list.erase(list.begin() + i);
			break;
		}
	}
}


vector<vector<OpenMesh::Vec3f>> polarMap::dividePointSet(MyMesh mesh,vector<pointFromTwoSource> pointSet){
	vector<vector<pointFromTwoSource>> groupSmall;
	vector<vector<OpenMesh::Vec3f>> *results=new vector<vector<OpenMesh::Vec3f>>();
	vector<pointFromTwoSource> pointSetS(pointSet);//make a pointSet copy
	/*cout << "pointSet length:" << pointSet.size();
	for each(pointFromTwoSource pfts in pointSet){
		cout << "point:" << pfts.pointPos[0] << "," << pfts.pointPos[1] << "," << pfts.pointPos[2] << " ";
	}*/
	contours = new vector<contour*>();
	while (pointSetS.size() > 0){
	/*bool found = false;
		
		vector<int> merge_idx;
		for each(vector<pointFromTwoSource> subgp in groupSmall){
			for each(pointFromTwoSource group_point in subgp){
				if (connect(mesh, point.sourceIdx_Small, group_point.sourceIdx_Small)){
					subgp.push_back(group_point);
					found = true;
				}
			}
		}
		if (found){
			groupSmall.push_back(vector<pointFromTwoSource>{point});
		}*/
		vector<pointFromTwoSource> newGroup;
		pointFromTwoSource nowpfts = pointSetS[0];
		int sourceIdx = nowpfts.sourceIdx_Small;
		MyMesh::VertexHandle vh(sourceIdx);
		MyMesh::HalfedgeHandle nowHalfEdge;
		for (MyMesh::VertexIHalfedgeCWIter vhe_it = mesh.vih_cwbegin(vh); vhe_it != mesh.vih_cwend(vh); vhe_it++){
			//cout << "from vertex" << mesh.from_vertex_handle(*vhe_it).idx();
			if (mesh.from_vertex_handle(*vhe_it).idx() == nowpfts.sourceIdx_Large){
				//cout << " get!!!";
				nowHalfEdge = *vhe_it;
				break;
			}
		}
		bool findNext = false;
		MyMesh::HalfedgeHandle firstHalfEdge = nowHalfEdge;//記錄起始半邊,用於遇到邊界的反向重組
		bool firstMeetBoundary = true;
		//cout << endl << "newGroup add:";

		do{
			findNext = false;
			newGroup.push_back(nowpfts);
			//if (mesh.is_boundary(mesh.edge_handle(nowHalfEdge))){//如果遇到了邊界,即過了這條邊不再有面存在
			if (mesh.is_boundary(mesh.opposite_halfedge_handle(nowHalfEdge))){
				//cout<<"找到boundary!" << endl;
				if (firstMeetBoundary){//如果是第一次遇到
					firstMeetBoundary = false;
					nowHalfEdge= mesh.opposite_halfedge_handle(firstHalfEdge);//設置當前半邊為起始半邊的反向,這樣會向一開始的反方向重組pfts
				}
				else{//遇到第二次代表重組完成了,一般來說只要有第一次就一定會有第二次
					break;
				}
			}
			//cout << "(" << nowpfts.sourceIdx_Small << "," << nowpfts.sourceIdx_Large << ")";
			removePFTS(pointSetS, nowpfts);
			MyMesh::FaceHandle nowFace = mesh.opposite_face_handle(nowHalfEdge);//找到當前halfEdge的反向面
			for (MyMesh::FaceHalfedgeCWIter fhe_it = mesh.fh_cwbegin(nowFace); fhe_it != mesh.fh_cwend(nowFace); fhe_it++){//找找這個面里有沒有另一半邊屬於已知的pointFromTwoSource
				int vidx1 = mesh.from_vertex_handle(*fhe_it).idx();
				int vidx2 = mesh.to_vertex_handle(*fhe_it).idx();
				int index1 = indexOf(pointSetS,vidx1,vidx2);//先試試from 放前面to 放後面能不能找到pointFormTwoSource
				if (index1 >= 0){//找到了
					findNext = true;
					nowpfts = pointSetS[index1];
					nowHalfEdge = *fhe_it;
					break;
				}
				else{//沒找到
					//試試看反過來找
					int temp = vidx1;
					int index2 = indexOf(pointSetS, vidx2, vidx1);
					if (index2 >= 0){
						findNext = true;
						nowpfts = pointSetS[index2];
						nowHalfEdge = *fhe_it;
						break;
					}
				}
			}


		} while (findNext);
		//int index = 0;
		//pointSetS.erase(pointSetS.begin());
		/*
		int pointSource_idx = pointSetS[0].sourceIdx_Small;//初始条件两个当前pointSource和之上一个pointSource 
		int lastSource_idx = pointSetS[0].sourceIdx_Small;
		int firstPoint_idx = -1;
		int nextPoint_index;
		bool getEndPointOfContour = false;
		do{
			nextPoint_index = -1;
			vector<pointFromTwoSource> pftsBuffers[3];
			int nowBuffer=0;
			bool last2next = true;
			//bool debug = (pointSource_idx == 2);
			OpenMesh::VertexHandle vh(pointSource_idx);
			vector<pointFromTwoSource> connecteds = getPFTSfromIdx(pointSetS, pointSource_idx, true);
			//先找出下一个source point,即和自身连接且距离最短的非上一个source point
			float minDist = 1000000;
			for (MyMesh::VertexVertexCWIter vv_cit = mesh.vv_cwbegin(vh); vv_cit != mesh.vv_cwend(vh); vv_cit++){
				int index = indexOf(pointSetS, (*vv_cit).idx(), true);
				if (index>=0){
					//OpenMesh::Vec3f pos = mesh.point(MyMesh::VertexHandle(pointSetS[index].sourceIdx_Small));
					//OpenMesh::Vec3f spos = mesh.point(vh);
					//float dist = (pos - spos).length();
					vector<pointFromTwoSource> vv_cSet = getPFTSfromIdx(pointSetS,(*vv_cit).idx(),true);
					for each(pointFromTwoSource pfts in vv_cSet){
						for each(pointFromTwoSource sub in connecteds){
							float dist = (pfts.pointPos - sub.pointPos).length();
								if (dist < minDist){
									if (firstPoint_idx == -1 && nextPoint_index != pointSetS[index].sourceIdx_Small){//如果这是这个圈的起始点
										lastSource_idx = nextPoint_index;//设置上个点索引值为第二小的点
									}
									nextPoint_index = pointSetS[index].sourceIdx_Small;
									minDist = dist;
								}
						}
					}
				}
			}
			if (nextPoint_index == -1){//找不到下一个相邻的source point 了意为这是最后一个点了
				getEndPointOfContour = true;
				nextPoint_index = firstPoint_idx;//把最后一个点的next idx设为起始点
			}
			for (MyMesh::VertexVertexCWIter vv_cit = mesh.vv_cwbegin(vh); vv_cit != mesh.vv_cwend(vh); vv_cit++){
				int index = indexOf(pointSetS, (*vv_cit).idx(), true);//搜寻small source point
				if ((*vv_cit).idx() == nextPoint_index){
					nowBuffer++;
					if (nowBuffer == 1){
						last2next = false;//先找到了next

						if (lastSource_idx == -1){//虽然不太可能发生,这种情况以为没有上一个节点,如果在圈的第一个点实在只能找到1个相连的sourcepoint就会出现这种情况
							lastSource_idx++;
						}
					}
				}
				else if ((*vv_cit).idx() == lastSource_idx){//因为上一个点已经从pointSetS中移除了,所以如果找到的是last source point,index一定小于0
					nowBuffer++;
					if (nowBuffer == 1){
						last2next = true;
						if (nextPoint_index < 0&&pftsBuffers[0].size()>0){//如果再也找不到下一个点
							nowBuffer++;
						}
					}
				}
				else{

					int c_index = indexOf(connecteds, (*vv_cit).idx(), false);
					if (c_index >= 0){//如果找到了一个traget point
						cout << "nowbuffer:"<<nowBuffer;
						pftsBuffers[nowBuffer].push_back(connecteds[c_index]);//把
					}
				}
			}
			if (last2next){//如果从last->next
				cout << "last"<<nextPoint_index<<"->next"<<lastSource_idx;
				for (int buffer = 0; buffer < 3; buffer++){
					for (int i = 0; i < pftsBuffers[buffer].size(); i++){
						cout << "(" << pftsBuffers[buffer][i].sourceIdx_Small << "," << pftsBuffers[buffer][i].sourceIdx_Large << ") ";
					}
					cout << "|";
				}
				cout << "        ";
				for (int b = pftsBuffers[0].size() - 1; b >= 0; b--){
					newGroup.push_back(pftsBuffers[0][b]);
					removePFTS(pointSetS, pftsBuffers[0][b]);
				}
				for (int f = 0; f < pftsBuffers[1].size(); f++){
					newGroup.push_back(pftsBuffers[1][f]);
					removePFTS(pointSetS, pftsBuffers[1][f]);
				}
				for (int b = pftsBuffers[2].size() - 1; b >= 0; b--){
					newGroup.push_back(pftsBuffers[2][b]);
					removePFTS(pointSetS, pftsBuffers[2][b]);
				}
			}
			else{
				cout << "next" << nextPoint_index<<"->last"<<lastSource_idx;
				for (int buffer = 0; buffer < 3; buffer++){
					for (int i = 0; i < pftsBuffers[buffer].size(); i++){
						cout << "(" << pftsBuffers[buffer][i].sourceIdx_Small << "," << pftsBuffers[buffer][i].sourceIdx_Large << ") ";
					}
					cout << "|";
				}
				cout << "        ";
				if (pointSource_idx == 2764){
					cout << "找到2764 下一个点是:" << nextPoint_index;
				}
				for (int f = 0; f < pftsBuffers[2].size(); f++){
					newGroup.push_back(pftsBuffers[2][f]);
					removePFTS(pointSetS, pftsBuffers[2][f]);
				}
				for (int b = pftsBuffers[1].size() - 1; b >= 0;b--){
					newGroup.push_back(pftsBuffers[1][b]);
					removePFTS(pointSetS, pftsBuffers[1][b]);
				}
				for (int f = 0; f < pftsBuffers[0].size(); f++){
					newGroup.push_back(pftsBuffers[0][f]);
					removePFTS(pointSetS, pftsBuffers[0][f]);
				}
			}
			if (firstPoint_idx == -1){
				firstPoint_idx = pointSource_idx;//记录第一个点的idx,也意味着解除第一次的标记
			}
			lastSource_idx = pointSource_idx;
			pointSource_idx = nextPoint_index;//设定pointSource_idx到下一个点的index
		}while (!getEndPointOfContour);// while (nextPoint_index >= 0);
		*/
		/*
		cout<<endl << "newGroup:";
		for each(pointFromTwoSource pfts in newGroup){
			//cout << "(" << pfts.pointPos[0] << "," << pfts.pointPos[1] << "," << pfts.pointPos[2] << ") ";
			cout << "idxs:(" << pfts.sourceIdx_Small << "," << pfts.sourceIdx_Large << ") ";
		}
		cout << endl;*/
		//cout << "-----------------------------------------------------------------" << endl;
		groupSmall.push_back(newGroup);
		contours->push_back(new contour(newGroup));
	}
	for each(vector<pointFromTwoSource> group in groupSmall){
		vector<OpenMesh::Vec3f> line;
		for each(pointFromTwoSource point in group){
			line.push_back(point.pointPos);
		}
		results->push_back(line);
	}
	return *results;
}
vector<vector<Vec3f_angle>> polarMap::dividePointSet(MyMesh mesh, vector<pointFromTwoSource> pointSet,vector<int> axis){
	vector<vector<pointFromTwoSource>> groupSmall;
	vector<vector<Vec3f_angle>> *results = new vector<vector<Vec3f_angle>>();
	vector<pointFromTwoSource> pointSetS(pointSet);//make a pointSet copy
	int counter = 0;
	while (pointSetS.size() > 0){
		vector<pointFromTwoSource> newGroup;
		pointFromTwoSource nowpfts = pointSetS[0];
		int sourceIdx = nowpfts.sourceIdx_Small;
		MyMesh::VertexHandle vh(sourceIdx);
		MyMesh::HalfedgeHandle nowHalfEdge;
		for (MyMesh::VertexIHalfedgeCWIter vhe_it = mesh.vih_cwbegin(vh); vhe_it != mesh.vih_cwend(vh); vhe_it++){
			//cout << "from vertex" << mesh.from_vertex_handle(*vhe_it).idx();
			if (mesh.from_vertex_handle(*vhe_it).idx() == nowpfts.sourceIdx_Large){
				//cout << " get!!!";
				nowHalfEdge = *vhe_it;
				break;
			}
		}
		bool findNext = false;
		MyMesh::HalfedgeHandle firstHalfEdge = nowHalfEdge;//記錄起始半邊,用於遇到邊界的反向重組
		bool firstMeetBoundary = true;
		//cout << endl << "newGroup add:";

		do{
			findNext = false;
			newGroup.push_back(nowpfts);
			//if (mesh.is_boundary(mesh.edge_handle(nowHalfEdge))){//如果遇到了邊界,即過了這條邊不再有面存在
			if (mesh.is_boundary(mesh.opposite_halfedge_handle(nowHalfEdge))){
				//cout<<"找到boundary!" << endl;
				if (firstMeetBoundary){//如果是第一次遇到
					firstMeetBoundary = false;
					nowHalfEdge = mesh.opposite_halfedge_handle(firstHalfEdge);//設置當前半邊為起始半邊的反向,這樣會向一開始的反方向重組pfts
				}
				else{//遇到第二次代表重組完成了,一般來說只要有第一次就一定會有第二次
					break;
				}
			}
			//cout << "(" << nowpfts.sourceIdx_Small << "," << nowpfts.sourceIdx_Large << ")";
			removePFTS(pointSetS, nowpfts);
			MyMesh::FaceHandle nowFace = mesh.opposite_face_handle(nowHalfEdge);//找到當前halfEdge的反向面
			for (MyMesh::FaceHalfedgeCWIter fhe_it = mesh.fh_cwbegin(nowFace); fhe_it != mesh.fh_cwend(nowFace); fhe_it++){//找找這個面里有沒有另一半邊屬於已知的pointFromTwoSource
				int vidx1 = mesh.from_vertex_handle(*fhe_it).idx();
				int vidx2 = mesh.to_vertex_handle(*fhe_it).idx();
				int index1 = indexOf(pointSetS, vidx1, vidx2);//先試試from 放前面to 放後面能不能找到pointFormTwoSource
				if (index1 >= 0){//找到了
					findNext = true;
					nowpfts = pointSetS[index1];
					nowHalfEdge = *fhe_it;
					break;
				}
				else{//沒找到
					//試試看反過來找
					int temp = vidx1;
					int index2 = indexOf(pointSetS, vidx2, vidx1);
					if (index2 >= 0){
						findNext = true;
						nowpfts = pointSetS[index2];
						nowHalfEdge = *fhe_it;
						break;
					}
				}
			}


		} while (findNext);
		groupSmall.push_back(newGroup);
		contours->push_back(new contour(newGroup));
	}
	//cout << ">>>>>>>>>>>>>>>counter:" << counter<<endl;
	for each(vector<pointFromTwoSource> group in groupSmall){//對於每個等高線
		//開始找等高線起始點
		int startIdx = -1;
		for (int i = 0; i < axis.size()-1; i++){
			for (int idx = 0; idx < group.size(); idx++){
				pointFromTwoSource pfts = group[idx];
				if ((pfts.sourceIdx_Small == axis[i] && pfts.sourceIdx_Large == axis[i + 1]) || (pfts.sourceIdx_Small == axis[i + 1] && pfts.sourceIdx_Large == axis[i])){
					startIdx = idx;
				}
			}
		}
		//從起點開始重組等高線
		if (startIdx == -1){//找不到起始點就直接GG
			//cout << "錯誤!找不到等高線的起點" << endl;
		}
		else{//按照環裝重組交點.
			//cout << "重組start:";
			vector<pointFromTwoSource> newGroup = vector<pointFromTwoSource>();
			for (int i = startIdx; i < group.size(); i++){
				newGroup.push_back(group[i]);
			}
			for (int i = 0; i < startIdx; i++){
				newGroup.push_back(group[i]);
			}
			//重組完成
			//計算等高線總長度
			float totalLength = 0;
			for (int i = 0; i < newGroup.size()-1; i++){
				OpenMesh::Vec3f toNext = newGroup[i + 1].pointPos - newGroup[i].pointPos;
				totalLength+= toNext.length();
			}
			totalLength +=(newGroup[0].pointPos - newGroup.back().pointPos).length();
			//cout << "totalLength:" << totalLength << "; ";
			//计算完成总长度
			vector<Vec3f_angle> *angleGroup =new vector<Vec3f_angle>();
			float nowpercent = 0;
			angleGroup->push_back(Vec3f_angle(newGroup[0].pointPos,nowpercent));//第一點一定是axis point 角度一定是0
			for (int i = 0; i < newGroup.size()-1; i++){
				nowpercent += ((newGroup[i+1].pointPos - newGroup[i].pointPos).length() / totalLength);
				//cout << "nowpercent:" << nowpercent << " ";
				angleGroup->push_back(Vec3f_angle(newGroup[i+1].pointPos, nowpercent));
			}
			results->push_back(*angleGroup);
			//cout << endl;
		}
	}
	counter++;
	return *results;
}
vector<PointCurve> polarMap::getContourLine(int line_number,vector<int> axis){
	cout << "OpenMesh::Vec3f size:" << sizeof(OpenMesh::Vec3f) << endl;
	cout << "记时开始" << endl;
	clock_t timeStart = clock();
	vector<PointCurve> curves;
	float interval = U[farthestIdx] / (line_number + 1);
	for (int i = 1; i <= line_number; i++){
		//PointCurve nowcurve;
		vector<pointFromTwoSource> pointSet;
		float nowLen = i*interval;
		//cout << "第" << i << "圈,nowLen為:" << nowLen<<endl;
		for (int v = 0; v < mesh.n_vertices();v++){
			if (U[v] == nowLen){
				//nowcurve.Add(Vec3f_angle(mesh.point(OpenMesh::VertexHandle(v)),theta[v]));
				pointSet.push_back(pointFromTwoSource(v,v,mesh.point(OpenMesh::VertexHandle(v))));
			}
			if (U[v] < nowLen){
				OpenMesh::VertexHandle vh(v);
				//cout << "nowLen:" << nowLen << " v:" << v << " ";
				for (MyMesh::VertexVertexCWIter vvc_it = mesh.vv_cwbegin(vh); vvc_it != mesh.vv_cwend(vh); vvc_it++){
					int v2 = (*vvc_it).idx();
					//cout << " v2:" << v2;
					if (U[v2] > nowLen){
						OpenMesh::Vec3f dir = mesh.point(*vvc_it)-mesh.point(vh);
						OpenMesh::Vec3f intersection = mesh.point(vh) + dir*((nowLen-U[v])/(U[v2]-U[v]));
						//nowcurve.Add(Vec3f_angle(intersection,theta[v]));
						//cout << "intersection:(" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ") ";
						pointSet.push_back(pointFromTwoSource(v, v2, intersection));
						//cout << "nowcurve add v:" << v << "theta:[v]" << theta[v] << endl;
					}
				}
			}
		}
		//cout << "pointSet length:" << pointSet.size();
		//for each(pointFromTwoSource pfts in pointSet){
			//cout << "point:" << pfts.pointPos[0] << "," << pfts.pointPos[1] << "," << pfts.pointPos[2] << " ";
		//}
		//cout << "-------------------------------------------------------"<<endl;
		//vector<vector<OpenMesh::Vec3f>> divided = dividePointSet(mesh, pointSet);
		if (axis.size() <= 0){
			vector<vector<OpenMesh::Vec3f>> divided = dividePointSet(mesh, pointSet);
			for each(vector<OpenMesh::Vec3f>points in divided){
				PointCurve newCurve;
				for each(OpenMesh::Vec3f point in points){
					newCurve.Add(Vec3f_angle(point, 0));//先不管角度
				}
				curves.push_back(newCurve);
			}
		}
		else{
			vector<vector<Vec3f_angle>> divided = dividePointSet(mesh, pointSet, axis);
			/*for each(vector<OpenMesh::Vec3f>points in divided){
				PointCurve newCurve;
				for each(OpenMesh::Vec3f point in points){
				newCurve.Add(Vec3f_angle(point, 0));//先不管角度
				}
				curves.push_back(newCurve);
				}*/
			for each(vector<Vec3f_angle> points in divided){
				PointCurve newCurve;
				for each(Vec3f_angle point in points){
					newCurve.Add(point);
				}
				curves.push_back(newCurve);
			}
		}
		
	}
	//divide Contour
	clock_t timeEnd = clock();
	double sec = (timeEnd - timeStart) / (double)(CLOCKS_PER_SEC);
	cout << "總花費時間:" << sec << endl;
	//cout << "sizeof curves is:" + sizeof(curves) << endl;
	return curves;
}
void polarMap::tryMeetNextContour(MyMesh::VertexHandle now_vh, contour current, vector<contour*> candiate, set<int> &meetContourIndex, set<int> &vertexHistory){

	vertexHistory.insert(now_vh.idx());
	for (MyMesh::VertexVertexCWIter vv_cwit = mesh.vv_cwbegin(now_vh); vv_cwit != mesh.vv_cwend(now_vh); vv_cwit++){
		if (vertexHistory.find((*vv_cwit).idx())!=vertexHistory.end())//*vv_cwit 已經在vertexHistory裡面了,就什麼都不做
		{
			continue;
		}
		if (current.contact_asSmall ((*vv_cwit).idx())){//如果接觸到當前等高線
			//啥都不做
			continue;
		}
		bool meetNext = false;
		for (int i = 0; i < candiate.size(); i++){
			if (candiate[i]->contact_asLarge((*vv_cwit).idx()))//如果跨過了下一條等高線
			{
				meetContourIndex.insert(i);//把索引值記錄在meetContourIndex裡面
				meetNext = true;
				break;
			}
		}
		if(!meetNext)
			tryMeetNextContour((*vv_cwit), current, candiate, meetContourIndex,vertexHistory);
	}
}
void polarMap::showContours(){
	cout << "= = =contours address:" << contours << endl;
	for (int i = 0; i < contours->size(); i++){
		cout << "第" << i << "條contour id:" << contours->operator[](i)->id << "address:" << contours->operator[](i) << " next:";
		for each (contour *c in  contours->operator[](i)->nextContours)
		{
			cout << c->id << ":"<<c<<"  ";
		}
		cout << endl;
	}
	//if (contours->operator[](6)->distance == contours->operator[](7)->distance){
	//	cout << contours->operator[](6)->distance << "==" << contours->operator[](7)->distance<<endl;
	//}

}
vector<PointCurve> polarMap::getContourLine(int line_number){
	//cout << "OpenMesh::Vec3f size:" << sizeof(OpenMesh::Vec3f) << endl;
	cout << "记时开始" << endl;
	clock_t timeStart = clock();
	vector<PointCurve> curves;
	totalLevel = line_number;
	intervalBtwLevels = U[farthestIdx] / (line_number + 1);
	pathLengthThreshold = 3 * intervalBtwLevels;
	cout << "設置pathLengthThreshold為:" << pathLengthThreshold << endl;
	vector<vector<contour*>> beforeSeries(line_number);
	int contourCounter = 0;
	for (int i = 1; i <= line_number; i++){
		vector<contour*> parallelGroups;
		float nowLen = i*intervalBtwLevels;
		vector<pointFromTwoSource> beforeDivided;
		//cout << "第" << i << "圈,nowLen為:" << nowLen<<endl;
		for (int v = 0; v < mesh.n_vertices(); v++){

			if (U[v] == nowLen){
				//nowcurve.Add(Vec3f_angle(mesh.point(OpenMesh::VertexHandle(v)),theta[v]));
				beforeDivided.push_back(pointFromTwoSource(v, v, mesh.point(OpenMesh::VertexHandle(v))));
			}
			if (U[v] < nowLen){
				OpenMesh::VertexHandle vh(v);
				//cout << "nowLen:" << nowLen << " v:" << v << " ";
				for (MyMesh::VertexVertexCWIter vvc_it = mesh.vv_cwbegin(vh); vvc_it != mesh.vv_cwend(vh); vvc_it++){
					int v2 = (*vvc_it).idx();
					//cout << " v2:" << v2;
					if (U[v2] > nowLen){
						OpenMesh::Vec3f dir = mesh.point(*vvc_it) - mesh.point(vh);
						OpenMesh::Vec3f intersection = mesh.point(vh) + dir*((nowLen - U[v]) / (U[v2] - U[v]));
						beforeDivided.push_back(pointFromTwoSource(v, v2, intersection));
					}
				}
			}
			
		}
		//pointSet.push_back(contour(beforeDivided));
		while (beforeDivided.size() > 0){//合併原dividePointSet在這裡
			vector<pointFromTwoSource> newGroup;
			pointFromTwoSource nowpfts = beforeDivided[0];
			int sourceIdx = nowpfts.sourceIdx_Small;
			MyMesh::VertexHandle vh(sourceIdx);
			MyMesh::HalfedgeHandle nowHalfEdge;
			for (MyMesh::VertexIHalfedgeCWIter vhe_it = mesh.vih_cwbegin(vh); vhe_it != mesh.vih_cwend(vh); vhe_it++){
				//cout << "from vertex" << mesh.from_vertex_handle(*vhe_it).idx();
				if (mesh.from_vertex_handle(*vhe_it).idx() == nowpfts.sourceIdx_Large){
					//cout << " get!!!";
					nowHalfEdge = *vhe_it;
					break;
				}
			}
			bool findNext = false;
			MyMesh::HalfedgeHandle firstHalfEdge = nowHalfEdge;//記錄起始半邊,用於遇到邊界的反向重組
			bool firstMeetBoundary = true;
			//cout << endl << "newGroup add:";

			do{
				findNext = false;
				newGroup.push_back(nowpfts);
				//if (mesh.is_boundary(mesh.edge_handle(nowHalfEdge))){//如果遇到了邊界,即過了這條邊不再有面存在
				if (mesh.is_boundary(mesh.opposite_halfedge_handle(nowHalfEdge))){
					//cout<<"找到boundary!" << endl;
					if (firstMeetBoundary){//如果是第一次遇到
						firstMeetBoundary = false;
						nowHalfEdge = mesh.opposite_halfedge_handle(firstHalfEdge);//設置當前半邊為起始半邊的反向,這樣會向一開始的反方向重組pfts
					}
					else{//遇到第二次代表重組完成了,一般來說只要有第一次就一定會有第二次
						break;
					}
				}
				//cout << "(" << nowpfts.sourceIdx_Small << "," << nowpfts.sourceIdx_Large << ")";
				removePFTS(beforeDivided, nowpfts);
				MyMesh::FaceHandle nowFace = mesh.opposite_face_handle(nowHalfEdge);//找到當前halfEdge的反向面
				for (MyMesh::FaceHalfedgeCWIter fhe_it = mesh.fh_cwbegin(nowFace); fhe_it != mesh.fh_cwend(nowFace); fhe_it++){//找找這個面里有沒有另一半邊屬於已知的pointFromTwoSource
					int vidx1 = mesh.from_vertex_handle(*fhe_it).idx();
					int vidx2 = mesh.to_vertex_handle(*fhe_it).idx();
					int index1 = indexOf(beforeDivided, vidx1, vidx2);//先試試from 放前面to 放後面能不能找到pointFormTwoSource
					if (index1 >= 0){//找到了
						findNext = true;
						nowpfts = beforeDivided[index1];
						nowHalfEdge = *fhe_it;
						break;
					}
					else{//沒找到
						//試試看反過來找
						int temp = vidx1;
						int index2 = indexOf(beforeDivided, vidx2, vidx1);
						if (index2 >= 0){
							findNext = true;
							nowpfts = beforeDivided[index2];
							nowHalfEdge = *fhe_it;
							break;
						}
					}
				}


			} while (findNext);
			contour *newContour=new contour(newGroup, nowLen);
			if (newContour->length>dropContourThreshold)
			{
				newContour->id = contourCounter++;
				parallelGroups.push_back(newContour);
			}
			else{

				cout << "drop newContour in contourCounter=" << contourCounter << endl;

			}
		}
		//cout << "line" << i << " ";
		beforeSeries[i-1]= parallelGroups;
	}//for line_number結束
	//series contour:指定contour的下一條contour是哪個
	contours = new vector<contour*>();
	//cout << "=====contours address:" << contours<<endl;
	for (int i = 0; i < beforeSeries.size()-1;i++){
		//cout << "series contour處理i=" << i << endl;
		if (beforeSeries[i].size()==0)//不合理的情況,當前有0條contour
		{
			cout << "錯誤!切割后第" << i << "條等高線數量為0!!!";
		}
		if (beforeSeries[i].size() == 1){//只有1條contour則i+1的所有contour都是下一條contour
			contour *now = beforeSeries[i][0];
			for (size_t j = 0; j < beforeSeries[i + 1].size(); j++)
			{
				now->nextContours.push_back(beforeSeries[i + 1][j]);
				beforeSeries[i + 1][j]->lastContours.push_back(now);
			}
			contours->push_back(now);
			now->calDivergence();
		}
		if (beforeSeries[i].size()>1)//有複數條contour時通過測試contour之間的連接性來分別找到每一條contour的下一條或幾條contour
		{

			for each (contour* c in beforeSeries[i])
			{
				MyMesh::VertexHandle startVertex = MyMesh::VertexHandle(c->pfts[0].sourceIdx_Large);
				set<int> result_index;
				set<int> allVertex;
				tryMeetNextContour(startVertex,*c,beforeSeries[i+1],result_index,allVertex);
				for each (int index in result_index)
				{
					c->nextContours.push_back(beforeSeries[i + 1][index]);
					beforeSeries[i + 1][index]->lastContours.push_back(c);
				}
				//cout << endl;
				contours->push_back(c);
				c->calDivergence();
			}
		}
	}
	for each (contour *c in beforeSeries.back())//最後一圈
	{
		contours->push_back(c);
	}
	initPathFinder();
	for each ( contour* c in *contours)
	{
		curves.push_back(c->toCurve());
	}
	removeIncorrectExtremumVertexs();
	clock_t timeEnd = clock();
	double sec = (timeEnd - timeStart) / (double)(CLOCKS_PER_SEC);
	cout << "getContourLine結束,總花費時間:" << sec << endl;
	//cout << "sizeof curves is:" + sizeof(curves) << endl;
	return curves;
}
float loopFindNearestCycle(float traget, float testAngle){
	int layer = floorf(traget);
	float nearestAngle =loopTransToLayer(testAngle,layer);
	if (abs(traget - loopTransToLayer(testAngle, layer - 1))<abs(traget-nearestAngle)){
		nearestAngle = loopTransToLayer(testAngle, layer - 1);
	}
	if (abs(traget - loopTransToLayer(testAngle, layer + 1)) < abs(traget - nearestAngle)){
		nearestAngle = loopTransToLayer(testAngle, layer + 1);
	}
	return nearestAngle;
}
float loopFindNearestBiggerValue(float traget, float testAngle){
	//cout << "biggerValue traget:" << traget << "testAngle:" << testAngle;
	int layer = floorf(traget);
	float value1 = loopTransToLayer(testAngle, layer);
	if (value1 > traget)
		return value1;
	float value2 = loopTransToLayer(testAngle, layer+1);
	if (value2 > traget)
		return value2;
	return value2;
}
float loopFindNearestSmallerValue(float traget, float testAngle){
	//cout << "biggerValue traget:" << traget << "testAngle:" << testAngle;
	int layer = floorf(traget);
	float value1 = loopTransToLayer(testAngle, layer);
	if (value1< traget)
		return value1;
	float value2 = loopTransToLayer(testAngle, layer -1);
	if (value2 < traget)
		return value2;
	return value2;
}
float loopTransToLayer(float value, int layer){
	float layer_value = floor(value);
	return value- (layer_value - (float)layer);
}
float loopClampTo01(float oriAngle){
	float result = oriAngle - (float)((int)oriAngle);
	if (result < 0){//判斷是不是負值且轉為正值,比如-3.3=>result = -3.3-(-3) =-0.3 
		result += 1;//-0.3+1 =0.7
	}
	return result;
}
float loopClampTo01(int index, int arraySize){
	if (index >= 0){
		return index%arraySize;
	}
	else{
		return arraySize + index;
	}
}
float loopInterpolateFromShort(float angle1, float angle2, float t){

}
vector<OpenMesh::Vec3f> polarMap::getPointBetween(int contour1_index,float at_angle,float percentage){
	vector<OpenMesh::Vec3f> points;
	contour contour1 = *contours->operator[](contour1_index);
	if (contour1.nextContours.size() <= 0){
		cout << "沒有下一個contour ";
		return points;
	}
	contour contour2 = *contour1.nextContours[0];
	pointFromTwoSource* pair1 = contour1.searchAngle(at_angle);
	//cout << "search angle之後:pair1[0] small:" << pair1[0].sourceIdx_Small << "large:" << pair1[0].sourceIdx_Large;
	OpenMesh::VertexHandle sv = OpenMesh::VertexHandle(pair1[0].sourceIdx_Small);
	OpenMesh::Vec3f nowPos = pointFromTwoSource::Interpolation(pair1[0], pair1[1], at_angle);//內插得到起始位置
	points.push_back(nowPos);
	//cout << "at_angle:" << at_angle << endl;
	OpenMesh::Vec3f endPos = contour2.getPosAtAngle(at_angle);
	glm::vec3 endPos_3d(endPos[0], endPos[1], endPos[2]);
	glm::vec3 oriPos_3d(nowPos[0], nowPos[1], nowPos[2]);
	MyMesh::FaceHandle nowFace;
	int updateFromEvidxp[2] = { -1, -1 };//nowPos衛位於的邊兩頂點的index
	OpenMesh::Vec3f mainVec=endPos - nowPos;
	cout << "now point:" << nowPos[0] << "," << nowPos[1] << "," << nowPos[2] << endl;
	cout << "起始點位於面:" << contour1.getFaceAtAngle(mesh, at_angle);
	cout << "end point:" << endPos[0] << "," << endPos[1] << "," << endPos[2] << endl;
	cout << "終點位於面:" << contour2.getFaceAtAngle(mesh, at_angle);
	for (MyMesh::VertexFaceIter vf_it = mesh.vf_begin(sv); vf_it != mesh.cvf_end(sv);vf_it++){//搜尋點所屬的所有面
		bool meetEdge1=false;
		bool meetEdge2 = false;
		//搜尋一個面,包含了所有pair1的所有點
		for (MyMesh::FaceHalfedgeIter fhe_it = mesh.fh_begin(*vf_it); fhe_it != mesh.fh_end(*vf_it); fhe_it++){
			MyMesh::HalfedgeHandle heh=*fhe_it;
			if ((mesh.from_vertex_handle(heh).idx() == pair1[0].sourceIdx_Small &&mesh.to_vertex_handle(heh).idx() == pair1[0].sourceIdx_Large) || (mesh.from_vertex_handle(heh).idx() == pair1[0].sourceIdx_Large &&mesh.to_vertex_handle(heh).idx() == pair1[0].sourceIdx_Small))
				meetEdge1 = true;
			else if ((mesh.from_vertex_handle(heh).idx() == pair1[1].sourceIdx_Small &&mesh.to_vertex_handle(heh).idx() == pair1[1].sourceIdx_Large) || (mesh.from_vertex_handle(heh).idx() == pair1[1].sourceIdx_Large &&mesh.to_vertex_handle(heh).idx() == pair1[1].sourceIdx_Small))
				meetEdge2 = true;
		}
		if (meetEdge1&&meetEdge2){
			nowFace = *vf_it;
			break;
		}
	}
	delete pair1;
	bool meetEndPoint = false;
	bool foundNext = false;
	do
	{
		foundNext = false;
		cout << "處理面:" << nowFace.idx()<<endl;
		glm::vec3 point3d[3];
		glm::vec2 point2d[3];
		int edge_vidx[3];
		int count = 0;
		for (MyMesh::FaceVertexCWIter fv_cwit = mesh.fv_cwbegin(nowFace); fv_cwit != mesh.fv_cwend(nowFace)&&count<3; fv_cwit++,count++)
		{
			OpenMesh::Vec3f pos = mesh.point(*fv_cwit);
			point3d[count] = glm::vec3(pos[0], pos[1], pos[2]);
			edge_vidx[count] = (*fv_cwit).idx();
		}
		glm::vec3 arrow1 = point3d[1] - point3d[0];
		glm::vec3 arrow2 = point3d[2] - point3d[0];
		/*cout << "point3d:" << endl;
		for (int i = 0; i < 3; i++){
			cout<<"idx:"<<edge_vidx[i]<< "(" << point3d[i][0] << "," << point3d[i][1] << "," << point3d[i][2] << ")" << endl;
		}*/
		
		float angle = acosf(glm::dot(arrow1, arrow2) / (glm::length(arrow1)*glm::length(arrow2)));
		//cout << "arrow1:(" << arrow1[0] << "," << arrow1[1] << "," << arrow1[2] << ") length:"<<length(arrow1)<<" arrow2:(" << arrow2[0] << "," << arrow2[1] << "," << arrow2[2] << ") length:"<<length(arrow2)<<" angle:"<<angle<<endl;
		//cout << "sin(angle)=" << sinf(angle) << " ,cos(angle)=" << cosf(angle)<<endl;
		//2d坐標系
		point2d[0] = glm::vec2(0, 0);//以p0為原點
		point2d[1] = glm::vec2(glm::length(arrow1)*cosf(0), glm::length(arrow1)*sinf(0));//p0->p1為坐標系x軸
		point2d[2] = glm::vec2(glm::length(arrow2)*cosf(angle),glm::length(arrow2)*sinf(angle));
		adapter_3dto2d trans(point3d, point2d);
		/*
		cout << "adapter_3dto2d matrix:" << endl;
		for (int y = 0; y < 3; y++){
			cout << "    " << trans.matrix[y][0] << "," << trans.matrix[y][1] << endl;
		}
	
		cout << "point2d:" << endl;
		for (int i = 0; i < 3; i++){
			cout << "(" << point2d[i][0] << "," << point2d[i][1] << ")";
		}*/
		glm::vec2 point2d_recal[3];
		//cout <<endl<< "重新算出來的point2d:";
		for (size_t i = 0; i < 3; i++)
		{
			point2d_recal[i]= trans.to2dPoint(point3d[i]);
			//cout << "(" << point2d_recal[i][0] << "," << point2d_recal[i][1] << ")";
		}
		//cout << endl;
		OpenMesh::Vec3f n = mesh.normal(nowFace);//法向量
		OpenMesh::Vec3f u = mainVec;// endPos - nowPos;
		OpenMesh::Vec3f p = u - (OpenMesh::dot(u, n) / n.length())*(n / n.length());
		//cout << "n:" << n << " u:" << u << " p:" << p<<endl;
		//2d點和向量----------------------------------------------------
		glm::vec2 pvec_2d = trans.to2dPoint(glm::vec3(p[0],p[1],p[2]));
		//cout << "pvec_2d:(" << pvec_2d[0]<<","<<pvec_2d[1]<<")"<<endl;
		float K_p = pvec_2d[1] / pvec_2d[0];//射線斜率 kp =y/x =>tan = sin/cos
		glm::vec3 nowPos_3d = glm::vec3(nowPos[0], nowPos[1], nowPos[2]);
		glm::vec2 pnow_2d = trans.to2dPoint(nowPos_3d);//轉化nowpos到2d空間
		//cout << "nowPos_2d(" << pnow_2d[0]<<","<<pnow_2d[1]<<")"<<endl;
		bool inEdge[3] {false,false,false};
		float timeRecord[3];
		
		for (int i = 0; i < 3; i++){//對於每一條邊測試一次有沒有相交,總共測試3次
			//cout << "第" << i << "次測試:" << endl;
			if ((updateFromEvidxp[0] == edge_vidx[(i + 1) % 3] && updateFromEvidxp[1] == edge_vidx[(i + 2) % 3]) || (updateFromEvidxp[1] == edge_vidx[(i + 1) % 3] && updateFromEvidxp[0] == edge_vidx[(i + 2) % 3])){
				//cout << "滿足邊條件,測試跳過";
				continue;
			}

			glm::vec2 edgePoint_2d[2];
			edgePoint_2d[0] = point2d[(i + 1) % 3];//i以外的兩個點作為邊
			edgePoint_2d[1] = point2d[(i + 2) % 3];



			glm::vec2 edgeVec_2d = edgePoint_2d[1] - edgePoint_2d[0];
			//cout << "edgeVec_2d:(" << edgeVec_2d[0] << "," << edgeVec_2d[1] << ") edgePoint_2d[0]:(" << edgePoint_2d[0][0] << "," << edgePoint_2d[0][1] << ") edgePoint_2d[1]:(" << edgePoint_2d[1][0] << "," << edgePoint_2d[1][1] << ")" << endl;
			float K_e = edgeVec_2d[1] / edgeVec_2d[0];//邊的斜率
			glm::mat2x3 problem(K_e, -1, K_e*edgePoint_2d[0][0] - edgePoint_2d[0][1],
								K_p, -1, K_p*pnow_2d[0]-pnow_2d[1]);
			/*cout << "problem 矩陣:" << endl;
			for (int y = 0; y < 2; y++){
				cout << "	" << problem[y][0] << " " << problem[y][1] << " " << problem[y][2] << endl;
			}*/
			if (problem[0][0] == 0){
				if (problem[1][0] == 0){
					cout << "錯誤!K_e與K_p同時等於0!將會觸發除0錯誤" << endl;
				}
				else//如果第二行滿足第一個不為0則交換行
				{
					glm::vec3 temp = problem[0];
					problem[0] = problem[1];
					problem[1] = temp;
				}
			}
			//解二項式
			problem[0] /= problem[0][0];
			problem[1] -= problem[1][0] * problem[0];
			problem[1] /= problem[1][1];
			problem[0] -= problem[0][1] * problem[1];
			/*cout << "解完后problem矩陣:" << endl;
			for (int y = 0; y < 2; y++){
				cout << "	" << problem[y][0] << " " << problem[y][1] << " " << problem[y][2] << endl;
			}*/

			glm::vec2 interaction(problem[0][2],problem[1][2]);
			//cout << "pvec_2d為:(" << pvec_2d[0] << "," << pvec_2d[1] << ")";
			//cout << "2d空間交點為:(" << interaction[0] << "," << interaction[1] << ")" << endl;
			
			glm::vec2 offset = interaction - pnow_2d;
			if (offset[0] / pvec_2d[0] < 0){//進入這裡說明是向射線的負方向
				//cout << "射線反方向!!!" << endl;
				continue;//忽略這個迴圈
			}
			glm::vec2 offset_e = interaction - edgePoint_2d[0];
			float time = offset_e[0] /edgeVec_2d[0];
			//cout << "time為:" << time<<" "<<endl;
			if (time < 0 || time>1){//超出邊線段的範圍
				//cout << "超出邊界範圍time為:" << time<<"offset_e:("<<offset_e[0]<<","<<offset_e[1]<<") edgeVec:("<<edgeVec_2d[0]<<","<<edgeVec_2d[1]<<")";
				inEdge[i] = true;
				timeRecord[i] = time;
				continue;//忽略這個迴圈
			}
			glm::vec3 interaction_3d = point3d[(i+1)%3] + (point3d[(i + 2) % 3] - point3d[(i + 1) % 3])*time;
			//cout << "nowPos:(" << nowPos[0] << "," << nowPos[1] << "," << nowPos[2] << ") 交點為:(" << interaction_3d[0] << "," << interaction_3d[1] << ","<<interaction_3d[2] << ")" << endl;
			//cout << "到交點距離:" << glm::length(interaction_3d - nowPos_3d) << " 到終點距離:" << glm::length(endPos_3d - nowPos_3d) << endl;
			if (glm::length(interaction_3d - oriPos_3d)>=glm::length(endPos_3d-oriPos_3d)){
			//if (glm::length(interaction_3d - nowPos_3d) >= glm::length(endPos_3d - nowPos_3d)){//交點到當前點的距離比終點到當前點的距離更遠,說明已經超過了endPos
				//cout << "meetEndPoint!" << endl;
				meetEndPoint = true;
				break;
			}
			else{//切到下一個面繼續
				points.push_back(OpenMesh::Vec3f(interaction_3d[0],interaction_3d[1],interaction_3d[2]));
				nowPos = points.back();
				//cout << "找到點(" << nowPos[0] << "," << nowPos[1] << "," << nowPos[2]<<")"<<endl;
				//cout << "當前面idx:" << nowFace.idx() << " fh_cwbegin:" << mesh.fh_cwbegin(nowFace);
				for (MyMesh::FaceHalfedgeCWIter fh_cwit = mesh.fh_cwbegin(nowFace); fh_cwit != mesh.fh_cwend(nowFace); fh_cwit++){
					int from_idx = mesh.from_vertex_handle(*fh_cwit).idx();
					int to_idx = mesh.to_vertex_handle(*fh_cwit).idx();
					//cout << "form_idx:" << from_idx << " to_idx:" << to_idx;
					if ((from_idx == edge_vidx[(i + 1) % 3] && to_idx == edge_vidx[(i + 2) % 3]) || (to_idx == edge_vidx[(i + 1) % 3] && from_idx == edge_vidx[(i + 2) % 3]))//找到交點位於哪條半邊
					{
						nowFace = mesh.opposite_face_handle(*fh_cwit);//設置nowFace到半邊的反面,即下一個面
						updateFromEvidxp[0] = edge_vidx[(i + 1) % 3];
						updateFromEvidxp[1] = edge_vidx[(i + 2) % 3];
						//cout << "更新nowFace到" << nowFace.idx()<<endl;
						foundNext = true;
						break;
					}
				}
				break;
			}
			//到了這裡一定是因為3個測試都沒有成功,可以算完全不能處理的case了
		}
		/*
		if (!foundNext){
			//cout << endl << "沒有找到任何可以繼續的邊,嘗試搶救一下" << endl;
			int mostIndex=-1;
			float smallestOffset=-1;
			for (int i = 0; i < 3; i++){
				if (inEdge[i]){
					float offset = 0;//計算offset,offset為計算time時,交點
					if (timeRecord[i] < 0){
						offset = timeRecord[i];
					}
					else//offset位於0到1之間的case不存在,所以這裡offset一定大於1
					{
						offset = timeRecord[i] - 1;
					}

					if (mostIndex < 0){
						mostIndex = i;
						smallestOffset = offset;
					}
					else
					{
						if (fabsf(offset) < fabsf(smallestOffset)){//比已知offset更小
							mostIndex = i;
							smallestOffset = offset;
						}
					}
				}
			}
			if (mostIndex >= 0){//有候選人
				//搶救一下,將交點最接近的邊的端點視為交點
				glm::vec3 fake_interaction;
				if (smallestOffset < 0){
					fake_interaction = point3d[(mostIndex + 1) % 3];
				}
				else{
					fake_interaction = point3d[(mostIndex + 2) % 3];
				}
				
				points.push_back(OpenMesh::Vec3f(fake_interaction[0], fake_interaction[1], fake_interaction[2]));
				nowPos = points.back();
				//cout << "找到點(" << nowPos[0] << "," << nowPos[1] << "," << nowPos[2] << ")" << endl;
				//cout << "當前面idx:" << nowFace.idx() << " fh_cwbegin:" << mesh.fh_cwbegin(nowFace);
				for (MyMesh::FaceHalfedgeCWIter fh_cwit = mesh.fh_cwbegin(nowFace); fh_cwit != mesh.fh_cwend(nowFace); fh_cwit++){
					int from_idx = mesh.from_vertex_handle(*fh_cwit).idx();
					int to_idx = mesh.to_vertex_handle(*fh_cwit).idx();
					//cout << "form_idx:" << from_idx << " to_idx:" << to_idx;
					if ((from_idx == edge_vidx[(mostIndex + 1) % 3] && to_idx == edge_vidx[(mostIndex + 2) % 3]) || (to_idx == edge_vidx[(mostIndex + 1) % 3] && from_idx == edge_vidx[(mostIndex + 2) % 3]))//找到交點位於哪條半邊
					{
						nowFace = mesh.opposite_face_handle(*fh_cwit);//設置nowFace到半邊的反面,即下一個面
						updateFromEvidxp[0] = edge_vidx[(mostIndex + 1) % 3];
						updateFromEvidxp[1] = edge_vidx[(mostIndex + 2) % 3];
						//cout << "更新nowFace到" << nowFace.idx()<<endl;
						foundNext = true;
						break;
					}
				}
			}
			else
			{
				cout << "搶救不了了!"<<endl;
			}
		}*/

	} while (!meetEndPoint&&foundNext);
	//cout << "getPointBetween 結束" << endl;
	return points;
}
glm::vec3 toGlmVec3(OpenMesh::Vec3f v){
	return glm::vec3(v[0],v[1],v[2]);
}
OpenMesh::Vec3f toMeshVec3f(glm::vec3 v){
	return OpenMesh::Vec3f(v[0], v[1], v[2]);
}
vector<OpenMesh::Vec3f> polarMap::getAxisAt(int contour1_index, float at_angle,float vec_length){
	vector<OpenMesh::Vec3f> points;
	contour contour1 = *contours->operator[](contour1_index);
	pointFromTwoSource* pair1 = contour1.searchAngle(at_angle);
	glm::vec3 startPoint = toGlmVec3(pointFromTwoSource::Interpolation(pair1[0], pair1[1], at_angle));
	MyMesh::FaceHandle face_start= contour1.getFaceAtAngle(mesh,at_angle);
	glm::vec3 norm_start = glm::normalize(toGlmVec3(mesh.normal(face_start)));
	//glm::vec3 spair_vec = glm::normalize(toGlmVec3(pair1[1].pointPos - pair1[0].pointPos));
	if ((contour1).nextContours.size() <= 0)
	{
		cout << "contour" << contour1.id << "沒有下一個contour.";
		return points;
	}
	contour contour2 = *((contours->operator[](contour1_index))->nextContours[0]);
	pointFromTwoSource* pair2 = contour2.searchAngle(at_angle);
	glm::vec3 endPoint = toGlmVec3(pointFromTwoSource::Interpolation(pair2[0],pair2[1],at_angle));
	MyMesh::FaceHandle face_end = contour2.getFaceAtAngle(mesh,at_angle);
	glm::vec3 norm_end = glm::normalize(toGlmVec3(mesh.normal(face_end)));
	//glm::vec3 epair_vec = glm::normalize(toGlmVec3(pair2[1].pointPos - pair2[0].pointPos));
	/*if (glm::dot(spair_vec, epair_vec) < 0){//說明夾角大於90度
		epair_vec = -epair_vec;
	}*/

	glm::vec3 axis_x = glm::normalize(endPoint - startPoint);
	cout << "axis_x:(" << axis_x[0] << "," << axis_x[1] << "," << axis_x[2] << ")" << endl;
	glm::vec3 temp_y =glm::normalize(toGlmVec3(pair1[1].pointPos) - startPoint);
	cout << "temp_y:(" << temp_y[0] << "," << temp_y[1] << "," << temp_y[2] << ")" << endl;
	glm::vec3 axis_z = glm::normalize(glm::cross(axis_x, temp_y));
	cout << "axis_z:(" << axis_z[0] << "," << axis_z[1] << "," << axis_z[2] << ")" << endl;
	glm::vec3 axis_y = glm::normalize(glm::cross(axis_z, axis_x));
	cout << "axis_y:" << axis_y[0] << "," << axis_y[1] << "," << axis_y[2] << ")" << endl;
	glm::vec3 mix_z = glm::normalize(norm_start + norm_end);
	glm::vec3 mix_y = glm::normalize(glm::cross(axis_x,mix_z));
	glm::vec3 new_z = glm::normalize(glm::cross(mix_y, axis_x));
	OpenMesh::Vec3f xp1 = OpenMesh::Vec3f(startPoint[0],startPoint[1],startPoint[2]);
	glm::vec3 xVec_end = startPoint + axis_x*vec_length;
	OpenMesh::Vec3f xp2 = OpenMesh::Vec3f(xVec_end[0], xVec_end[1], xVec_end[2]);
	points.push_back(xp1);
	points.push_back(xp2);

	OpenMesh::Vec3f yp1 = OpenMesh::Vec3f(startPoint[0], startPoint[1], startPoint[2]);
	glm::vec3 yVec_end = startPoint + axis_y*vec_length;
	OpenMesh::Vec3f yp2 = OpenMesh::Vec3f(yVec_end[0], yVec_end[1], yVec_end[2]);
	//points.push_back(yp1);
	//points.push_back(yp2);
	
	glm::vec3 tyVec_end = startPoint + temp_y*vec_length;
	OpenMesh::Vec3f typ2 = OpenMesh::Vec3f(tyVec_end[0],tyVec_end[1],tyVec_end[2]);
	//points.push_back(yp1);
	//points.push_back(typ2);

	glm::vec3 zVec_end = startPoint + axis_z*vec_length;
	OpenMesh::Vec3f zp2 = OpenMesh::Vec3f(zVec_end[0], zVec_end[1], zVec_end[2]);
	//points.push_back(yp1);
	//points.push_back(zp2);
	glm::vec3 mixZVec_end = startPoint + mix_z*vec_length;
	OpenMesh::Vec3f mzp2 = OpenMesh::Vec3f(mixZVec_end[0], mixZVec_end[1], mixZVec_end[2]);
	//points.push_back(yp1);
	//points.push_back(mzp2);

	glm::vec3 mixYVec_end = startPoint + mix_y*vec_length;
	OpenMesh::Vec3f myp2 = OpenMesh::Vec3f(mixYVec_end[0], mixYVec_end[1], mixYVec_end[2]);
	points.push_back(yp1);
	points.push_back(myp2);

	glm::vec3 newZVec_end = startPoint + new_z*vec_length;
	OpenMesh::Vec3f nzp2 = OpenMesh::Vec3f(newZVec_end[0], newZVec_end[1], newZVec_end[2]);
	points.push_back(yp1);
	points.push_back(nzp2);
	return points;
}

vector<OpenMesh::Vec3f> polarMap::getPointsBtw(int contour1_index, float at_angle){
	const int max_iter_num = 100;
	vector<OpenMesh::Vec3f> points;
	contour contour1 = *contours->operator[](contour1_index);
	if ((contour1).nextContours.size() <= 0)
	{
		cout << "contour" << contour1.id << "沒有下一個contour.";
		return points;
	}
	pointFromTwoSource* pair1 = contour1.searchAngle(at_angle);
	glm::vec3 startPoint=toGlmVec3(pointFromTwoSource::Interpolation(pair1[0], pair1[1], at_angle));
	MyMesh::FaceHandle face_start = contour1.getFaceAtAngle(mesh, at_angle);
	glm::vec3 norm_start = glm::normalize(toGlmVec3(mesh.normal(face_start)));

	contour contour2 = *((contours->operator[](contour1_index))->nextContours[0]);
	glm::vec3 endPoint= toGlmVec3(contour2.getPosAtAngle(at_angle));
	MyMesh::FaceHandle face_end = contour2.getFaceAtAngle(mesh, at_angle);
	glm::vec3 norm_end = glm::normalize(toGlmVec3(mesh.normal(face_end)));

	int endFaceIdx = contour2.getFaceAtAngle(mesh, at_angle).idx();
	points.push_back(OpenMesh::Vec3f(startPoint[0], startPoint[1], startPoint[2]));
	//cout << "getPointsBtw " << contour1.id << " angle:" << at_angle << ":-------------------------------------------" << endl;
	/*
	glm::vec3 test_x =glm::normalize(glm::vec3(2.19, 0, 0));
	cout << "test_x:(" << test_x[0] << "," << test_x[1] << "," << test_x[2] << ")" << endl;
	glm::vec3 test_y = glm::normalize(glm::vec3(1,1,0));
	cout << "test_y:(" << test_y[0] << "," << test_y[1] << "," << test_y[2] << ")" << endl;
	glm::vec3 crossResult = glm::cross(test_x, test_y);
	cout << "cross result:(" << crossResult[0] << "," << crossResult[1] << "," << crossResult[2] << ")" << endl;
	glm::vec3 test_z = glm::normalize(crossResult);
	cout << "test_z:(" << test_z[0] << "," << test_z[1] << "," << test_z[2] << ")" << endl;
	test_y = glm::cross(test_z,test_x);
	cout << "test_y:(" << test_y[0] << "," << test_y[1] << "," << test_y[2] << ")" << endl;

	glm::vec3 testPoint(0.79, -11, 0.881);
	cout << "測試(0.79,-11,0.881):";
	float x = glm::dot(testPoint, test_x);
	float y = glm::dot(testPoint, test_y);
	float z = glm::dot(testPoint, test_z);
	cout << "(" << x << "," << y << "," << z << ")"<<endl;
	cout << "----------------------------------------------------------------------" << endl;*/
	glm::vec3 axis_x = glm::normalize(endPoint - startPoint);
	//cout << "axis_x:(" << axis_x[0] << "," << axis_x[1] << "," << axis_x[2] << ")" << endl;
	/*glm::vec3 temp_y = glm::normalize(toGlmVec3(pair1[1].pointPos)-startPoint);
	cout << "temp_y:(" << temp_y[0] << "," << temp_y[1] << "," << temp_y[2] << ")" << endl;
	glm::vec3 axis_z = glm::normalize(glm::cross(axis_x,temp_y));
	cout << "axis_z:(" << axis_z[0] << "," << axis_z[1] << "," << axis_z[2] << ")" << endl;
	glm::vec3 axis_y = glm::normalize(glm::cross(axis_z,axis_x));
	cout << "axis_y:" << axis_y[0] << "," << axis_y[1] << "," << axis_y[2] << ")" << endl;*/
	glm::vec3 temp_z = glm::normalize(norm_start+norm_end);
	//cout << "temp_z:(" << temp_z[0] << "," << temp_z[1] << "," << temp_z[2] << ")" << endl;
	glm::vec3 axis_y = glm::normalize(glm::cross(axis_x, temp_z));
	//cout << "axis_y:(" << axis_y[0] << "," << axis_y[1] << "," << axis_y[2] << ")" << endl;
	glm::vec3 axis_z = glm::normalize(glm::cross(axis_y, axis_x));
	//cout << "axis_z:(" << axis_z[0] << "," << axis_z[1] << "," << axis_z[2] << ")" << endl;

	coverAxis trans = coverAxis(startPoint, axis_x, axis_y, axis_z);
	MyMesh::FaceHandle nowFace = contour1.getFaceAtAngle(mesh, at_angle);
	//cout << "起始面為:" << nowFace.idx() << endl;
	glm::vec3 p1 = trans.cover(toGlmVec3(pair1[0].pointPos));
	glm::vec3 p2 = trans.cover(toGlmVec3(pair1[1].pointPos));
	//cout << "covered pair1[0].pointpos:(" << p1[0] << "," << p1[1] << "," << p1[2] << ")"<<endl;
	//cout << "covered pair1[1].pointpos:(" << p2[0] << "," << p2[1] << "," << p2[2] << ")" << endl;
	glm::vec3 pe = trans.cover(endPoint);
	//cout << "covered endPoint:(" << pe[0] << "," << pe[1] << "," << pe[2] << ")";

	MyMesh::HalfedgeHandle mostPossible_he;
	OpenMesh::Vec3f most_interaction;
	float  most_avgX=-9999;
	for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(nowFace); fh_it != mesh.fh_end(nowFace); fh_it++){
		MyMesh::VertexHandle from = mesh.from_vertex_handle(*fh_it);
		MyMesh::VertexHandle to = mesh.to_vertex_handle(*fh_it);
		glm::vec3 point1 = trans.cover(toGlmVec3(mesh.point(from)));
		glm::vec3 point2 = trans.cover(toGlmVec3(mesh.point(to)));
		//cout << "起始邊:"<<from.idx()<<" (" << point1[0] << "," << point1[1] << "," << point1[2] << "); "<<to.idx()<<"(" << point2[0] << "," << point2[1] << "," << point2[2] << ");" << endl;
		bool isInteraction = (point1[1] >= 0 && point2[1] <= 0) || (point1[1] <= 0 && point2[1] >= 0);//如果有相交的話兩個點會有一個端點y值為正另一個端點y值為負
		if (isInteraction){
			glm::vec3 offset = point2 - point1;
			float time = (0 - point1[1]) / offset[1];
			glm::vec3 coverDomain_intersection = point1 + time*offset;
			float avgX = coverDomain_intersection[0];
			if (avgX>most_avgX){
				mostPossible_he = *fh_it;
				most_avgX = avgX;
				//glm::vec3 P1toP2 = point2 - point1;
				//float time = (0 - point1[1]) / P1toP2[1];
				OpenMesh::Vec3f interaction = mesh.point(from) + (mesh.point(to) - mesh.point(from))*time;
				most_interaction = interaction;
			}
		}
		//points.push_back(OpenMesh::Vec3f(interaction[0], interaction[1], interaction[2]));
	}
	//cout << "搜尋點開始:" << endl;
	points.push_back(most_interaction);
	MyMesh::HalfedgeHandle now_he = mesh.opposite_halfedge_handle(mostPossible_he);
	nowFace = mesh.opposite_face_handle(mostPossible_he);
	//cout << "終點位於面:" << endFaceIdx << endl;
	//cout << "第一個面為:" << nowFace.idx()<<endl;
	int count = 0;
	while (nowFace.idx() != endFaceIdx){
		count++;
		if (count >= max_iter_num){//強制終止
			//cout << "第一百次了!!!";
			vector<OpenMesh::Vec3f> failResult(2);
			failResult[0] = toMeshVec3f(startPoint);
			failResult[1] = toMeshVec3f(endPoint);
			return failResult;
		}
		//cout << "處理面 " << nowFace.idx()<<":";
		for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(nowFace); fh_it != mesh.fh_end(nowFace); fh_it++){
			
			if ((*fh_it) != now_he){//忽略更新的來源邊
				//cout << "半邊:" << (*fh_it).idx() << "反向半邊為:" << mesh.opposite_halfedge_handle(*fh_it).idx() << "反向面為:" << mesh.opposite_face_handle(*fh_it);
				OpenMesh::Vec3f from = mesh.point(mesh.from_vertex_handle(*fh_it));
				OpenMesh::Vec3f to = mesh.point(mesh.to_vertex_handle(*fh_it));
				glm::vec3 point1 = trans.cover(toGlmVec3(from));
				//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
				glm::vec3 point2 = trans.cover(toGlmVec3(to));
				//cout << "point2 idx"<<mesh.to_vertex_handle(*fh_it).idx()<<":(" << point2[0] << "," << point2[1] << "," << point2[2] << ") ";
				if ((point1[1] >= 0 && point2[1] <= 0) || (point1[1] <= 0 && point2[1] >= 0)){//找到一個符合兩端的各在x軸一邊的邊
					//cout << "從半邊:" << (*fh_it).idx();
					glm::vec3 P1toP2 = point2 - point1;
					float time = (0 - point1[1]) / P1toP2[1];
					//glm::vec3 interaction = point1 + P1toP2*time;
					OpenMesh::Vec3f intersection = from + (to - from)*time;
					points.push_back(intersection);
					//cout << "找到點(" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")";
					
					now_he = mesh.opposite_halfedge_handle(*fh_it);//反向當前半邊
					//cout << "更新半邊到:" << now_he.idx() << endl;
					nowFace = mesh.face_handle(now_he);//推進面
					break;
				}
				//cout << ";     ";
			}
		}
		//cout << endl;
	}
	points.push_back(OpenMesh::Vec3f(endPoint[0], endPoint[1], endPoint[2]));
	//cout << "getPointsBtw end----------------------------------------------"<<endl;
	return points;
}
vector<OpenMesh::Vec3f> polarMap::getPointsBtw(int contour1_index, float at_angle,bool &fail){
	const int max_iter_num = 100;
	vector<OpenMesh::Vec3f> points;
	contour contour1 = *contours->operator[](contour1_index);
	if ((contour1).nextContours.size() <= 0)
	{
		cout << "contour" << contour1.id << "沒有下一個contour.";
		fail = false;
		return points;
	}
	pointFromTwoSource* pair1 = contour1.searchAngle(at_angle);
	glm::vec3 startPoint = toGlmVec3(pointFromTwoSource::Interpolation(pair1[0], pair1[1], at_angle));
	MyMesh::FaceHandle face_start = contour1.getFaceAtAngle(mesh, at_angle);
	glm::vec3 norm_start = glm::normalize(toGlmVec3(mesh.normal(face_start)));

	contour contour2 = *((contours->operator[](contour1_index))->nextContours[0]);
	glm::vec3 endPoint = toGlmVec3(contour2.getPosAtAngle(at_angle));
	MyMesh::FaceHandle face_end = contour2.getFaceAtAngle(mesh, at_angle);
	glm::vec3 norm_end = glm::normalize(toGlmVec3(mesh.normal(face_end)));

	int endFaceIdx = contour2.getFaceAtAngle(mesh, at_angle).idx();
	points.push_back(OpenMesh::Vec3f(startPoint[0], startPoint[1], startPoint[2]));
	//cout << "getPointsBtw " << contour1.id << " angle:" << at_angle << ":-------------------------------------------" << endl;
	/*
	glm::vec3 test_x =glm::normalize(glm::vec3(2.19, 0, 0));
	cout << "test_x:(" << test_x[0] << "," << test_x[1] << "," << test_x[2] << ")" << endl;
	glm::vec3 test_y = glm::normalize(glm::vec3(1,1,0));
	cout << "test_y:(" << test_y[0] << "," << test_y[1] << "," << test_y[2] << ")" << endl;
	glm::vec3 crossResult = glm::cross(test_x, test_y);
	cout << "cross result:(" << crossResult[0] << "," << crossResult[1] << "," << crossResult[2] << ")" << endl;
	glm::vec3 test_z = glm::normalize(crossResult);
	cout << "test_z:(" << test_z[0] << "," << test_z[1] << "," << test_z[2] << ")" << endl;
	test_y = glm::cross(test_z,test_x);
	cout << "test_y:(" << test_y[0] << "," << test_y[1] << "," << test_y[2] << ")" << endl;

	glm::vec3 testPoint(0.79, -11, 0.881);
	cout << "測試(0.79,-11,0.881):";
	float x = glm::dot(testPoint, test_x);
	float y = glm::dot(testPoint, test_y);
	float z = glm::dot(testPoint, test_z);
	cout << "(" << x << "," << y << "," << z << ")"<<endl;
	cout << "----------------------------------------------------------------------" << endl;*/
	glm::vec3 axis_x = glm::normalize(endPoint - startPoint);
	//cout << "axis_x:(" << axis_x[0] << "," << axis_x[1] << "," << axis_x[2] << ")" << endl;
	/*glm::vec3 temp_y = glm::normalize(toGlmVec3(pair1[1].pointPos)-startPoint);
	cout << "temp_y:(" << temp_y[0] << "," << temp_y[1] << "," << temp_y[2] << ")" << endl;
	glm::vec3 axis_z = glm::normalize(glm::cross(axis_x,temp_y));
	cout << "axis_z:(" << axis_z[0] << "," << axis_z[1] << "," << axis_z[2] << ")" << endl;
	glm::vec3 axis_y = glm::normalize(glm::cross(axis_z,axis_x));
	cout << "axis_y:" << axis_y[0] << "," << axis_y[1] << "," << axis_y[2] << ")" << endl;*/
	glm::vec3 temp_z = glm::normalize(norm_start + norm_end);
	//cout << "temp_z:(" << temp_z[0] << "," << temp_z[1] << "," << temp_z[2] << ")" << endl;
	glm::vec3 axis_y = glm::normalize(glm::cross(axis_x, temp_z));
	//cout << "axis_y:(" << axis_y[0] << "," << axis_y[1] << "," << axis_y[2] << ")" << endl;
	glm::vec3 axis_z = glm::normalize(glm::cross(axis_y, axis_x));
	//cout << "axis_z:(" << axis_z[0] << "," << axis_z[1] << "," << axis_z[2] << ")" << endl;

	coverAxis trans = coverAxis(startPoint, axis_x, axis_y, axis_z);
	MyMesh::FaceHandle nowFace = contour1.getFaceAtAngle(mesh, at_angle);
	//cout << "起始面為:" << nowFace.idx() << endl;
	glm::vec3 p1 = trans.cover(toGlmVec3(pair1[0].pointPos));
	glm::vec3 p2 = trans.cover(toGlmVec3(pair1[1].pointPos));
	//cout << "covered pair1[0].pointpos:(" << p1[0] << "," << p1[1] << "," << p1[2] << ")"<<endl;
	//cout << "covered pair1[1].pointpos:(" << p2[0] << "," << p2[1] << "," << p2[2] << ")" << endl;
	glm::vec3 pe = trans.cover(endPoint);
	//cout << "covered endPoint:(" << pe[0] << "," << pe[1] << "," << pe[2] << ")";

	MyMesh::HalfedgeHandle mostPossible_he;
	cout << "起始面Id:" << nowFace.idx();
	OpenMesh::Vec3f most_interaction;
	float  most_avgX = -9999;
	for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(nowFace); fh_it != mesh.fh_end(nowFace); fh_it++){
		MyMesh::VertexHandle from = mesh.from_vertex_handle(*fh_it);
		MyMesh::VertexHandle to = mesh.to_vertex_handle(*fh_it);
		glm::vec3 point1 = trans.cover(toGlmVec3(mesh.point(from)));
		glm::vec3 point2 = trans.cover(toGlmVec3(mesh.point(to)));
		//cout << "起始邊:"<<from.idx()<<" (" << point1[0] << "," << point1[1] << "," << point1[2] << "); "<<to.idx()<<"(" << point2[0] << "," << point2[1] << "," << point2[2] << ");" << endl;
		bool isInteraction = (point1[1] >= 0 && point2[1] <= 0) || (point1[1] <= 0 && point2[1] >= 0);//如果有相交的話兩個點會有一個端點y值為正另一個端點y值為負
		if (isInteraction){
			glm::vec3 offset = point2 - point1;
			float time = (0 - point1[1]) / offset[1];
			glm::vec3 coverDomain_intersection = point1 + time*offset;
			float avgX = coverDomain_intersection[0];
			if (avgX>most_avgX){
				mostPossible_he = *fh_it;
				most_avgX = avgX;
				//glm::vec3 P1toP2 = point2 - point1;
				//float time = (0 - point1[1]) / P1toP2[1];
				OpenMesh::Vec3f interaction = mesh.point(from) + (mesh.point(to) - mesh.point(from))*time;
				most_interaction = interaction;
				cout << "most_interaction v1:" << from.idx() << " v2:" << to.idx();
			}
		}
		//points.push_back(OpenMesh::Vec3f(interaction[0], interaction[1], interaction[2]));
	}
	//cout << "搜尋點開始:" << endl;
	points.push_back(most_interaction);
	MyMesh::HalfedgeHandle now_he = mesh.opposite_halfedge_handle(mostPossible_he);
	nowFace = mesh.opposite_face_handle(mostPossible_he);
	//cout << "終點位於面:" << endFaceIdx << endl;
	//cout << "第一個面為:" << nowFace.idx()<<endl;
	int count = 0;
	while (nowFace.idx() != endFaceIdx){
		count++;
		if (count >= max_iter_num){//強制終止
			//cout << "第一百次了!!!";
			fail = true;
			vector<OpenMesh::Vec3f> failResult(2);
			failResult[0] = toMeshVec3f(startPoint);
			failResult[1] = toMeshVec3f(endPoint);
			return failResult;
		}
		//cout << "處理面 " << nowFace.idx()<<":";
		for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(nowFace); fh_it != mesh.fh_end(nowFace); fh_it++){

			if ((*fh_it) != now_he){//忽略更新的來源邊
				//cout << "半邊:" << (*fh_it).idx() << "反向半邊為:" << mesh.opposite_halfedge_handle(*fh_it).idx() << "反向面為:" << mesh.opposite_face_handle(*fh_it);
				OpenMesh::Vec3f from = mesh.point(mesh.from_vertex_handle(*fh_it));
				OpenMesh::Vec3f to = mesh.point(mesh.to_vertex_handle(*fh_it));
				glm::vec3 point1 = trans.cover(toGlmVec3(from));
				//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
				glm::vec3 point2 = trans.cover(toGlmVec3(to));
				//cout << "point2 idx"<<mesh.to_vertex_handle(*fh_it).idx()<<":(" << point2[0] << "," << point2[1] << "," << point2[2] << ") ";
				if ((point1[1] >= 0 && point2[1] <= 0) || (point1[1] <= 0 && point2[1] >= 0)){//找到一個符合兩端的各在x軸一邊的邊
					cout << "interaction v1:" << mesh.from_vertex_handle(*fh_it) << " v2:" << mesh.to_vertex_handle(*fh_it);
					//cout << "從半邊:" << (*fh_it).idx();
					glm::vec3 P1toP2 = point2 - point1;
					float time = (0 - point1[1]) / P1toP2[1];
					//glm::vec3 interaction = point1 + P1toP2*time;
					OpenMesh::Vec3f intersection = from + (to - from)*time;
					points.push_back(intersection);
					//cout << "找到點(" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")";

					now_he = mesh.opposite_halfedge_handle(*fh_it);//反向當前半邊
					//cout << "更新半邊到:" << now_he.idx() << endl;
					nowFace = mesh.face_handle(now_he);//推進面
					break;
				}
				//cout << ";     ";
			}
		}
		//cout << endl;
	}
	points.push_back(OpenMesh::Vec3f(endPoint[0], endPoint[1], endPoint[2]));
	//cout << "getPointsBtw end----------------------------------------------"<<endl;
	fail = false;
	return points;
}
vector<OpenMesh::Vec3f> polarMap::getPointsBtw(int contour1_index, float at_angle1,float at_angle2, bool &fail){
	//cout << "測試-9.99 loopClampTo01:" << loopClampTo01(-9.99) << endl;
	//cout << "測試 100.0 loopClampTo01:" << loopClampTo01(1000.0) << endl;
	//cout << "測試- 1.478 loopClampTo01:" << loopClampTo01(-1.478) << endl;
	/*if (contour1_index == 1){
		cout << at_angle1 << ">>" << at_angle2 << ";";
	}*/
	at_angle1 = loopClampTo01(at_angle1);
	at_angle2 = loopClampTo01(at_angle2);
	//cout << "clamp to 01后 at_angle1:" << at_angle1 << " at_angle2:" << at_angle2<<endl;
	const int max_iter_num = 100;
	vector<OpenMesh::Vec3f> points;
	contour contour1 = *contours->operator[](contour1_index);
	if ((contour1).nextContours.size() <= 0)
	{
		//cout << "contour" << contour1.id << "沒有下一個contour.";
		fail = false;
		return points;
	}

	pointFromTwoSource* pair1 = contour1.searchAngle(at_angle1);
	glm::vec3 startPoint = toGlmVec3(pointFromTwoSource::Interpolation(pair1[0], pair1[1], at_angle1));
	//cout << "startPoint:(" << startPoint[0] << "," << startPoint[1] << "," << startPoint[2] << ")";
	MyMesh::FaceHandle face_start = contour1.getFaceAtAngle(mesh, at_angle1);
	glm::vec3 norm_start = glm::normalize(toGlmVec3(mesh.normal(face_start)));
	//cout << "起始點為:(" << startPoint[0] << "," << startPoint[1] << "," << startPoint[2] << ") 起始面:"<<face_start.idx()<<endl;
	contour contour2 = *((contours->operator[](contour1_index))->nextContours[0]);
	glm::vec3 endPoint = toGlmVec3(contour2.getPosAtAngle(at_angle2));
	MyMesh::FaceHandle face_end = contour2.getFaceAtAngle(mesh, at_angle2);
	glm::vec3 norm_end = glm::normalize(toGlmVec3(mesh.normal(face_end)));

	int endFaceIdx = contour2.getFaceAtAngle(mesh, at_angle2).idx();
	//cout << "終點fidx:" << endFaceIdx;
	points.push_back(OpenMesh::Vec3f(startPoint[0], startPoint[1], startPoint[2]));
	//cout << "終點為:(" << endPoint[0] << "," << endPoint[1] << "," << endPoint[2] << ") 終點面:" << face_end.idx()<<endl;
	//cout << "getPointsBtw " << contour1.id << " angle:" << at_angle << ":-------------------------------------------" << endl;

	glm::vec3 axis_x = glm::normalize(endPoint - startPoint);

	glm::vec3 temp_z = glm::normalize(norm_start + norm_end);
	//cout << "temp_z:(" << temp_z[0] << "," << temp_z[1] << "," << temp_z[2] << ")" << endl;
	glm::vec3 axis_y = glm::normalize(glm::cross(axis_x, temp_z));
	//cout << "axis_y:(" << axis_y[0] << "," << axis_y[1] << "," << axis_y[2] << ")" << endl;
	glm::vec3 axis_z = glm::normalize(glm::cross(axis_y, axis_x));
	//cout << "axis_z:(" << axis_z[0] << "," << axis_z[1] << "," << axis_z[2] << ")" << endl;

	coverAxis trans = coverAxis(startPoint, axis_x, axis_y, axis_z);
	MyMesh::FaceHandle nowFace = contour1.getFaceAtAngle(mesh, at_angle1);
	//cout << " 起點fidx:" << nowFace.idx();
	if (nowFace.idx() == endFaceIdx){//如果起點和終點在同一個面上,就不需要找了直接回傳起終點就好了
		points.push_back(OpenMesh::Vec3f(endPoint[0], endPoint[1], endPoint[2]));
		return points;
	}
	//cout << "起始面為:" << nowFace.idx() << endl;
	glm::vec3 p1 = trans.cover(toGlmVec3(pair1[0].pointPos));
	glm::vec3 p2 = trans.cover(toGlmVec3(pair1[1].pointPos));
	//cout << "covered pair1[0].pointpos:(" << p1[0] << "," << p1[1] << "," << p1[2] << ")"<<endl;
	//cout << "covered pair1[1].pointpos:(" << p2[0] << "," << p2[1] << "," << p2[2] << ")" << endl;
	glm::vec3 pe = trans.cover(endPoint);
	lastCover = trans;
	//cout << "covered endPoint:(" << pe[0] << "," << pe[1] << "," << pe[2] << ")";

	MyMesh::HalfedgeHandle mostPossible_he;
	OpenMesh::Vec3f most_interaction;
	float  most_avgX = -9999;
	//float coverDomain_y[2];
	for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(nowFace); fh_it != mesh.fh_end(nowFace); fh_it++){
		MyMesh::VertexHandle from = mesh.from_vertex_handle(*fh_it);
		MyMesh::VertexHandle to = mesh.to_vertex_handle(*fh_it);
		glm::vec3 point1 = trans.cover(toGlmVec3(mesh.point(from)));
		if (abs(point1[1])<getPathYAsZeroThreshold)//應該是因為精度問題當點太靠近轉換拿空間的y軸時,數值有可能會變成極小的複數導致本來找的到的交點變成找不到
		{
			point1[1] = 0;//為了讓演算法能運作,強制設置為0
		}
		glm::vec3 point2 = trans.cover(toGlmVec3(mesh.point(to)));
		if (abs(point2[1])<getPathYAsZeroThreshold)
		{
			point2[1] = 0;//為了讓演算法能運作,強制設置為0
		}
		//if (contour1_index==7)
			//cout << "起始point1:(" << point1[0] << "," << point1[1] << "," << point1[2] << ") point2:(" << point2[0] << "," << point2[1] << "," << point2[2] << ")";
		//if (contour1_index == 7)
			//cout << "起始邊:"<<from.idx()<<" (" << point1[0] << "," << point1[1] << "," << point1[2] << "); "<<to.idx()<<"(" << point2[0] << "," << point2[1] << "," << point2[2] << ");" << endl;
		bool isInteraction = (point1[1] >= 0 && point2[1] <= 0) || (point1[1] <= 0 && point2[1] >= 0);//如果有相交的話兩個點會有一個端點y值為正另一個端點y值為負
		if (isInteraction){
			glm::vec3 offset = point2 - point1;
			float time = (0 - point1[1]) / offset[1];
			glm::vec3 coverDomain_intersection = point1 + time*offset;//轉換空間的交點
			float avgX = coverDomain_intersection[0];//轉換空間交點的x值
			if (avgX>most_avgX){//我們取x最遠的轉換空間交點做為第一個交點
				mostPossible_he = *fh_it;
				most_avgX = avgX;
				//glm::vec3 P1toP2 = point2 - point1;
				//float time = (0 - point1[1]) / P1toP2[1];
				OpenMesh::Vec3f interaction = mesh.point(from) + (mesh.point(to) - mesh.point(from))*time;
				//if (contour1_index == 7)
				//	cout << "most_interaction v1:" << from.idx() << " v2:" << to.idx();
				most_interaction = interaction;
				//coverDomain_y[0] = point1.y;
				//coverDomain_y[1] = point2.y;
			}
		}
		//points.push_back(OpenMesh::Vec3f(interaction[0], interaction[1], interaction[2]));
	}
	//cout << "搜尋點開始:" << endl;
	points.push_back(most_interaction);
	MyMesh::HalfedgeHandle now_he = mesh.opposite_halfedge_handle(mostPossible_he);

	nowFace = mesh.opposite_face_handle(mostPossible_he);
	float totalPathLength = 0;
	while (nowFace.idx() != endFaceIdx){
		//cout << "半邊:" << (*fh_it).idx() << "反向半邊為:" << mesh.opposite_halfedge_handle(*fh_it).idx() << "反向面為:" << mesh.opposite_face_handle(*fh_it);

		//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
		//glm::vec3 point2 = trans.cover(toGlmVec3(to));
		MyMesh::VertexHandle now_he_p1 = mesh.from_vertex_handle(now_he);
		//if (nowFace.idx() == 7297)
		//cout << "now_he_p1:" << now_he_p1.idx();
		glm::vec3 now_he_p1_cover = trans.cover(toGlmVec3(mesh.point(now_he_p1)));
		if (abs(now_he_p1_cover[1]) < getPathYAsZeroThreshold)
		{
			now_he_p1_cover[1] = 0;//為了讓演算法能運作,強制設置為0
		}
		MyMesh::VertexHandle now_he_p2 = mesh.to_vertex_handle(now_he);
		//if ( nowFace.idx()==7297)
		//cout << "nowFace:" << nowFace.idx();
		//cout << " now_he_p2:" << now_he_p2.idx();
		//cout << " vec:" << mesh.point(now_he_p2);
		glm::vec3 now_he_p2_cover = trans.cover(toGlmVec3(mesh.point(now_he_p2)));
		//cout << "now_he_p1:" << now_he_p1 << "now_he_p2:" << now_he_p2<<endl;
		if (abs(now_he_p2_cover[1]) < getPathYAsZeroThreshold)
		{
			now_he_p2_cover[1] = 0;//為了讓演算法能運作,強制設置為0
		}
		OpenMesh::Vec3f lastPoint = points.back();
		glm::vec3 lastPointCovered = trans.cover(toGlmVec3(lastPoint));
		if (now_he_p1_cover[1] == 0)//如果上一個點是某個halfEdge的端點的話,就不能走普通的路線,因為當交點位於端點時,以這個點為中心向外的任何一條半邊都會滿足條件(因為其中一個點y==0,另一個點不論大於還是小於都會成立),最終導致交點一直卡在那個端點上
		{
			MyMesh::VertexHandle now_vh = now_he_p1;
			vector<MyMesh::VertexHandle> vNeighbors;
			for (MyMesh::VertexVertexCWIter vvcw_it = mesh.vv_cwbegin(now_vh); vvcw_it!=mesh.vv_cwend(now_vh); vvcw_it++)
			{
				vNeighbors.push_back(*vvcw_it);
			}
			float max_x=-9999;
			MyMesh::HalfedgeHandle mostPossible_he;
			OpenMesh::Vec3f most_intersection;
			//glm::vec3 most_coverDomain_intersection;
			for (size_t i = 0; i < vNeighbors.size(); i++)
			{
				OpenMesh::Vec3f from = mesh.point(vNeighbors[i]);
				OpenMesh::Vec3f to = mesh.point(vNeighbors[(i + 1) % vNeighbors.size()]);
				
				glm::vec3 point1 = trans.cover(toGlmVec3(from));
				if (abs(point1[1]) < getPathYAsZeroThreshold)//應該是因為精度問題當點太靠近轉換拿空間的y軸時,數值有可能會變成極小的複數導致本來找的到的交點變成找不到
				{
					point1[1] = 0;//為了讓演算法能運作,強制設置為0
				}
				//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
				glm::vec3 point2 = trans.cover(toGlmVec3(to));
				if (abs(point2[1]) < getPathYAsZeroThreshold)
				{
					point2[1] = 0;//為了讓演算法能運作,強制設置為0
				}
				//if (nowFace.idx() == 7297)
				//cout << "vNeighbors" << vNeighbors[i] << ":(" << point1[0] << "," << point1[1] << ") vNeighbors" << vNeighbors[(i + 1) % vNeighbors.size()] << ":(" << point2[0] << "," << point2[1] << ")";
				if ((point1[1]<=0&&point2[1]>=0)||(point1[1]>=0&&point2[1]<=0))
				{
					glm::vec3 P1toP2 = point2 - point1;
					float time = (0 - point1[1]) / P1toP2[1];
					glm::vec3 coverDomain_intersection = point1 + time*P1toP2;
					//if (nowFace.idx() == 7297)
					//cout << "coverd:(" << coverDomain_intersection[0] << "," << coverDomain_intersection[1] << ")";
					//glm::vec3 interaction = point1 + P1toP2*time;
					if (coverDomain_intersection[0]>max_x)
					{
						//if (nowFace.idx() == 7297)
						//cout << " set!";
						most_intersection = from + (to - from)*time;
						max_x = coverDomain_intersection[0];
						mostPossible_he = mesh.find_halfedge(vNeighbors[i], vNeighbors[(i + 1) % vNeighbors.size()]);
					}
				}
				//cout << endl;
			}

			if (max_x>pe[0])//如果延伸距離超過終點所在的距離
			{
				float time = (pe[0] - lastPointCovered[0]) / (max_x - lastPointCovered[0]);
				if (time < 0){//早已超過endPoint
					fail = true;
					return points;
				}
				OpenMesh::Vec3f possibleEndPoint = time*(most_intersection - lastPoint) + lastPoint;
				if (glm::length(toGlmVec3(possibleEndPoint) - endPoint)<matchEndPointThreshold){//如果可能的交點非常靠近終點了,就認為它是終點
					points.push_back(possibleEndPoint);
					break;
				}
			}
			now_he = mostPossible_he;//這裡不用反向
			nowFace = mesh.face_handle(now_he);
			points.push_back(most_intersection);
			//cout << "over" << endl;
		}
		else if (now_he_p2_cover[1] == 0){
			MyMesh::VertexHandle now_vh = now_he_p2;
			vector<MyMesh::VertexHandle> vNeighbors;
			for (MyMesh::VertexVertexCWIter vvcw_it = mesh.vv_cwbegin(now_vh); vvcw_it != mesh.vv_cwend(now_vh); vvcw_it++)
			{
				vNeighbors.push_back(*vvcw_it);
			}
			float max_x = -9999;
			MyMesh::HalfedgeHandle mostPossible_he;
			OpenMesh::Vec3f most_intersection;
			//glm::vec3 most_coverDomain_intersection;
			for (size_t i = 0; i < vNeighbors.size(); i++)
			{
				OpenMesh::Vec3f from = mesh.point(vNeighbors[i]);
				OpenMesh::Vec3f to = mesh.point(vNeighbors[(i + 1) % vNeighbors.size()]);

				glm::vec3 point1 = trans.cover(toGlmVec3(from));
				if (abs(point1[1]) < getPathYAsZeroThreshold)//應該是因為精度問題當點太靠近轉換拿空間的y軸時,數值有可能會變成極小的負數導致本來找的到的交點變成找不到
				{
					point1[1] = 0;//為了讓演算法能運作,強制設置為0
				}
				//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
				glm::vec3 point2 = trans.cover(toGlmVec3(to));
				if (abs(point2[1]) < getPathYAsZeroThreshold)
				{
					point2[1] = 0;//為了讓演算法能運作,強制設置為0
				}
				//if (nowFace.idx() == 7297)
				//cout << "vNeighbors" << vNeighbors[i] << ":(" << point1[0] << "," << point1[1] << ") vNeighbors" << vNeighbors[(i + 1) % vNeighbors.size()] << ":(" << point2[0] << "," << point2[1] << ")";
				if ((point1[1] <= 0 && point2[1] >= 0) || (point1[1] >= 0 && point2[1] <= 0))
				{
					glm::vec3 P1toP2 = point2 - point1;
					float time = (0 - point1[1]) / P1toP2[1];
					glm::vec3 coverDomain_intersection = point1 + time*P1toP2;
					//if (nowFace.idx() == 7297)
					//cout << "coverd:(" << coverDomain_intersection[0] << "," << coverDomain_intersection[1] << ")";
					//glm::vec3 interaction = point1 + P1toP2*time;
					if (coverDomain_intersection[0]>max_x)
					{
						//if (nowFace.idx() == 7297)
						//cout << " set!";
						most_intersection = from + (to - from)*time;
						max_x = coverDomain_intersection[0];
						mostPossible_he = mesh.find_halfedge(vNeighbors[i], vNeighbors[(i + 1) % vNeighbors.size()]);
					}
				}
				//cout << endl;
			}

			if (max_x>pe[0])//如果延伸距離超過終點所在的距離
			{
				float time = (pe[0] - lastPointCovered[0]) / (max_x - lastPointCovered[0]);
				if (time < 0){//早已超過endPoint
					fail = true;
					return points;
				}
				OpenMesh::Vec3f possibleEndPoint = time*(most_intersection - lastPoint) + lastPoint;
				if (glm::length(toGlmVec3(possibleEndPoint) - endPoint)<matchEndPointThreshold){//如果可能的交點非常靠近終點了,就認為它是終點
					points.push_back(possibleEndPoint);
					break;
				}
			}
			now_he = mostPossible_he;//這裡不用反向
			nowFace = mesh.face_handle(now_he);
			points.push_back(most_intersection);
			//cout << "over" << endl;
		}
		else
		{
			for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(nowFace); fh_it != mesh.fh_end(nowFace); fh_it++){

				if ((*fh_it) != now_he){//忽略更新的來源邊
					//cout << "半邊:" << (*fh_it).idx() << "反向半邊為:" << mesh.opposite_halfedge_handle(*fh_it).idx() << "反向面為:" << mesh.opposite_face_handle(*fh_it);
					OpenMesh::Vec3f from = mesh.point(mesh.from_vertex_handle(*fh_it));
					OpenMesh::Vec3f to = mesh.point(mesh.to_vertex_handle(*fh_it));
					glm::vec3 point1 = trans.cover(toGlmVec3(from));
					if (abs(point1[1]) < getPathYAsZeroThreshold)//應該是因為精度問題當點太靠近轉換拿空間的y軸時,數值有可能會變成極小的複數導致本來找的到的交點變成找不到
					{
						point1[1] = 0;//為了讓演算法能運作,強制設置為0
					}
					//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
					glm::vec3 point2 = trans.cover(toGlmVec3(to));
					if (abs(point2[1]) < getPathYAsZeroThreshold)
					{
						point2[1] = 0;//為了讓演算法能運作,強制設置為0
					}
					//if (contour1_index == 7)
					//	cout << "interaction v1:" << mesh.from_vertex_handle(*fh_it) << ">(" << point1[0] << "," << point1[1] << ") v2:" << mesh.to_vertex_handle(*fh_it) << ">(" << point2[0] << "," << point2[1] << ")";

					//cout << "point2 idx"<<mesh.to_vertex_handle(*fh_it).idx()<<":(" << point2[0] << "," << point2[1] << "," << point2[2] << ") ";

					if ((point1[1] >= 0 && point2[1] <= 0) || (point1[1] <= 0 && point2[1]>= 0)){//找到一個符合兩端的各在x軸一邊的邊
						//cout << "從半邊:" << (*fh_it).idx();
						glm::vec3 P1toP2 = point2 - point1;
						float time = (0 - point1[1]) / P1toP2[1];
						//glm::vec3 interaction = point1 + P1toP2*time;
						OpenMesh::Vec3f intersection = from + (to - from)*time;
						//glm::vec3 intersection_trans = point1 + (point2 - point1)*time;
						//totalPathLength += (intersection - points.back()).length();

						points.push_back(intersection);
						/*if (totalPathLength > pathLengthThreshold){//說明超過了終點所在的範圍,在雙層結構下很可能會發生,終點在上層,路徑卻走下層
							fail = true;
							return points;
						}*/
						glm::vec3 coverDomain_intersection = point1 + time*P1toP2;
						if (coverDomain_intersection[0] - pe[0]){
							float time = (pe[0] - lastPointCovered[0]) / (coverDomain_intersection[0] - lastPointCovered[0]);
							if (time < 0){//早已超過endPoint
								fail = true;
								return points;
							}
						}
						//cout << "找到點(" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")";

						now_he = mesh.opposite_halfedge_handle(*fh_it);//反向當前半邊
						nowFace = mesh.face_handle(now_he);//推進面
						break;
					}
				}
				//cout << ";     ";
				//cout << endl;
			}
		}
	}
		//cout << endl;
	points.push_back(OpenMesh::Vec3f(endPoint[0], endPoint[1], endPoint[2]));
	//cout << "getPointsBtw end----------------------------------------------"<<endl;
	fail = false;
	return points;
}
vector<OpenMesh::Vec3f> polarMap::getPointsBtw(int contour1_index,int next_contour_index, float at_angle1, float at_angle2, bool &fail){//next_contour_index 並不是下一條contour的index,而是下一條contour位於nextContours的索引值
	//cout << "測試-9.99 loopClampTo01:" << loopClampTo01(-9.99) << endl;
	//cout << "測試 100.0 loopClampTo01:" << loopClampTo01(1000.0) << endl;
	//cout << "測試- 1.478 loopClampTo01:" << loopClampTo01(-1.478) << endl;
	at_angle1 = loopClampTo01(at_angle1);
	at_angle2 = loopClampTo01(at_angle2);
	//cout << "clamp to 01后 at_angle1:" << at_angle1 << " at_angle2:" << at_angle2<<endl;
	const int max_iter_num = 100;
	vector<OpenMesh::Vec3f> points;
	contour contour1 = *contours->operator[](contour1_index);
	if ((contour1).nextContours.size() <= 0)
	{
		cout << "contour" << contour1.id << "沒有下一個contour.";
		fail = false;
		return points;
	}
	pointFromTwoSource* pair1 = contour1.searchAngle(at_angle1);
	glm::vec3 startPoint = toGlmVec3(pointFromTwoSource::Interpolation(pair1[0], pair1[1], at_angle1));
	MyMesh::FaceHandle face_start = contour1.getFaceAtAngle(mesh, at_angle1);
	glm::vec3 norm_start = glm::normalize(toGlmVec3(mesh.normal(face_start)));
	//cout << "起始點為:(" << startPoint[0] << "," << startPoint[1] << "," << startPoint[2] << ") 起始面:"<<face_start.idx()<<endl;
	//cout << "getPointsBtw: contour1_index:" << contour1_index << "next_contour_index:" << next_contour_index << endl;
	contour contour2 = *((contours->operator[](contour1_index))->nextContours[next_contour_index]);
	glm::vec3 endPoint = toGlmVec3(contour2.getPosAtAngle(at_angle2));
	MyMesh::FaceHandle face_end = contour2.getFaceAtAngle(mesh, at_angle2);
	glm::vec3 norm_end = glm::normalize(toGlmVec3(mesh.normal(face_end)));

	int endFaceIdx = contour2.getFaceAtAngle(mesh, at_angle2).idx();
	points.push_back(OpenMesh::Vec3f(startPoint[0], startPoint[1], startPoint[2]));
	//cout << "終點為:(" << endPoint[0] << "," << endPoint[1] << "," << endPoint[2] << ") 終點面:" << face_end.idx()<<endl;
	//cout << "getPointsBtw " << contour1.id << " angle:" << at_angle << ":-------------------------------------------" << endl;

	glm::vec3 axis_x = glm::normalize(endPoint - startPoint);

	glm::vec3 temp_z = glm::normalize(norm_start + norm_end);
	//cout << "temp_z:(" << temp_z[0] << "," << temp_z[1] << "," << temp_z[2] << ")" << endl;
	glm::vec3 axis_y = glm::normalize(glm::cross(axis_x, temp_z));
	//cout << "axis_y:(" << axis_y[0] << "," << axis_y[1] << "," << axis_y[2] << ")" << endl;
	glm::vec3 axis_z = glm::normalize(glm::cross(axis_y, axis_x));
	//cout << "axis_z:(" << axis_z[0] << "," << axis_z[1] << "," << axis_z[2] << ")" << endl;

	coverAxis trans = coverAxis(startPoint, axis_x, axis_y, axis_z);
	MyMesh::FaceHandle nowFace = contour1.getFaceAtAngle(mesh, at_angle1);
	//cout << "起始面為:" << nowFace.idx() << endl;
	glm::vec3 p1 = trans.cover(toGlmVec3(pair1[0].pointPos));
	glm::vec3 p2 = trans.cover(toGlmVec3(pair1[1].pointPos));
	//cout << "covered pair1[0].pointpos:(" << p1[0] << "," << p1[1] << "," << p1[2] << ")"<<endl;
	//cout << "covered pair1[1].pointpos:(" << p2[0] << "," << p2[1] << "," << p2[2] << ")" << endl;
	glm::vec3 pe = trans.cover(endPoint);
	lastCover = trans;
	//cout << "covered endPoint:(" << pe[0] << "," << pe[1] << "," << pe[2] << ")";

	MyMesh::HalfedgeHandle mostPossible_he;
	OpenMesh::Vec3f most_interaction;
	float  most_avgX = -9999;
	for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(nowFace); fh_it != mesh.fh_end(nowFace); fh_it++){
		MyMesh::VertexHandle from = mesh.from_vertex_handle(*fh_it);
		MyMesh::VertexHandle to = mesh.to_vertex_handle(*fh_it);
		glm::vec3 point1 = trans.cover(toGlmVec3(mesh.point(from)));
		glm::vec3 point2 = trans.cover(toGlmVec3(mesh.point(to)));
		//cout << "起始邊:"<<from.idx()<<" (" << point1[0] << "," << point1[1] << "," << point1[2] << "); "<<to.idx()<<"(" << point2[0] << "," << point2[1] << "," << point2[2] << ");" << endl;
		bool isInteraction = (point1[1] >= 0 && point2[1] <= 0) || (point1[1] <= 0 && point2[1] >= 0);//如果有相交的話兩個點會有一個端點y值為正另一個端點y值為負
		if (isInteraction){
			glm::vec3 offset = point2 - point1;
			float time = (0 - point1[1]) / offset[1];
			glm::vec3 coverDomain_intersection = point1 + time*offset;
			float avgX = coverDomain_intersection[0];
			if (avgX>most_avgX){
				mostPossible_he = *fh_it;
				most_avgX = avgX;
				//glm::vec3 P1toP2 = point2 - point1;
				//float time = (0 - point1[1]) / P1toP2[1];
				OpenMesh::Vec3f interaction = mesh.point(from) + (mesh.point(to) - mesh.point(from))*time;
				most_interaction = interaction;
			}
		}
		//points.push_back(OpenMesh::Vec3f(interaction[0], interaction[1], interaction[2]));
	}
	//cout << "搜尋點開始:" << endl;
	points.push_back(most_interaction);
	MyMesh::HalfedgeHandle now_he = mesh.opposite_halfedge_handle(mostPossible_he);
	nowFace = mesh.opposite_face_handle(mostPossible_he);

	float totalPathLength = 0;
	/*while (nowFace.idx() != endFaceIdx){
		for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(nowFace); fh_it != mesh.fh_end(nowFace); fh_it++){

			if ((*fh_it) != now_he){//忽略更新的來源邊
				//cout << "半邊:" << (*fh_it).idx() << "反向半邊為:" << mesh.opposite_halfedge_handle(*fh_it).idx() << "反向面為:" << mesh.opposite_face_handle(*fh_it);
				OpenMesh::Vec3f from = mesh.point(mesh.from_vertex_handle(*fh_it));
				OpenMesh::Vec3f to = mesh.point(mesh.to_vertex_handle(*fh_it));
				glm::vec3 point1 = trans.cover(toGlmVec3(from));
				//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
				glm::vec3 point2 = trans.cover(toGlmVec3(to));
				//cout << "point2 idx"<<mesh.to_vertex_handle(*fh_it).idx()<<":(" << point2[0] << "," << point2[1] << "," << point2[2] << ") ";
				if ((point1[1] >= 0 && point2[1] <= 0) || (point1[1] <= 0 && point2[1] >= 0)){//找到一個符合兩端的各在x軸一邊的邊
					//cout << "從半邊:" << (*fh_it).idx();
					glm::vec3 P1toP2 = point2 - point1;
					float time = (0 - point1[1]) / P1toP2[1];
					//glm::vec3 interaction = point1 + P1toP2*time;
					OpenMesh::Vec3f intersection = from + (to - from)*time;
					//glm::vec3 intersection_trans = point1 + (point2 - point1)*time;
					totalPathLength += (intersection - points.back()).length();

					points.push_back(intersection);
					//cout << "(" << intersection_trans[0] << "," << intersection_trans[1] << "," << intersection_trans[2] << ")";

					if (totalPathLength>pathLengthThreshold){//說明超過了終點所在的範圍,在雙層結構下很可能會發生,終點在上層,路徑卻走下層
						fail = true;
						return points;
					}
					//cout << "找到點(" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")";

					now_he = mesh.opposite_halfedge_handle(*fh_it);//反向當前半邊
					//cout << "更新半邊到:" << now_he.idx() << endl;
					nowFace = mesh.face_handle(now_he);//推進面
					break;
				}
			}
		}
	}*/
	while (nowFace.idx() != endFaceIdx){
		MyMesh::VertexHandle now_he_p1 = mesh.from_vertex_handle(now_he);
		//if (nowFace.idx() == 7297)
		//cout << "now_he_p1:" << now_he_p1.idx();
		glm::vec3 now_he_p1_cover = trans.cover(toGlmVec3(mesh.point(now_he_p1)));
		if (abs(now_he_p1_cover[1]) < getPathYAsZeroThreshold)
		{
			now_he_p1_cover[1] = 0;//為了讓演算法能運作,強制設置為0
		}
		MyMesh::VertexHandle now_he_p2 = mesh.to_vertex_handle(now_he);
		glm::vec3 now_he_p2_cover = trans.cover(toGlmVec3(mesh.point(now_he_p2)));
		//cout << "now_he_p1:" << now_he_p1 << "now_he_p2:" << now_he_p2<<endl;
		if (abs(now_he_p2_cover[1]) < getPathYAsZeroThreshold)
		{
			now_he_p2_cover[1] = 0;//為了讓演算法能運作,強制設置為0
		}
		OpenMesh::Vec3f lastPoint = points.back();
		glm::vec3 lastPointCovered = trans.cover(toGlmVec3(lastPoint));
		if (now_he_p1_cover[1] == 0)//如果上一個點是某個halfEdge的端點的話,就不能走普通的路線,因為當交點位於端點時,以這個點為中心向外的任何一條半邊都會滿足條件(因為其中一個點y==0,另一個點不論大於還是小於都會成立),最終導致交點一直卡在那個端點上
		{
			MyMesh::VertexHandle now_vh = now_he_p1;
			vector<MyMesh::VertexHandle> vNeighbors;
			for (MyMesh::VertexVertexCWIter vvcw_it = mesh.vv_cwbegin(now_vh); vvcw_it != mesh.vv_cwend(now_vh); vvcw_it++)
			{
				vNeighbors.push_back(*vvcw_it);
			}
			float max_x = -9999;
			MyMesh::HalfedgeHandle mostPossible_he;
			OpenMesh::Vec3f most_intersection;
			//glm::vec3 most_coverDomain_intersection;
			for (size_t i = 0; i < vNeighbors.size(); i++)
			{
				OpenMesh::Vec3f from = mesh.point(vNeighbors[i]);
				OpenMesh::Vec3f to = mesh.point(vNeighbors[(i + 1) % vNeighbors.size()]);

				glm::vec3 point1 = trans.cover(toGlmVec3(from));
				if (abs(point1[1]) < getPathYAsZeroThreshold)//應該是因為精度問題當點太靠近轉換拿空間的y軸時,數值有可能會變成極小的複數導致本來找的到的交點變成找不到
				{
					point1[1] = 0;//為了讓演算法能運作,強制設置為0
				}
				//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
				glm::vec3 point2 = trans.cover(toGlmVec3(to));
				if (abs(point2[1]) < getPathYAsZeroThreshold)
				{
					point2[1] = 0;//為了讓演算法能運作,強制設置為0
				}
				//if (nowFace.idx() == 7297)
				//cout << "vNeighbors" << vNeighbors[i] << ":(" << point1[0] << "," << point1[1] << ") vNeighbors" << vNeighbors[(i + 1) % vNeighbors.size()] << ":(" << point2[0] << "," << point2[1] << ")";
				if ((point1[1] <= 0 && point2[1] >= 0) || (point1[1] >= 0 && point2[1] <= 0))
				{
					glm::vec3 P1toP2 = point2 - point1;
					float time = (0 - point1[1]) / P1toP2[1];
					glm::vec3 coverDomain_intersection = point1 + time*P1toP2;
					if (coverDomain_intersection[0]>max_x)
					{
						//if (nowFace.idx() == 7297)
						//cout << " set!";
						most_intersection = from + (to - from)*time;
						max_x = coverDomain_intersection[0];
						mostPossible_he = mesh.find_halfedge(vNeighbors[i], vNeighbors[(i + 1) % vNeighbors.size()]);
					}
				}
				//cout << endl;
			}

			if (max_x>pe[0])//如果延伸距離超過終點所在的距離
			{
				float time = (pe[0] - lastPointCovered[0]) / (max_x - lastPointCovered[0]);
				if (time < 0){//早已超過endPoint
					fail = true;
					return points;
				}
				OpenMesh::Vec3f possibleEndPoint = time*(most_intersection - lastPoint) + lastPoint;
				if (glm::length(toGlmVec3(possibleEndPoint) - endPoint)<matchEndPointThreshold){//如果可能的交點非常靠近終點了,就認為它是終點
					points.push_back(possibleEndPoint);
					break;
				}
			}
			now_he = mostPossible_he;//這裡不用反向
			nowFace = mesh.face_handle(now_he);
			points.push_back(most_intersection);
			//cout << "over" << endl;
		}
		else if (now_he_p2_cover[1] == 0){
			MyMesh::VertexHandle now_vh = now_he_p2;
			vector<MyMesh::VertexHandle> vNeighbors;
			for (MyMesh::VertexVertexCWIter vvcw_it = mesh.vv_cwbegin(now_vh); vvcw_it != mesh.vv_cwend(now_vh); vvcw_it++)
			{
				vNeighbors.push_back(*vvcw_it);
			}
			float max_x = -9999;
			MyMesh::HalfedgeHandle mostPossible_he;
			OpenMesh::Vec3f most_intersection;
			//glm::vec3 most_coverDomain_intersection;
			for (size_t i = 0; i < vNeighbors.size(); i++)
			{
				OpenMesh::Vec3f from = mesh.point(vNeighbors[i]);
				OpenMesh::Vec3f to = mesh.point(vNeighbors[(i + 1) % vNeighbors.size()]);

				glm::vec3 point1 = trans.cover(toGlmVec3(from));
				if (abs(point1[1]) < getPathYAsZeroThreshold)//應該是因為精度問題當點太靠近轉換拿空間的y軸時,數值有可能會變成極小的負數導致本來找的到的交點變成找不到
				{
					point1[1] = 0;//為了讓演算法能運作,強制設置為0
				}
				//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
				glm::vec3 point2 = trans.cover(toGlmVec3(to));
				if (abs(point2[1]) < getPathYAsZeroThreshold)
				{
					point2[1] = 0;//為了讓演算法能運作,強制設置為0
				}
				//if (nowFace.idx() == 7297)
				//cout << "vNeighbors" << vNeighbors[i] << ":(" << point1[0] << "," << point1[1] << ") vNeighbors" << vNeighbors[(i + 1) % vNeighbors.size()] << ":(" << point2[0] << "," << point2[1] << ")";
				if ((point1[1] <= 0 && point2[1] >= 0) || (point1[1] >= 0 && point2[1] <= 0))
				{
					glm::vec3 P1toP2 = point2 - point1;
					float time = (0 - point1[1]) / P1toP2[1];
					glm::vec3 coverDomain_intersection = point1 + time*P1toP2;
					//if (nowFace.idx() == 7297)
					//cout << "coverd:(" << coverDomain_intersection[0] << "," << coverDomain_intersection[1] << ")";
					//glm::vec3 interaction = point1 + P1toP2*time;
					if (coverDomain_intersection[0]>max_x)
					{
						//if (nowFace.idx() == 7297)
						//cout << " set!";
						most_intersection = from + (to - from)*time;
						max_x = coverDomain_intersection[0];
						mostPossible_he = mesh.find_halfedge(vNeighbors[i], vNeighbors[(i + 1) % vNeighbors.size()]);
					}
				}
				//cout << endl;
			}

			if (max_x>pe[0])//如果延伸距離超過終點所在的距離
			{
				float time = (pe[0] - lastPointCovered[0]) / (max_x - lastPointCovered[0]);
				if (time < 0){//早已超過endPoint
					fail = true;
					return points;
				}
				OpenMesh::Vec3f possibleEndPoint = time*(most_intersection - lastPoint) + lastPoint;
				if (glm::length(toGlmVec3(possibleEndPoint) - endPoint)<matchEndPointThreshold){//如果可能的交點非常靠近終點了,就認為它是終點
					points.push_back(possibleEndPoint);
					break;
				}
			}
			now_he = mostPossible_he;//這裡不用反向
			nowFace = mesh.face_handle(now_he);
			points.push_back(most_intersection);
			//cout << "over" << endl;
		}
		else
		{
			for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(nowFace); fh_it != mesh.fh_end(nowFace); fh_it++){

				if ((*fh_it) != now_he){//忽略更新的來源邊
					//cout << "半邊:" << (*fh_it).idx() << "反向半邊為:" << mesh.opposite_halfedge_handle(*fh_it).idx() << "反向面為:" << mesh.opposite_face_handle(*fh_it);
					OpenMesh::Vec3f from = mesh.point(mesh.from_vertex_handle(*fh_it));
					OpenMesh::Vec3f to = mesh.point(mesh.to_vertex_handle(*fh_it));
					glm::vec3 point1 = trans.cover(toGlmVec3(from));
					if (abs(point1[1]) < getPathYAsZeroThreshold)//應該是因為精度問題當點太靠近轉換拿空間的y軸時,數值有可能會變成極小的複數導致本來找的到的交點變成找不到
					{
						point1[1] = 0;//為了讓演算法能運作,強制設置為0
					}
					//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
					glm::vec3 point2 = trans.cover(toGlmVec3(to));
					if (abs(point2[1]) < getPathYAsZeroThreshold)
					{
						point2[1] = 0;//為了讓演算法能運作,強制設置為0
					}
					//if (contour1_index == 7)
					//	cout << "interaction v1:" << mesh.from_vertex_handle(*fh_it) << ">(" << point1[0] << "," << point1[1] << ") v2:" << mesh.to_vertex_handle(*fh_it) << ">(" << point2[0] << "," << point2[1] << ")";

					//cout << "point2 idx"<<mesh.to_vertex_handle(*fh_it).idx()<<":(" << point2[0] << "," << point2[1] << "," << point2[2] << ") ";

					if ((point1[1] >= 0 && point2[1] <= 0) || (point1[1] <= 0 && point2[1] >= 0)){//找到一個符合兩端的各在x軸一邊的邊
						//cout << "從半邊:" << (*fh_it).idx();
						glm::vec3 P1toP2 = point2 - point1;
						float time = (0 - point1[1]) / P1toP2[1];
						//glm::vec3 interaction = point1 + P1toP2*time;
						OpenMesh::Vec3f intersection = from + (to - from)*time;
						//glm::vec3 intersection_trans = point1 + (point2 - point1)*time;
						//totalPathLength += (intersection - points.back()).length();

						points.push_back(intersection);
						/*if (totalPathLength > pathLengthThreshold){//說明超過了終點所在的範圍,在雙層結構下很可能會發生,終點在上層,路徑卻走下層
						fail = true;
						return points;
						}*/
						glm::vec3 coverDomain_intersection = point1 + time*P1toP2;
						if (coverDomain_intersection[0] - pe[0]){
							float time = (pe[0] - lastPointCovered[0]) / (coverDomain_intersection[0] - lastPointCovered[0]);
							if (time < 0){//早已超過endPoint
								fail = true;
								return points;
							}
						}
						//cout << "找到點(" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")";

						now_he = mesh.opposite_halfedge_handle(*fh_it);//反向當前半邊
						nowFace = mesh.face_handle(now_he);//推進面
						break;
					}
				}
				//cout << ";     ";
				//cout << endl;
			}
		}
	}
	points.push_back(OpenMesh::Vec3f(endPoint[0], endPoint[1], endPoint[2]));
	//cout << "getPointsBtw end----------------------------------------------"<<endl;
	fail = false;
	return points;
}
OpenMesh::Vec3f polarMap::getPointAtBtw(int contour1_index, float at_angle, float percent){
	vector<OpenMesh::Vec3f> points = getPointsBtw(contour1_index, at_angle);
	float totalLength = 0;
	for (int i = 0; i < points.size() - 1; i++){
		totalLength += (points[i + 1] - points[i]).length();
	}
	float nowLength = 0;
	for (int i = 0; i < points.size() - 1; i++){
		float len = (points[i + 1] - points[i]).length();
		if ((nowLength + len) / totalLength >= percent){//點就在這個範圍內
			float sub = ((percent*totalLength) - nowLength) / len;
			return points[i]+ (points[i + 1] - points[i])*sub;
		}
		else
		{
			nowLength += len;
		}
	}
	cout << "percent:" << percent << "搜尋不到點!"<<endl;
	return points.back();
}
OpenMesh::Vec3f polarMap::getPointAtBtw(int contour1_index, float at_angle, float percent, bool &fail){
	vector<OpenMesh::Vec3f> points = getPointsBtw(contour1_index, at_angle,fail);
	float totalLength = 0;
	for (int i = 0; i < points.size() - 1; i++){
		totalLength += (points[i + 1] - points[i]).length();
	}
	float nowLength = 0;
	for (int i = 0; i < points.size() - 1; i++){
		float len = (points[i + 1] - points[i]).length();
		if ((nowLength + len) / totalLength >= percent){//點就在這個範圍內
			float sub = ((percent*totalLength) - nowLength) / len;
			return points[i] + (points[i + 1] - points[i])*sub;
		}
		else
		{
			nowLength += len;
		}
	}
	cout << "percent:" << percent << "搜尋不到點!" << endl;
	return points.back();
}
OpenMesh::Vec3f polarMap::getPointAtBtw(int contour1_index, float at_angle1, float at_angle2,float percent, bool &fail){
	vector<OpenMesh::Vec3f> points = getPointsBtw(contour1_index, at_angle1, at_angle2, fail);
	if (fail){//如果失敗了就直接結束了
		//points = getgdcPointsBtw(contour1_index, at_angle1, at_angle2);
		if (contour1_index ==3)
		{
			pathRecord = points;
		}
		return OpenMesh::Vec3f(0, 0, 0);
	}
	float totalLength = 0;
	for (int i = 0; i < points.size() - 1; i++){
		totalLength += (points[i + 1] - points[i]).length();
	}
	float nowLength = 0;
	for (int i = 0; i < points.size() - 1; i++){
		float len = (points[i + 1] - points[i]).length();
		if ((nowLength + len) / totalLength >= percent){//點就在這個範圍內
			float sub = ((percent*totalLength) - nowLength) / len;
			return points[i] + (points[i + 1] - points[i])*sub;
		}
		else
		{
			nowLength += len;
		}
	}
	cout << "percent:" << percent << "搜尋不到點!" << endl;
	return points.back();
}
OpenMesh::Vec3f polarMap::getPointAtBtw(int contour1_index, int next_contour_index, float at_angle1, float at_angle2, float percent, bool &fail){
	vector<OpenMesh::Vec3f> points = getPointsBtw(contour1_index,next_contour_index,at_angle1, at_angle2, fail);
	if (fail){//如果失敗了就直接結束了
		//points = getgdcPointsBtw(contour1_index, at_angle1, at_angle2);
		return OpenMesh::Vec3f(0, 0, 0);
	}
	float totalLength = 0;
	for (int i = 0; i < points.size() - 1; i++){
		totalLength += (points[i + 1] - points[i]).length();
	}
	float nowLength = 0;
	for (int i = 0; i < points.size() - 1; i++){
		float len = (points[i + 1] - points[i]).length();
		if ((nowLength + len) / totalLength >= percent){//點就在這個範圍內
			float sub = ((percent*totalLength) - nowLength) / len;
			return points[i] + (points[i + 1] - points[i])*sub;
		}
		else
		{
			nowLength += len;
		}
	}
	cout << "percent:" << percent << "搜尋不到點!" << endl;
	return points.back();
}
void calVertexDist(int now_vidx,MyMesh mesh, vector<int> domain, vector<int> &v_dist, vector<int> &v_last_vidx){

}
vector<OpenMesh::Vec3f> polarMap::getPointsBtw_gdcPath(int contour1_index, float at_angle1, float at_angle2, vector<int> &shortestPath){
	explorePossibleBoundary.clear();
	contour *c = contours->operator[](contour1_index);
	explorePossibleBoundary.push_back(c);
	for each (contour* next in c->nextContours)
	{
		explorePossibleBoundary.push_back(next);
	}
	contour *c2 = c->nextContours[0];
	//找到離起始點最近範圍內的頂點
	pointFromTwoSource* startPair = c->searchAngle(at_angle1);
	OpenMesh::Vec3f startPos = c->getPosAtAngle(at_angle1);
	float oriDist;
	if (startPair[0].sourceIdx_Large == startPair[1].sourceIdx_Large){
		exploreOrigenVIdx = startPair[0].sourceIdx_Large;
		oriDist = (mesh.point(OpenMesh::VertexHandle(exploreOrigenVIdx))-startPos).length();
	}
	else{
		OpenMesh::Vec3f v1 = mesh.point(OpenMesh::VertexHandle(startPair[0].sourceIdx_Large))-startPos;
		OpenMesh::Vec3f v2 = mesh.point(OpenMesh::VertexHandle(startPair[1].sourceIdx_Large))-startPos;
		if (v1.length()<v2.length())
		{
			exploreOrigenVIdx = startPair[0].sourceIdx_Large;
			oriDist = v1.length();
		}
		else{
			exploreOrigenVIdx = startPair[1].sourceIdx_Large;
			oriDist = v2.length();
		}
	}
	//找到離終點最近範圍內的頂點
	pointFromTwoSource *endPair = c2->searchAngle(at_angle2);
	OpenMesh::Vec3f endPos = c2->getPosAtAngle(at_angle2);
	int endVertexIdx;
	float endDist;
	if (endPair[0].sourceIdx_Large == endPair[1].sourceIdx_Large){
		endVertexIdx = endPair[0].sourceIdx_Small;
		endDist = (mesh.point(OpenMesh::VertexHandle(endVertexIdx)) - endPos).length();
	}
	else
	{
		OpenMesh::Vec3f v1 = mesh.point(OpenMesh::VertexHandle(endPair[0].sourceIdx_Small));
		OpenMesh::Vec3f v2 = mesh.point(OpenMesh::VertexHandle(endPair[1].sourceIdx_Small));
		if (v1.length() < v2.length()){
			endVertexIdx = endPair[0].sourceIdx_Small;
			endDist = v1.length();
		}
		else{
			endVertexIdx = endPair[0].sourceIdx_Small;
			endDist = v2.length();
		}
	}
	vector<int> vIdx;
	float* vidx_dist;
	int* last_vidx;
	set<int> cIdx;
	exploreVertex(exploreOrigenVIdx, vIdx, cIdx,false);
	vidx_dist = new float[mesh.n_vertices()];
	for (int i = 0; i <mesh.n_vertices(); i++){//初始化距離矩陣
		vidx_dist[i] = BIG_NUMBER;
	}
	vidx_dist[exploreOrigenVIdx] = 0;
	last_vidx = new int[mesh.n_vertices()];
	vector<OpenMesh::VertexHandle> candiate;
	candiate.push_back(OpenMesh::VertexHandle(exploreOrigenVIdx));
	/*cout << "exploreOrigenVIdx:" << exploreOrigenVIdx << endl;
	cout << "endVertexIdx:" << endVertexIdx << endl;
	for each (int vidx  in vIdx)
	{
		cout << vidx << ",";
	}
	if (contain(vIdx, endVertexIdx)){
		cout << "終止點位於vIdx";
	}
	else{
		cout << "終止點不在vIdx";
	}*/
	while (candiate.size() > 0){
		calRingVertex(candiate, vidx_dist, last_vidx, vIdx);
	}
	shortestPath.clear();
	bool meetStart = false;
	int now_vidx = endVertexIdx;
	shortestPath.push_back(now_vidx);
	while (true){
		//cout << "path insert:" << last_vidx[now_vidx]<<",";
		shortestPath.insert(shortestPath.begin(), last_vidx[now_vidx]);
		now_vidx = last_vidx[now_vidx];
		if (now_vidx == exploreOrigenVIdx){
			shortestPath.insert(shortestPath.begin(), exploreOrigenVIdx);
			break;
		}
	}
	OpenMesh::Vec3f totalVec(0,0,0);
	for each(int vidx in shortestPath){
		OpenMesh::Vec3f vpos = mesh.point(OpenMesh::VertexHandle(vidx));
		totalVec += (vpos-startPos);
	}
	vector<OpenMesh::Vec3f> points;
	glm::vec3 temp_axis_z = toGlmVec3(totalVec.normalize());
	glm::vec3 axis_x = toGlmVec3((endPos- startPos).normalize());
	glm::vec3 axis_y = glm::normalize(glm::cross(temp_axis_z, axis_x));
	glm::vec3 axis_z = glm::normalize(glm::cross(axis_x, axis_y));
	coverAxis trans = coverAxis(toGlmVec3(startPos), axis_x, axis_y, axis_z);
	lastCover = trans;
	lastCover.tempAxis = temp_axis_z;
	MyMesh::FaceHandle nowFace = c->getFaceAtAngle(mesh, at_angle1);
	pointFromTwoSource* pair1 = c->searchAngle(at_angle1);
	//cout << "起始面為:" << nowFace.idx() << endl;
	glm::vec3 p1 = trans.cover(toGlmVec3(pair1[0].pointPos));
	glm::vec3 p2 = trans.cover(toGlmVec3(pair1[1].pointPos));
	//cout << "covered pair1[0].pointpos:(" << p1[0] << "," << p1[1] << "," << p1[2] << ")"<<endl;
	//cout << "covered pair1[1].pointpos:(" << p2[0] << "," << p2[1] << "," << p2[2] << ")" << endl;
	glm::vec3 pe = trans.cover(toGlmVec3(endPos));
	//cout << "covered endPoint:(" << pe[0] << "," << pe[1] << "," << pe[2] << ")";

	MyMesh::HalfedgeHandle mostPossible_he;
	OpenMesh::Vec3f most_interaction;
	float  most_avgX = -9999;
	for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(nowFace); fh_it != mesh.fh_end(nowFace); fh_it++){
		MyMesh::VertexHandle from = mesh.from_vertex_handle(*fh_it);
		MyMesh::VertexHandle to = mesh.to_vertex_handle(*fh_it);
		glm::vec3 point1 = trans.cover(toGlmVec3(mesh.point(from)));
		glm::vec3 point2 = trans.cover(toGlmVec3(mesh.point(to)));
		//cout << "起始邊:"<<from.idx()<<" (" << point1[0] << "," << point1[1] << "," << point1[2] << "); "<<to.idx()<<"(" << point2[0] << "," << point2[1] << "," << point2[2] << ");" << endl;
		bool isInteraction = (point1[1] >= 0 && point2[1] <= 0) || (point1[1] <= 0 && point2[1] >= 0);//如果有相交的話兩個點會有一個端點y值為正另一個端點y值為負
		if (isInteraction){
			glm::vec3 offset = point2 - point1;
			float time = (0 - point1[1]) / offset[1];
			glm::vec3 coverDomain_intersection = point1 + time*offset;
			float avgX = coverDomain_intersection[0];
			if (avgX>most_avgX){
				mostPossible_he = *fh_it;
				most_avgX = avgX;
				//glm::vec3 P1toP2 = point2 - point1;
				//float time = (0 - point1[1]) / P1toP2[1];
				OpenMesh::Vec3f interaction = mesh.point(from) + (mesh.point(to) - mesh.point(from))*time;
				most_interaction = interaction;
			}
		}
		//points.push_back(OpenMesh::Vec3f(interaction[0], interaction[1], interaction[2]));
	}
	//cout << "搜尋點開始:" << endl;
	points.push_back(most_interaction);
	MyMesh::HalfedgeHandle now_he = mesh.opposite_halfedge_handle(mostPossible_he);
	nowFace = mesh.opposite_face_handle(mostPossible_he);
	int endFaceIdx = c2->getFaceAtAngle(mesh, at_angle2).idx();
	const int max_iter_num = 100;
	//cout << "終點位於面:" << endFaceIdx << endl;
	//cout << "第一個面為:" << nowFace.idx()<<endl;
	int count = 0;
	while (nowFace.idx() != endFaceIdx){
		count++;
		if (count >= max_iter_num){//強制終止
			cout << "第一百次了!!!";
			//fail = true;
			/*vector<OpenMesh::Vec3f> failResult(2);
			failResult[0] = startPos;
			failResult[1] = endPos;
			return failResult;*/
			return points;
		}
		//cout << "處理面 " << nowFace.idx()<<":";
		for (MyMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(nowFace); fh_it != mesh.fh_end(nowFace); fh_it++){

			if ((*fh_it) != now_he){//忽略更新的來源邊
				//cout << "半邊:" << (*fh_it).idx() << "反向半邊為:" << mesh.opposite_halfedge_handle(*fh_it).idx() << "反向面為:" << mesh.opposite_face_handle(*fh_it);
				OpenMesh::Vec3f from = mesh.point(mesh.from_vertex_handle(*fh_it));
				OpenMesh::Vec3f to = mesh.point(mesh.to_vertex_handle(*fh_it));
				glm::vec3 point1 = trans.cover(toGlmVec3(from));
				//cout << "point1 idx" << mesh.from_vertex_handle(*fh_it).idx() << ":(" << point1[0] << "," << point1[1] << "," << point1[2] << ") ";
				glm::vec3 point2 = trans.cover(toGlmVec3(to));
				//cout << "point2 idx"<<mesh.to_vertex_handle(*fh_it).idx()<<":(" << point2[0] << "," << point2[1] << "," << point2[2] << ") ";
				if ((point1[1] >= 0 && point2[1] <= 0) || (point1[1] <= 0 && point2[1] >= 0)){//找到一個符合兩端的各在x軸一邊的邊
					//cout << "從半邊:" << (*fh_it).idx();
					glm::vec3 P1toP2 = point2 - point1;
					float time = (0 - point1[1]) / P1toP2[1];
					//glm::vec3 interaction = point1 + P1toP2*time;
					OpenMesh::Vec3f intersection = from + (to - from)*time;
					points.push_back(intersection);
					//cout << "找到點(" << intersection[0] << "," << intersection[1] << "," << intersection[2] << ")";

					now_he = mesh.opposite_halfedge_handle(*fh_it);//反向當前半邊
					//cout << "更新半邊到:" << now_he.idx() << endl;
					nowFace = mesh.face_handle(now_he);//推進面
					break;
				}
				//cout << ";     ";
			}
		}
		//cout << endl;
	}
	points.push_back(endPos);
	return points;
}
OpenMesh::Vec3f polarMap::getgdcPointAtBtw(int contour1_index, float at_angle1, float at_angle2, float percent){
	vector<OpenMesh::Vec3f> points = getgdcPointsBtw(contour1_index, at_angle1, at_angle2);
	float totalLength = 0;
	for (int i = 0; i < points.size() - 1; i++){
		totalLength += (points[i + 1] - points[i]).length();
	}
	float nowLength = 0;
	for (int i = 0; i < points.size() - 1; i++){
		float len = (points[i + 1] - points[i]).length();
		if ((nowLength + len) / totalLength >= percent){//點就在這個範圍內
			float sub = ((percent*totalLength) - nowLength) / len;
			return points[i] + (points[i + 1] - points[i])*sub;
		}
		else
		{
			nowLength += len;
		}
	}
	cout << "percent:" << percent << "搜尋不到點!" << endl;
	return points.back();
}
vector<vector<OpenMesh::Vec3f>> polarMap::getIntersections(vector<PointCurve> contours,vector<int> pointNum){
	vector<vector<OpenMesh::Vec3f>> *points = new vector<vector<OpenMesh::Vec3f>>();
	points->resize(contours.size());
	if (contours.size() != pointNum.size()){
		return *points;
	}

	for (int i = 0; i < contours.size(); i++){//第i条等高线
		//cout << "getIntersections 第" << i << "条等高线";
		int pnum = pointNum[i];
		float interval = 1.0 / pnum;
		int nowCIdx = 0;
		for (int p = 0; p < pnum; p++){//等高线上第p个等角度点
			float nowAngle = p*interval;
			//cout << "nowAngle:" << nowAngle<<"|";
			while (true){//每個搜尋,這個搜尋會持續到找到該點應該處於的間隔為止
				if (nowCIdx != contours[i].points.size() - 1){//如果還沒有達到最後一個點
					float lastAngle = contours[i].points[nowCIdx].angle;
					float nextAngle = contours[i].points[nowCIdx + 1].angle;
					//cout << "lastAngle:" << lastAngle << " nextAngle:" << nextAngle;
					if (lastAngle <= nowAngle && nextAngle > nowAngle){//找到了点
						OpenMesh::Vec3f toNext = contours[i].points[nowCIdx+1].pos - contours[i].points[nowCIdx].pos;//从上个点到下个点的向量
						float percentage = (nowAngle - lastAngle) / (nextAngle - lastAngle);
						//cout << "percentage:" << percentage << " ";
						OpenMesh::Vec3f now_toNext = toNext*percentage;
						points->operator[](i).push_back(contours[i].points[nowCIdx].pos + now_toNext);
						break;

					}
					else{
						nowCIdx++;
					}
				}
				else{
					float lastAngle = contours[i].points[nowCIdx].angle;
					float nextAngle = 1;
					OpenMesh::Vec3f toNext = contours[i].points[0].pos - contours[i].points[nowCIdx].pos;//从上个点到下个点的向量
					float percentage = (nowAngle - lastAngle) / (nextAngle - lastAngle);
					OpenMesh::Vec3f now_toNext = toNext*percentage;
					points->operator[](i).push_back(contours[i].points[nowCIdx].pos + now_toNext);
					break;
				}
			}
			//cout << "| "<<endl;
		}
		//cout << endl;
	}
	return *points;
}
float cal_t(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3){
	return (x3*x2 - x3*x1 - x2*x2 + x2*x1
		+ y3*y2 - y3*y1 - y2*y2 + y2*y1
		+ z3*z2 - z3*z1 - z2*z2 + z2*z1) /
		-(x3*x3 -2*x3*x2 +x2*x2
		+y3*y3 -2*y3*y2 +y2*y2
		+z3*z2 -2*z3*z2 +z2*z2);
}
pointAtContour polarMap::closetPointInNext(int contourIdx, float angle){
	contour nowContour = *(contours->operator[](contourIdx)); 
	pointAtContour closet;
	float closet_dist=1000000;
	OpenMesh::Vec3f nowPos = nowContour.getPosAtAngle(angle);
	for each(contour* next in nowContour.nextContours){
		for (int i = 0; i < next->pfts.size(); i++){
			if ((next->pfts[i].pointPos - nowPos).length()<closet_dist)
			{
				OpenMesh::Vec3f p1 = nowPos;
				OpenMesh::Vec3f p2 = next->pfts[i].pointPos;
				if (i - 1 > 0){//試試前一個點->當前點的邊
					OpenMesh::Vec3f p3 = next->pfts[i - 1].pointPos;
					float t = cal_t(p1[0],p1[1],p1[2], p2[0],p2[1],p2[2], p3[0],p3[1],p3[2]);
					if (t >= 0 && t <= 1){//找到極值點,這個問題的極值點一定是極小值點
						OpenMesh::Vec3f pe = (p3 - p2)*t + p2;
						closet_dist = (pe - p1).length();
						float angle = (next->pfts[i - 1].angleByCal - next->pfts[i].angleByCal)*t + next->pfts[i].angleByCal;
						closet = pointAtContour(next->id,angle, pe);
					}
				}
				if (i + 1 <= next->pfts.size())
				{
					OpenMesh::Vec3f p3 = next->pfts[i + 1].pointPos;
					float t = cal_t(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
					if (t >= 0 && t <= 1){//找到極值點,這個問題的極值點一定是極小值點
						OpenMesh::Vec3f pe = (p3 - p2)*t + p2;
						if ((pe - p1).length() < closet_dist){
							closet_dist = (pe - p1).length();
							float angle = (next->pfts[i + 1].angleByCal - next->pfts[i].angleByCal)*t + next->pfts[i].angleByCal;
							closet = pointAtContour(next->id , angle, pe);
						}
					}
				}
			}
		}
		return closet;
	}
}
float distance_perp2proj(OpenMesh::Vec3f n, OpenMesh::Vec3f v){//計算與"v投影于n的向量"合成v的向量的長度
	n = n.normalize();
	float dot_nv = n[0] * v[0] + n[1] * v[1] + n[2] * v[2];
	return (v - n*dot_nv).length();
}
polarMap_angleProjection::polarMap_angleProjection(MyMesh mesh, int poleIndex):polarMap(mesh,poleIndex){
}
vector<vector<int>> polarMap_angleProjection::calBaseAxis(){
	distFormVertex = new float[mesh.n_vertices()]{BIG_NUMBER};//初始化距離矩陣
	cout << "初始化calBaseAxis:";
	for (int i = 0; i < mesh.n_vertices(); i++){

		distFormVertex[i] = BIG_NUMBER;

	}
	distFormVertex[poleIdx] = 0;
	int* lastVertexIndex = new int[mesh.n_vertices()];//初始化上一個點索引
	for (int idx = 0; idx < mesh.n_vertices(); idx++){
		lastVertexIndex[idx] = idx;
	}
	vector<OpenMesh::VertexHandle> candiates;
	candiates.push_back(OpenMesh::VertexHandle(poleIdx));
	while (candiates.size() > 0){//計算所有點到原點(poleIdx)的距離和路徑,這裡不是測地距離而是通過頂點之間的直接連接來就算
		calRingVertex(candiates, distFormVertex, lastVertexIndex);
	}
	for (int i = 0; i < mesh.n_vertices(); i++){
		if (distFormVertex[i] > farestFromVertex){
			farestFromVertex = distFormVertex[i];
		}
	}
	vector<vector<int>> axises;
	vector<candiateAxis> candiate_axis;

	cout << "localExtremum_idx:";
	for (int i = 0; i < localExtremum_idx.size(); i++){
		cout << localExtremum_idx[i] << ",";
	}

	for each (int traget_vidx in localExtremum_idx)//對於每一個局部極值點
	{//組裝一條從其實點到局部極值點的路徑
		//cout << "局部極值點來自:";
		vector<int> now_axis;
		now_axis.push_back(traget_vidx);
		int next_vidx = traget_vidx;
		while (next_vidx != poleIdx)
		{
			//cout << next_vidx << "<=";
			next_vidx = lastVertexIndex[next_vidx];
			//now_axis.push_back(next_vidx);
			now_axis.insert(now_axis.begin(), next_vidx);
		}
		now_axis.insert(now_axis.begin(), poleIdx);
		float total_length = 0;
		for (int i = 0; i < now_axis.size() - 1; i++){
			float len = (mesh.point(OpenMesh::VertexHandle(now_axis[i + 1])) - mesh.point(OpenMesh::VertexHandle(now_axis[i]))).length();
			total_length += len;
		}
		//cout << endl;
		candiate_axis.push_back(candiateAxis(now_axis, total_length));
		axises.push_back(now_axis);
	}
	bool calContourAngles = true;
	if(calContourAngles){
		vector<contour*> wasteContour;
		for each(contour *c in *contours){//每一條等高線
			vector<OpenMesh::Vec3f> pointRing;
			for each (pointFromTwoSource p in c->pfts)
			{
				pointRing.push_back(p.pointPos);
			}

			PlaneFit plane(pointRing);
			//cout << "contour" << c->id << "size:" << pointRing.size() << "c->size():" << c->pfts.size() << endl;
			if (plane.is_fit){
				OpenMesh::Vec3f n = plane.normal();
				bool debug =false;
				vector<floatInt_pair> candiate_startPoint;
				if (debug)
					cout << "contour" << c->id << ":";
				for each(candiateAxis axis in candiate_axis){

					for (int i = 0; i < axis.vidx.size() - 1; i++){
						int index = c->indexOfTwoVidx(axis.vidx[i], axis.vidx[i + 1]);
						if (index >= 0)
						{
							if (debug)
								cout << "找到一個可能起始點index為:" << index;
							candiate_startPoint.push_back(floatInt_pair(axis.weight, index));
						}
					}
				}
				if (debug)
					cout << endl;
				floatInt_pair max_pair = candiate_startPoint[0];

				for (int i = 1; i < candiate_startPoint.size(); i++)
				{
					if (candiate_startPoint[i].weight > max_pair.weight){
						max_pair = candiate_startPoint[i];
					}
				}

				if (debug)
					cout << "contour 長度" << c->length << endl;
				//計算contour的平行于平面的長度,直接計算點之間的歐幾里得距離的長度不一樣
				float contourPlaneLength = 0;
				for (int i = 0; i < c->pfts.size() - 1; i++){
					OpenMesh::Vec3f v = c->pfts[i + 1].pointPos - c->pfts[i].pointPos;
					contourPlaneLength += distance_perp2proj(n, v);
				}
				OpenMesh::Vec3f v = c->pfts.back().pointPos - c->pfts[0].pointPos;
				contourPlaneLength += distance_perp2proj(n, v);
				int startIdx = max_pair.vidx;
				float now_angle = 0;
				c->pfts[startIdx].angleByCal = now_angle;
				//計算contour每一個點的角度
				for (int i = startIdx; i < c->pfts.size() - 1; i++){
					OpenMesh::Vec3f v = c->pfts[i + 1].pointPos - c->pfts[i].pointPos;
					now_angle += distance_perp2proj(n, v) / contourPlaneLength;
					if (debug)
						cout << i + 1 << ":" << "nowangle:" << now_angle << "| ";
					c->pfts[i + 1].angleByCal = now_angle;
				}
				if (debug)
					cout << endl << "到達contour尾端" << endl;
				now_angle += (distance_perp2proj(n, c->pfts[c->pfts.size() - 1].pointPos - c->pfts[0].pointPos) / contourPlaneLength);
				if (debug)
					cout << 0 << ":" << "nowangle:" << now_angle << "| ";
				c->pfts[0].angleByCal = now_angle;
				for (int i = 0; i < startIdx; i++){
					OpenMesh::Vec3f v = c->pfts[i + 1].pointPos - c->pfts[i].pointPos;
					now_angle += (distance_perp2proj(n, v) / contourPlaneLength);
					if (debug)
						cout << i + 1 << ":" << "dist:" << (c->pfts[i + 1].pointPos - c->pfts[i].pointPos).length() << "_nowangle:" << now_angle << "| ";
					c->pfts[i + 1].angleByCal = now_angle;
				}
				c->pfts[startIdx].angleByCal = 1;
			}
			else
			{
				cout << "警告!polarMap_angleProjection calBaseAxis contour" << c->id << "無法擬合平面!" << endl;
				//system("pause");
				bool debug = c->id == 9;
				vector<floatInt_pair> candiate_startPoint;
				if (debug)
					cout << "contour" << c->id << ":";
				cout << "candiate_axis size:" << candiate_axis.size();
				for each(candiateAxis axis in candiate_axis){

					for (int i = 0; i < axis.vidx.size() - 1; i++){
						int index = c->indexOfTwoVidx(axis.vidx[i], axis.vidx[i + 1]);
						if (index >= 0)
						{
							if (debug)
								cout << "找到一個可能起始點index為:" << index;
							candiate_startPoint.push_back(floatInt_pair(axis.weight, index));
						}
					}
				}
				if (debug)
					cout << "candiate_startPoint size:" << candiate_startPoint.size() << endl;
				floatInt_pair max_pair = candiate_startPoint[0];
				for (int i = 1; i < candiate_startPoint.size(); i++)
				{
					if (candiate_startPoint[i].weight > max_pair.weight){
						max_pair = candiate_startPoint[i];
					}
				}

				if (debug)
					cout << "contour 長度" << c->length << endl;
				int startIdx = max_pair.vidx;
				float now_angle = 0;
				c->pfts[startIdx].angleByCal = now_angle;
				//計算contour每一個點的角度
				for (int i = startIdx; i < c->pfts.size() - 1; i++){
					now_angle += ((c->pfts[i + 1].pointPos - c->pfts[i].pointPos).length() / c->length);
					if (debug)
						cout << i + 1 << ":" << "nowangle:" << now_angle << "| ";
					c->pfts[i + 1].angleByCal = now_angle;
				}
				if (debug)
					cout << endl << "到達contour尾端" << endl;
				now_angle += ((c->pfts[c->pfts.size() - 1].pointPos - c->pfts[0].pointPos).length() / c->length);
				if (debug)
					cout << 0 << ":" << "nowangle:" << now_angle << "| ";
				c->pfts[0].angleByCal = now_angle;
				for (int i = 0; i < startIdx; i++){
					now_angle += ((c->pfts[i + 1].pointPos - c->pfts[i].pointPos).length() / c->length);
					if (debug)
						cout << i + 1 << ":" << "dist:" << (c->pfts[i + 1].pointPos - c->pfts[i].pointPos).length() << "_nowangle:" << now_angle << "| ";
					c->pfts[i + 1].angleByCal = now_angle;
				}
				c->pfts[startIdx].angleByCal = 1;
				printf("pfts[startIdx].angleByCal:%e", c->pfts[startIdx].angleByCal);
				cout << endl;
				//system("pause");
			}
		}
	}
	/*
	cout << endl;
	for (int i = 0; i < (*contours).size();i++)
	{
	for (float angle = 0.1; angle <= 1; angle += 0.1){
	cout << "contour:" << (*contours)[i]->id<< " at angle:" << angle;
	getPointBetween(i, angle, angle);
	}
	}
	*/
	return axises;

}
void polarMap::initPathFinder(){
	vector<float> vertexList(mesh.n_vertices() * 3);
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		OpenMesh::VertexHandle v(i);
		OpenMesh::Vec3f pos = mesh.point(v);
		vertexList[i * 3 + 0] = pos[0];
		vertexList[i * 3 + 1] = pos[1];
		vertexList[i * 3 + 2] = pos[2];
	}
	vector<int> faceList(mesh.n_faces() * 3);
	for (int i = 0; i < mesh.n_faces(); i++){
		OpenMesh::FaceHandle f(i);
		int num = 0;
		for (MyMesh::FaceVertexCWIter fv_cwit = mesh.fv_cwbegin(f); fv_cwit != mesh.fv_cwend(f); fv_cwit++)
		{
			faceList[i * 3 + num++] = (*fv_cwit).idx();
		}
	}
	pFinder.initMesh(vertexList, faceList, contours->size());
	asynFinder.init(contours->size(), meshPath);
	for (int c = 0; c < contours->size(); c++){
		contour* now = contours->operator[](c);
		vector<int> vBoundary;
		for each(pointFromTwoSource pfts in now->pfts){
			if (!contain(vBoundary, pfts.sourceIdx_Small)){
				vBoundary.push_back(pfts.sourceIdx_Small);
			}
		}
		//cout << "contour small size:" << vBoundary.size();
		for each (contour* next in now->nextContours)//每一個下一條等高線也是邊界
		{
			vector<int> vBoundary_next;
			for each (pointFromTwoSource pfts in next->pfts)
			{
				if (!contain(vBoundary_next, pfts.sourceIdx_Large)){
					vBoundary_next.push_back(pfts.sourceIdx_Large);
				}
			}
			//cout << " contour large size:" << vBoundary_next.size();
			vBoundary.insert(vBoundary.begin(), vBoundary_next.begin(), vBoundary_next.end());
		}
		/*cout << "vBoundart size:" << vBoundary.size();
		for each (int vidx in vBoundary)
		{
			cout << vidx << ",";
		}
		cout << endl;*/
		pFinder.initSubFinder(c, vBoundary);
		asynFinder.initBoundaryRecord(c, vBoundary);
	}
	asynFinder.writeBoundaryData();
	//cout << "寫入boundary 完畢" << endl;
}
vector<OpenMesh::Vec3f> polarMap::getgdcPointsBtw(int contour1_index, float at_angle1, float at_angle2){
	
	contour *ct1 =contours->operator[](contour1_index);
	contour *ct2 = ct1->nextContours[0];
	OpenMesh::Vec3f ps = ct1->getPosAtAngle(at_angle1);
	int fs = ct1->getFaceAtAngle(mesh,at_angle1).idx();
	float pos_start[3]{ps[0], ps[1], ps[2]};
	OpenMesh::Vec3f pe = ct2->getPosAtAngle(at_angle2);
	int fe = ct2->getFaceAtAngle(mesh, at_angle2).idx();
	float pos_end[3]{pe[0], pe[1], pe[2]};
	vector<float*> path = pFinder.getPathInSubArea(contour1_index, pos_start, fs, pos_end, fe);
	//cout << "getgdcPointsBtw 找到最短路徑size:" << path.size()<<endl;
	vector<OpenMesh::Vec3f> result(path.size());
	for (int i = 0; i < path.size(); i++){
		result[i] = OpenMesh::Vec3f(path[i][0], path[i][1], path[i][2]);
		delete []path[i];
	}
	//cout << "轉換成vector后size:" << result.size();
	return result;
}
vector<OpenMesh::Vec3f> polarMap::getgdcPointsBtw(int contour1_index,int next_contour_index, float at_angle1, float at_angle2){

	contour *ct1 = contours->operator[](contour1_index);
	contour *ct2 = ct1->nextContours[next_contour_index];
	OpenMesh::Vec3f ps = ct1->getPosAtAngle(at_angle1);
	int fs = ct1->getFaceAtAngle(mesh, at_angle1).idx();
	float pos_start[3]{ps[0], ps[1], ps[2]};
	OpenMesh::Vec3f pe = ct2->getPosAtAngle(at_angle2);
	int fe = ct2->getFaceAtAngle(mesh, at_angle2).idx();
	float pos_end[3]{pe[0], pe[1], pe[2]};
	vector<float*> path = pFinder.getPathInSubArea(contour1_index, pos_start, fs, pos_end, fe);
	//cout << "getgdcPointsBtw 找到最短路徑size:" << path.size()<<endl;
	vector<OpenMesh::Vec3f> result(path.size());
	for (int i = 0; i < path.size(); i++){
		result[i] = OpenMesh::Vec3f(path[i][0], path[i][1], path[i][2]);
		delete[]path[i];
	}
	//cout << "轉換成vector后size:" << result.size();
	return result;
}
void polarMap::recordGdcPointLineBtw(vector<pathRequestPoint> requstLine){
	for each (pathRequestPoint data in requstLine)
	{
		contour *ct1 = contours->operator[](data.contour_index[0]);
		OpenMesh::Vec3f ps = ct1->getPosAtAngle(data.contour_angle[0]);
		int fs = ct1->getFaceAtAngle(mesh, data.contour_angle[0]).idx();
		
		contour *ct2 = contours->operator[](data.contour_index[1]);
		OpenMesh::Vec3f pe = ct2->getPosAtAngle(data.contour_angle[1]);
		int fe = ct2->getFaceAtAngle(mesh, data.contour_angle[1]).idx();
		asynFinder.requestPath(data.contour_index[0], fs, ps.data(), fe, pe.data());
	}
	if (requestTable.size() == 0){
		requestTable.push_back(requstLine);
	}
	else
	{//將requestLine插入到最後一次找到和其相同contour_index[0]的組後面
		int lastFind = -1;
		for (size_t i = 0; i < requestTable.size(); i++)
		{
			if (requestTable[i][0].contour_index[0] <= requstLine[0].contour_index[0]){
				lastFind = i;
			}
		}
		requestTable.insert(requestTable.begin() + lastFind + 1, requstLine);
	}
}
OpenMesh::Vec3f getPointAtPath(vector<OpenMesh::Vec3f> path,float percent){
	//cout << "getPointAtPath path.size():" << path.size() << " percent:" << percent << ";";

	if (path.size() <= 0){
		cout << "getPointAtPath 錯誤path size為0!";
		return OpenMesh::Vec3f(0, 0, 0);
	}
	if (percent < 0){
		return path.front();
	}
	if (percent > 1){
		return path.back();
	}
	float totalLen = 0;
	for (size_t i = 0; i < path.size()-1; i++)
	{
		totalLen += (path[i + 1] - path[i]).length();
	}
	float percentCount = 0;
	for (size_t i = 1; i < path.size(); i++)
	{
		float percent1 = percentCount;
		percentCount += (path[i] - path[i - 1]).length() /totalLen;
		//cout << "i=" << i << " percent:" << percentCount;
		if (percentCount >= percent){
			OpenMesh::Vec3f vec = path[i] - path[i - 1];
			return path[i-1]+ vec*((percent - percent1) / (percentCount - percent1));
		}
	}
	cout << "到結尾也沒有找到點!percentCount:" << percentCount << " tragetpercent:" << percent;
	return path.back();
}
vector<vector<pathRequestPoint>> polarMap::RespondsGdcPoints(){
	int groupIndex = 0;
	int lineIndex = 0;
	bool reachEnd = false;
	//cout << "印出requestTable sizes:";
	int totalRespond = 0;
	for (int i = 0; i < requestTable.size(); i++){
		//cout << i << ":" << requestTable[i].size() << ",";
		totalRespond += requestTable[i].size();
	}
	cout << "totalRespond:" << totalRespond;
	int respondCount = 0;
	do{
		reachEnd= asynFinder.writeToMissionFile(RESPOND_CYCLE);
		asynFinder.doIt();
		vector<vector<OpenMesh::Vec3f>> result= asynFinder.readPathResult();
		//cout << "實際收到" << result.size() << "條路徑";
		/*cout << "印出路徑:" << endl;
		for (int i = 0; i < result.size(); i++){
			cout << "path" << i << ":";
			for each (OpenMesh::Vec3f pt in result[i])
			{
				cout << "(" << pt[0] << "," << pt[1] << "," << pt[2] << "), ";
			}
			cout << endl;
		}*/
		int totalNum = result.size();
		respondCount += totalNum;
		int count = 0;
		int time = 0;

		while (totalNum>0 && groupIndex<requestTable.size())//totalNum沒有消耗完並且groupIndex沒有指向超越尾端
		{
			//cout << "第" << time++ << "次迭代:";
			if (totalNum < (requestTable[groupIndex].size() - lineIndex)){
				for (int i = 0; i < totalNum; i++){
					int rIdx = result.size() - (totalNum-i);
					float fix_y = 1 - requestTable[groupIndex][lineIndex + i].y;
					requestTable[groupIndex][lineIndex + i].pointPos = getPointAtPath(result[rIdx],fix_y);
					OpenMesh::Vec3f pos = requestTable[groupIndex][lineIndex + i].pointPos;
					//cout << "rIdx:" << rIdx << ">>(" <<pos[0]<<","<<pos[1]<<","<<pos[2]<<") " ;
					count++;
				}
				lineIndex += totalNum;
				//cout << "最後寫入" << totalNum << "點.";
				totalNum = 0;
			}
			else
			{
				//cout << "requestTable.size:" << requestTable.size()<<" result.size:"<<result.size();
				for (size_t i = lineIndex; i < requestTable[groupIndex].size(); i++)
				{
					int rIdx = result.size() - (totalNum - (i - lineIndex));
					//cout <<"i:"<<i<< " groupIndex:" << groupIndex << " lineIndex:" << lineIndex << "requestTable[groupIndex].size():" << requestTable[groupIndex].size()<<endl;
					//cout << "requestTable[groupIndex][i].y:" << requestTable[groupIndex][i].y;
					//cout << "result[rIdx].size():" << result[rIdx].size();
					float fix_y = 1 - requestTable[groupIndex][lineIndex + i].y;
					requestTable[groupIndex][i].pointPos = getPointAtPath(result[rIdx], fix_y);
					OpenMesh::Vec3f pos = requestTable[groupIndex][i].pointPos;
					//cout << "rIdx:" << rIdx << ">>(" << pos[0] << "," << pos[1] << "," << pos[2] << ") ";
					//cout << endl;
					count++;
				}
				//cout << "寫入" << requestTable[groupIndex].size() - lineIndex << "點,";
				totalNum -= requestTable[groupIndex].size() - lineIndex;
				groupIndex++;
				lineIndex = 0;
				
			}
		}
		//cout << "實際寫入" << count << "個點" << endl;
		//cout << "完成度" << respondCount << "/" << totalRespond << endl;
		//cout << "groupIndex:" << groupIndex << "lineIndex:" << lineIndex;
		/*int convNum = 0;
		for (int i = 0; i < groupIndex; i++){
			convNum += requestTable[i].size();
		}
		convNum += lineIndex;
		cout << "換算后的數量:" << convNum<<endl;*/
	} while (!reachEnd);
    /*if (groupIndex == requestTable.size() - 1 && lineIndex == requestTable.back().size() - 1)
		cout << "完美完成回應!";
	else
		cout << "fuck! groupIndex:"<<groupIndex<<" requestTable.size():"<<requestTable.size()<<" line Index:"<<lineIndex<<"  requestTable.back().size():"<<requestTable.back().size();
	*/
	return requestTable;
}
vector<vector<Vec3f_angle>> combineSample(vector<contour*> contourUnion, int sampleNum)
{
	float totalLen = 0;
	vector<vector<Vec3f_angle>> results;
	for each(contour* contour in contourUnion)
	{
		totalLen += contour->length;
	}
	for each (contour* contour in contourUnion)
	{
		int num = sampleNum*(contour->length / totalLen);
		float decimal = (float)sampleNum*(contour->length / totalLen) - (int)sampleNum*(contour->length / totalLen);
		if (decimal > 0.5f){//四捨五入
			num += 1;
		}
		vector<Vec3f_angle> ring(num);
		Vec3f_angle* result = contour->sample(num);
		for (int i = 0; i < num; i++){
			ring[i] = result[i];
		}
		delete []result;
		results.push_back(ring);
	}
	return results;
}
vector<vector<OpenMesh::Vec3f>> combineSample_withoutAngle(vector<contour*> contourUnion, int sampleNum){
	float totalLen = 0;
	vector<vector<OpenMesh::Vec3f>> results;
	for each(contour* contour in contourUnion)
	{
		totalLen += contour->length;
	}
	for each (contour* contour in contourUnion)
	{
		int num = sampleNum*(contour->length / totalLen);
		float decimal = (float)sampleNum*(contour->length / totalLen) - (int)sampleNum*(contour->length / totalLen);
		if (decimal > 0.5f){//四捨五入
			num += 1;
		}
		vector<OpenMesh::Vec3f> result = contour->sample_withoutAngle(num);
		results.push_back(result);
	}
	return results;
}
