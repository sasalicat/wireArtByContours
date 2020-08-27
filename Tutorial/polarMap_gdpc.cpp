#include "stdafx.h"
#include "polarMap_gdpc.h"

omMesh normalizeMesh(omMesh input){
	DGPC::Vector3<double> meshCenter = DGPC::Vector3<double>(0, 0, 0);

	//遍历所有模型顶点加总获得平均值
	for (int i = 0; i < input.n_vertices(); i++){
		omMesh::VertexHandle vh = omMesh::VertexHandle(i);
		DGPC::Vector3<double> point = input.point(vh);
		meshCenter = meshCenter+point;
	}
	meshCenter = DGPC::Vector3<double>(meshCenter.x() / input.n_vertices(), meshCenter.y() / input.n_vertices(),meshCenter.z()/input.n_vertices());
	cout << "avg  mesh center:(" << meshCenter[0] << "," << meshCenter[1] << "," << meshCenter[2] << ")" << endl;
	//更新所有的顶点到相对于中心点的位置
	float maxLength = 0;
	for (int i = 0; i < input.n_vertices(); i++){
		omMesh::VertexHandle vh = omMesh::VertexHandle(i);
		DGPC::Vector3<double> point = input.point(vh);
		//MyMesh::Point point = input.point(vh);
		float len = (point - meshCenter).length();
		if (len > maxLength)
		{
			maxLength = len;
		}
	}
	for (int i = 0; i < input.n_vertices(); i++){
		omMesh::VertexHandle vh = omMesh::VertexHandle(i);
		DGPC::Vector3<double> offset = input.point(vh)-meshCenter;
		
		input.set_point(vh, DGPC::Vector3<double>(offset.x()/maxLength,offset.y()/maxLength,offset.z()/maxLength));
		//input.set_point(vh, (point - meshCenter) / maxLength);
	}

	input.update_normals();
	return input;
}

polarMap_gdpc::polarMap_gdpc(string path, int polarIdx) :polarMap(MyMesh(),polarIdx){
	mesh_om.openOBJ(path.c_str());
	mesh_om = normalizeMesh(mesh_om);
	//OpenMesh::ArrayKernel mesh = mesh_om;
	//cout << "mesh_om n_faces():" << mesh.n_faces();
}

polarMap_gdpc::~polarMap_gdpc()
{
}
void polarMap_gdpc::findExtremumVertexs(){
	localExtremum_idx.clear();
	cout << "印出localExtremum開始---------------------------------------------------" << endl;

	for (int i = 0; i<mesh_om.n_vertices(); i++){
		//cout << "i=" << i<<":";
		MyMesh::VertexHandle vh = MyMesh::VertexHandle(i);
		bool isExtremum = true;//極值點的旗標
		for (MyMesh::VertexVertexCWIter vv_cwit = mesh_om.vv_cwbegin(vh); vv_cwit != mesh_om.vv_cwend(vh); vv_cwit++){
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
			cout << "U[" << i << "]:" << U[i] << "=>";
			for (MyMesh::VertexVertexCWIter vv_cwit = mesh_om.vv_cwbegin(vh); vv_cwit != mesh_om.vv_cwend(vh); vv_cwit++){
				int vidx = (*vv_cwit).idx();
				cout << "idx " << vidx << ":" << U[vidx];
			}
			cout << endl;
			
		}
	}
	cout << "印出localExtremum結束---------------------------------------------------"<<endl;
}
int polarMap_gdpc::pointNum(){
	return mesh_om.n_vertices();
}
void polarMap_gdpc::createMap(){
	generator = new DGPCgenerator(mesh_om);
	generator->setNodeSource(poleIdx);
	endIdx = generator->run();
	cout << "polar map gdpcЧΘ!"<<endl;
	int vnum= pointNum();
	U = new float[vnum];
	theta = new float[vnum];
	for (int i = 0; i <vnum; i++){
		polarPoint pointi = getPointFrom(i);
		U[i] = pointi.distance;
		theta[i] = pointi.angle;
	}
	findExtremumVertexs();
}
polarPoint polarMap_gdpc::getPointFrom(int vidx){
	//cout << "gdpc get point from;";
	Point p = mesh_om.point(omMesh::VertexHandle(vidx));
	OpenMesh::Vec3f vf = OpenMesh::Vec3f(p[0],p[1],p[2]);
	return polarPoint(vidx,vf,(float)generator->getDistance(vidx),(float)generator->getAngle(vidx));
}
void polarMap_gdpc::drawMeshToImage(){
	cout << "gdpc drawMeshToImage" << endl;
	painter = new draw2dConnectMap(1001, 1001);
	//璸衡border
	OpenMesh::Vec2f *euclidList = new OpenMesh::Vec2f[mesh_om.n_vertices()];
	float border[4] = { 0, 0, 0, 0 };
	for (int v = 0; v < mesh_om.n_vertices(); v++){
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
	//睰フ
	float center[2] = { (border[0] + border[2]) / 2, (border[1] + border[3]) / 2 };
	float radius[2] = { border[2] - center[0], border[3] - center[1] };
	radius[0] *= 1.1;//糴籹硑0.1フ
	radius[1] *= 1.1;//蔼籹硑0.1フ
	border[0] = center[0] - radius[0];//эborder
	border[1] = center[1] - radius[1];
	border[2] = center[0] + radius[0];
	border[3] = center[1] + radius[1];
	for (int f = 0; f < mesh_om.n_faces(); f++){
		OpenMesh::FaceHandle fh(f);
		vector<float> list;
		bool containEP = false;
		for (MyMesh::FaceVertexIter fv_it = mesh_om.fv_begin(fh); fv_it != mesh_om.fv_end(fh); fv_it++){
			if ((*fv_it).idx() == endIdx)
			{
				containEP = true;
				break;
			}
			//OpenMesh::Vec2f pos= getPointFrom((*fv_it).idx()).toEuclid2d();
			OpenMesh::Vec2f pos = euclidList[(*fv_it).idx()];
			pos = normalizePos(border, pos);//盢翴x,y耴て01ぇ丁
			list.push_back(pos[0]);
			list.push_back(pos[1]);
		}
		if (!containEP){
			float* ldata = list.data();
			painter->drawTriangle(threePoints(ldata));
		}
	}
	//OpenMesh::Vec2f traget_point = normalizePos(border, euclidList[zero_anix_vertexIdx]);
	//painter->drawPoint(traget_point[0], traget_point[1]);
	painter->showMap();
	delete euclidList;
}
int polarMap_gdpc::getFasthestIdx(){
	return endIdx;
}