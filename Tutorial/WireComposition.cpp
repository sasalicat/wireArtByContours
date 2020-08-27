#include "stdafx.h"
#include "WireComposition.h"

WireComposition::WireComposition(GraphArray *graph,vector<int> *startNodes,vector<int> *endNodes)
{
	/*this->graph = graph;
	mappings.resize(graph->array.size());
	//解析graph,把以node為主的GraphArray分解為以path為主
	for (size_t i = 0; i < graph->array.size(); i++)
	{
		for each (pathNode *node in graph->array[i])
		{
			pathMapping newMapping;
			for (size_t n = 0; n < node->connects.size(); n++)
			{
				pathNode* next = node->connects[n];
				
				newMapping.intersectionNode = next;
				path *p_in= new path();
				p_in->start_2d[0] = node->x_2d;
				p_in->start_2d[1] = node->y_2d;
				p_in->end_2d[0] = next->x_2d;
				p_in->end_2d[1] = next->y_2d;
				p_in->point3ds.push_back(node->position);
				for each (OpenMesh::Vec3f point in node->path[n])
				{
					p_in->point3ds.push_back(point);
				}
				p_in->point3ds.push_back(next->position);
				newMapping.path_in = p_in;
				int p_length = p_in->point3ds.size();
				OpenMesh::Vec3f vec1 = p_in->point3ds[p_length-1]-p_in->point3ds[p_length-2];
				for (int l = 0; l < next->connects.size(); l++){
					
					pathNode *nextnext = next->connects[l];
					if (nextnext!=node)
					{
						path *p_out = new path();
						p_out->start_2d[0] = next->x_2d;
						p_out->start_2d[1] = next->y_2d;
						p_out->end_2d[0] = nextnext->x_2d;
						p_out->end_2d[1] = nextnext->y_2d;
						p_out->point3ds.push_back(next->position);
						for each (OpenMesh::Vec3f point in next->path[l])
						{
							p_out->point3ds.push_back(point);
						}
						p_out->point3ds.push_back(nextnext->position);

						int p_length = p_out->point3ds.size();
						newMapping.path_out.push_back(p_out);
						OpenMesh::Vec3f vec2 = p_out->point3ds[1] - p_out->point3ds[0];
						float w= OpenMesh::dot(vec1.normalize(), vec2.normalize());
						newMapping.weight.push_back(w);
					}
				
				}
			}
			int index = floorf(newMapping.intersectionNode->y_2d);
			
			mappings[index].push_back(newMapping);
		}
	}*/
	/*NodeWeights.resize(graph->nodeNum);
	for each (vector<pathNode*> level in graph->array)
	{
		for each (pathNode* node in level)
		{
			vector<float> weights;
			for (size_t i = 0; i < node->connects.size(); i++)//i為來方向的索引值
			{
				for (size_t j = 0; j < node->connects.size(); j++)//j為去方向的索引值
				{
					if (i!=j){
						OpenMesh::Vec3f vec_from =node->position- node->path[i].back();
						OpenMesh::Vec3f vec_to =node->path[j].front()- node->position;
						float weight = OpenMesh::dot(vec_from, vec_to);
						weights.push_back(weight);
					}
				}
			}
		}
	}*/
	this->startNodeids = startNodes;
	this->endNodeids = endNodes;
	for (size_t i = 0; i < graph->array.size(); i++)
	{
		for (size_t j = 0; j < graph->array[i].size(); j++)
		{
			pathNode *node = graph->array[i][j];
			if (!contain(*startNodeids, node->id) && !contain(*endNodeids, node->id)){
				cout << "node:" << node->id;
				nodeMapping newMapping;
				newMapping.init(node);
				if (mappings.size() == 5){
					cout << "創造第5個mapping!";
				}
				mappings.push_back(newMapping);
			}
		}
	}
}
bool WireComposition::Solve(){
	IloEnv env;
	printf("IloEnv!");
	IloExpr expressions(env);
	IloRangeArray constraints(env);
	for (size_t i = 0; i < mappings.size(); i++)//i代表node,每個node迭代一次
	{
		mappings[i].initOptValues(env);
		int upPointNum = mappings[i].upNodeIndex.size();
		int lowPointNum = mappings[i].lowNodeIndex.size();

			for (size_t j = 0; j < mappings[i].upNodeIndex.size(); j++)//j代表node的上級點,每個node的每個上級點迭代一次
			{
				for (size_t k = 0; k < mappings[i].lowNodeIndex.size(); k++)//k代表上級點與下級點的對應關係,在每個上級點的基礎下,每個下級點對應迭代一次
				{
					expressions += mappings[i].weights[j][k] * mappings[i].optValues[j]->operator[](k);
				}
			}
		
		//cout << "expression:" << expressions;

		if (upPointNum == lowPointNum)
		{
			int halfPathNum = mappings[i].weights.size();
			for (size_t j = 0; j < halfPathNum; j++)//j代表node的上級點,每個node的每個上級點迭代一次
			{
				IloExpr sum_up_choose(env);//每個上級點必須且只能選擇一個下級點
				//cout << "j=" << j << " mappings[i].optValues[j].size():" << mappings[i].optValues[j]->getSize();
				for (size_t k = 0; k < halfPathNum; k++)//k代表上級點與下級點的對應關係,在每個上級點的基礎下,每個下級點對應迭代一次
				{
					//cout << "k=" << k;
					//IloNumVarArray *array = mappings[i].optValues[j];
					sum_up_choose += mappings[i].optValues[j]->operator[](k);
					//cout << "add";
				}
				//cout << "sum_up_choose:" << sum_up_choose<<endl;
				constraints.add(sum_up_choose == 1);
			}
			
			for (size_t j = 0; j < halfPathNum; j++)//j代表node的上級點,每個node的每個上級點迭代一次
			{
				IloExpr sum_low_bechoosen(env);//且每個下級點只能被一個上級點選擇
				for (size_t k = 0; k < halfPathNum; k++)//k代表上級點與下級點的對應關係,在每個上級點的基礎下,每個下級點對應迭代一次
				{
					sum_low_bechoosen += mappings[i].optValues[k]->operator[](j);
				}
				//cout << "sum_low_bechoosen:" << sum_low_bechoosen << endl;
				constraints.add(sum_low_bechoosen == 1);
			}
		}
		else if (upPointNum>lowPointNum)
		{
			for (int j = 0; j < upPointNum; j++){
				IloExpr sum_up_choose(env);
				for (int k = 0; k < lowPointNum; k++){
					sum_up_choose += mappings[i].optValues[j]->operator[](k);
				}
				constraints.add(sum_up_choose <= 1);//上級點不一定要選擇下級點
			}
			for (size_t k = 0; k < lowPointNum; k++)
			{
				IloExpr sum_low_bechoosen(env);
				for (size_t j = 0; j < upPointNum; j++)
				{
					sum_low_bechoosen += mappings[i].optValues[j]->operator[](k);
				}
				constraints.add(sum_low_bechoosen == 1);//下級點必須要被選中一次
			}
		}
		else//如果下級點比上級點多
		{
			for (int j = 0; j < upPointNum; j++)
			{
				IloExpr sum_up_choose(env);
				for (int k = 0; k < lowPointNum; k++){
					sum_up_choose += mappings[i].optValues[j]->operator[](k);
				}
				constraints.add(sum_up_choose == 1);//上級點一定要選擇一個下級點
			}
			for (size_t k = 0; k < lowPointNum; k++)
			{
				IloExpr sum_low_bechoosen(env);
				for (size_t j = 0; j < upPointNum; j++)
				{
					sum_low_bechoosen += mappings[i].optValues[j]->operator[](k);
				}
				constraints.add(sum_low_bechoosen <= 1);//下級點不一定要被選中
			}
		}
		cout << "mapping[" << i << "]結束.";
	}
	IloModel model(env);
	model.add(IloMinimize(env, expressions));
	model.add(constraints);

	IloCplex cplex(model);
	cplex.setOut(env.getNullStream());
	cout << "開始solve:";
	if (!cplex.solve()){
		cout << "Fail to optimize the model." << endl;
		return false;
	}
	else
	{
		cout << "optimize success!" << endl;
	}
	for (int i = 0; i < mappings.size(); i++){
		cout << "第i="<<i<<"個mapping:";
		int upPointNum = mappings[i].upNodeIndex.size();
		int lowPointNum = mappings[i].lowNodeIndex.size();
		vector<pathNode*> hiddenNextPoints;
		if (upPointNum < lowPointNum){
			for each (int index in  mappings[i].lowNodeIndex)
			{
				hiddenNextPoints.push_back(mappings[i].centerNode->connects[index]);
			}
		}
		for (int j = 0; j < mappings[i].optValues.size(); j++)
		{
			int index1 = mappings[i].upNodeIndex[j];
			cout << "第j=" << j << "個 index1:"<<index1<<"upNode result:{";
			IloNumVarArray *varArray = mappings[i].optValues[j];
			IloNumArray result(env);
			cplex.getValues(*varArray, result);

			for (int k = 0; k < result.getSize(); k++){
				cout << "weight:" << mappings[i].weights[j][k] << ">>" << result[k] << ", ";
				if (result[k]==1)
				{
					int index2 = mappings[i].lowNodeIndex[k];
					cout << "index2:" << index2;
					mappings[i].centerNode->directedConnects[mappings[i].centerNode->connects[index1]] = mappings[i].centerNode->connects[index2];
					for (int h = 0; h < hiddenNextPoints.size(); h++)
					{
						if (hiddenNextPoints[h] == mappings[i].centerNode->connects[index2]){
							hiddenNextPoints.erase(hiddenNextPoints.begin() + h);
							break;
						}
					}
				}
			}
			cout << "}";
		}
		cout << "印出node:" << mappings[i].centerNode->id << " directedConnects:";
		for each(pair<pathNode*, pathNode*>pair in mappings[i].centerNode->directedConnects){
			cout << pair.first << ">>" << pair.second << "; ";
		}
		cout << endl;
		if (upPointNum < lowPointNum){
			for each (pathNode* node in hiddenNextPoints)
			{
				specialStartNodes.push_back(pair<pathNode*, pathNode*>(mappings[i].centerNode, node));
			}
		}
	}
	return true;
}

void WireComposition::SolveOneByOne(){
	IloEnv env;
	for (int i = 0; i < mappings.size(); i++)
	{
		cout << "i=" << i << "開始處理";
		mappings[i].initOptValues(env);
		int upPointNum = mappings[i].upNodeIndex.size();
		int lowPointNum = mappings[i].lowNodeIndex.size();
		IloExpr expressions(env);
		for (size_t j = 0; j < mappings[i].upNodeIndex.size(); j++)//j代表node的上級點,每個node的每個上級點迭代一次
		{
			for (size_t k = 0; k < mappings[i].lowNodeIndex.size(); k++)//k代表上級點與下級點的對應關係,在每個上級點的基礎下,每個下級點對應迭代一次
			{
				expressions += mappings[i].weights[j][k] * mappings[i].optValues[j]->operator[](k);
			}
		}
		IloRangeArray constraints(env);
		if (upPointNum == lowPointNum)
		{
			int halfPathNum = mappings[i].weights.size();
			for (size_t j = 0; j < halfPathNum; j++)//j代表node的上級點,每個node的每個上級點迭代一次
			{
				IloExpr sum_up_choose(env);//每個上級點必須且只能選擇一個下級點
				//cout << "j=" << j << " mappings[i].optValues[j].size():" << mappings[i].optValues[j]->getSize();
				for (size_t k = 0; k < halfPathNum; k++)//k代表上級點與下級點的對應關係,在每個上級點的基礎下,每個下級點對應迭代一次
				{
					//cout << "k=" << k;
					//IloNumVarArray *array = mappings[i].optValues[j];
					sum_up_choose += mappings[i].optValues[j]->operator[](k);
					//cout << "add";
				}
				//cout << "sum_up_choose:" << sum_up_choose<<endl;
				constraints.add(sum_up_choose == 1);
			}

			for (size_t j = 0; j < halfPathNum; j++)//j代表node的上級點,每個node的每個上級點迭代一次
			{
				IloExpr sum_low_bechoosen(env);//且每個下級點只能被一個上級點選擇
				for (size_t k = 0; k < halfPathNum; k++)//k代表上級點與下級點的對應關係,在每個上級點的基礎下,每個下級點對應迭代一次
				{
					sum_low_bechoosen += mappings[i].optValues[k]->operator[](j);
				}
				//cout << "sum_low_bechoosen:" << sum_low_bechoosen << endl;
				constraints.add(sum_low_bechoosen == 1);
			}
		}
		else if (upPointNum>lowPointNum)
		{
			for (int j = 0; j < upPointNum; j++){
				IloExpr sum_up_choose(env);
				for (int k = 0; k < lowPointNum; k++){
					sum_up_choose += mappings[i].optValues[j]->operator[](k);
				}
				constraints.add(sum_up_choose <= 1);//上級點不一定要選擇下級點
			}
			for (size_t k = 0; k < lowPointNum; k++)
			{
				IloExpr sum_low_bechoosen(env);
				for (size_t j = 0; j < upPointNum; j++)
				{
					sum_low_bechoosen += mappings[i].optValues[j]->operator[](k);
				}
				constraints.add(sum_low_bechoosen == 1);//下級點必須要被選中一次
			}
		}
		else//如果下級點比上級點多
		{
			for (int j = 0; j < upPointNum; j++)
			{
				IloExpr sum_up_choose(env);
				for (int k = 0; k < lowPointNum; k++){
					sum_up_choose += mappings[i].optValues[j]->operator[](k);
				}
				constraints.add(sum_up_choose == 1);//上級點一定要選擇一個下級點
			}
			for (size_t k = 0; k < lowPointNum; k++)
			{
				IloExpr sum_low_bechoosen(env);
				for (size_t j = 0; j < upPointNum; j++)
				{
					sum_low_bechoosen += mappings[i].optValues[j]->operator[](k);
				}
				constraints.add(sum_low_bechoosen <= 1);//下級點不一定要被選中
			}
		}

		IloModel model(env);
		model.add(IloMinimize(env, expressions));
		model.add(constraints);

		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		if (!cplex.solve()){
			cout << "mapping[" <<i<<"] solve失敗.";
			solveFailIndexRecord.push_back(i);
			cout << "upPointNum:" << upPointNum << "lowPointNum:" << lowPointNum;
			cout << "expression:" << expressions << endl;
			cout << "constraints:" << constraints;
			cout << endl;
		}
		else{
			bool debug = false;
			if (upPointNum > lowPointNum){
				cout << "upPointNum:"<<upPointNum<<" > lowPointNum:"<<lowPointNum;
				debug = true;
			}
			else if (upPointNum<lowPointNum)
			{
				cout << "upPointNum:"<<upPointNum<<"< lowPointNum:"<<lowPointNum;
				debug = true;
			}
			vector<pathNode*> hiddenNextPoints;
			if (upPointNum < lowPointNum){
				for each (int index in  mappings[i].lowNodeIndex)
				{
					hiddenNextPoints.push_back(mappings[i].centerNode->connects[index]);
				}
			}
			for (int j = 0; j < mappings[i].optValues.size(); j++)
			{
				int index1 = mappings[i].upNodeIndex[j];
				//cout << "第j=" << j << "個 index1:" << index1 << "upNode result:{";
				IloNumVarArray *varArray = mappings[i].optValues[j];
				IloNumArray result(env);
				cplex.getValues(*varArray, result);
				
				for (int k = 0; k < result.getSize(); k++){
					if (debug)
						cout << "weight:" << mappings[i].weights[j][k] << ">>" << result[k] << ", ";
					if (result[k] == 1)
					{
						int index2 = mappings[i].lowNodeIndex[k];
						//cout << "index2:" << index2;
						mappings[i].centerNode->directedConnects[mappings[i].centerNode->connects[index1]] = mappings[i].centerNode->connects[index2];
						for (int h = 0; h < hiddenNextPoints.size();h++)
						{
							if (hiddenNextPoints[h] == mappings[i].centerNode->connects[index2]){
								hiddenNextPoints.erase(hiddenNextPoints.begin()+h);
								break;
							}
						}
					}
				}

				//cout << "}";
			}
			if (upPointNum < lowPointNum){
				for each (pathNode* node in hiddenNextPoints)
				{
					specialStartNodes.push_back(pair<pathNode*, pathNode*>(mappings[i].centerNode, node));
				}
			}
			if (debug)
			{
				cout << endl;
			}
		}
	}
}
WireComposition::~WireComposition()
{
}

void WireComposition::test(){
	printf("WireComposition 開始!");
	IloEnv env;
	printf("IloEnv!");
	IloExpr expressions(env);
	IloRangeArray constraints(env);
	IloNumVarArray x(env);
	for (int i = 0; i < 100000; i++){
		x.add(IloNumVar(env));
		expressions += (x[i]) *(x[i]);
	}
	IloModel model(env);
	model.add(IloMinimize(env, expressions));
	model.add(constraints);
	IloCplex cplex(model);
	printf("setOut開始");
	cplex.setOut(env.getNullStream());
	printf("solve開始");
	bool result = cplex.solve();
	printf("solve結果:%d", result);
	IloNumArray result_x(env);
	cplex.getValues(x, result_x);
	for (size_t i = 0; i < 1000; i++)
	{
		printf("i=%d x:%f",i,result_x[i]);
	}
	
}
MyMTSP::MyMTSP(){
}
struct pixelPathMarker
{
	int pixelNode1[2];
	int pixelNode2[2];
	vector<OpenMesh::Vec3f> path3d;
	pixelPathMarker(graphLine* line){
		pixelNode1[0] = line->points->front().x;
		pixelNode1[1] = line->points->front().y;
		pixelNode2[0] = line->points->back().x;
		pixelNode2[1] = line->points->back().y;
		for (size_t i = 0; i < line->points->size(); i++)
		{
			pixelPoint pt= line->points->operator[](i);
			path3d.push_back(OpenMesh::Vec3f(pt.percentage_x, pt.percentage_y, 0));
		}
	}
};
bool MyMTSP::Init(controlGraph pixelGraph){
	vector<pixelPathMarker> paths;
	cout << "ngraph node 數量:" << pixelGraph.mainGraphs.size();
	for (size_t i = 0; i < pixelGraph.mainGraphs.size(); i++)
	{
		paths.push_back(pixelPathMarker(pixelGraph.mainGraphs[i]));
		if (pixelGraph.mainGraphs[i]->points->front().percentage_y==0)
		{
			dominant_startIdx.push_back(1 + i * 2);
		}
		else if (pixelGraph.mainGraphs[i]->points->front().percentage_y == 1)
		{
			dominant_endIdx.push_back(1 + i * 2);
		}
		if (pixelGraph.mainGraphs[i]->points->back().percentage_y == 0)
		{
			dominant_startIdx.push_back(1 + i * 2+1);
		}
		else if (pixelGraph.mainGraphs[i]->points->back().percentage_y ==1)
		{
			dominant_endIdx.push_back(1 + i * 2 + 1);
		}
	}
	cout << "印出dominant_startIdx:";
	for each (int idx in dominant_startIdx)
	{
		cout << idx << ",";
	}
	cout << endl;
	cout << "印出dominant_endIdx:";
	for each (int idx in dominant_endIdx)
	{
		cout << idx << ",";
	}
	cout << endl;
	cout << "總共" << paths.size() << "條path.";
	weight.clear();
	segments_endpoints_idx.clear();
	connectness.clear();
	solvable = false;
	for each (pixelPathMarker path in paths)
	{
		if (path.path3d.size() < 2){
			printf("Not a line segment.\n");
			return solvable;
		}
	}
	num_node = 2 * paths.size() + 1;
	segments_endpoints_idx.resize(paths.size());
	for (int i = 0; i < paths.size(); ++i){
		segments_endpoints_idx[i].resize(2);
		segments_endpoints_idx[i][0] = 2 * i + 1;
		segments_endpoints_idx[i][1] = 2 * i + 2;
	}
	std::vector<bool> connectness_element(num_node, true);
	connectness.resize(num_node, connectness_element);
	for (int i = 0; i < paths.size(); ++i){
		//int front_idx = paths[i].node1->id;
		//int back_idx = paths[i].node2->id;
		//int endpoint_idx_current[2] = { front_idx, back_idx };
		int posId_current[2][2] = { { paths[i].pixelNode1[0], paths[i].pixelNode1[1] }, { paths[i].pixelNode2[0], paths[i].pixelNode2[1] } };
		
		for (int j = i + 1; j < paths.size(); ++j){
			//front_idx = paths[j].node1->id;
			//back_idx = paths[j].node2->id;
			//int endpoint_idx_compare[2] = { front_idx, back_idx };
			int posId_compare[2][2] = { { paths[j].pixelNode1[0], paths[j].pixelNode1[1] }, { paths[j].pixelNode2[0], paths[j].pixelNode2[1] } };
			for (int m = 0; m < 2; ++m){
				for (int n = 0; n < 2; ++n){
					if (posId_current[m][0] != posId_compare[n][0] || posId_current[m][1] != posId_compare[n][1]){
						connectness[segments_endpoints_idx[i][m]][segments_endpoints_idx[j][n]] = false;
						connectness[segments_endpoints_idx[j][n]][segments_endpoints_idx[i][m]] = false;
					}
				}
			}
		}
	}
	double u = 1;
	std::vector<double> weight_element(num_node, 0);
	weight.resize(num_node, weight_element);
	for (int i = 0; i < paths.size(); ++i){
		int n_point = paths[i].path3d.size();
		OpenMesh::Vec3f point_front1 = paths[i].path3d[0];
		OpenMesh::Vec3f point_front2 = paths[i].path3d[1];
		OpenMesh::Vec3f point_back1 = paths[i].path3d[n_point - 1];
		OpenMesh::Vec3f point_back2 = paths[i].path3d[n_point - 2];
		OpenMesh::Vec3f endpoint_current[2] = { point_front1, point_back1 };
		OpenMesh::Vec3f tangentangle_current[2] = { (point_front1 - point_front2).normalize(), (point_back1 - point_back2).normalize() };
		for (int j = i + 1; j < paths.size(); ++j){
			n_point = paths[j].path3d.size();
			point_front1 = paths[j].path3d[0];
			point_front2 = paths[j].path3d[1];
			point_back1 = paths[j].path3d[n_point - 1];
			point_back2 = paths[j].path3d[n_point - 2];
			OpenMesh::Vec3f endpoint_compare[2] = { point_front1, point_back1 };
			OpenMesh::Vec3f tangentangle_compare[2] = { (point_front1 - point_front2).normalize(), (point_back1 - point_back2).normalize() };
			for (int m = 0; m < 2; ++m){
				for (int n = 0; n < 2; ++n){
					double distance = (endpoint_current[m] - endpoint_compare[n]).length();
					double w = distance + u * (1.0 + OpenMesh::dot(tangentangle_current[m],tangentangle_compare[n]))*0.5;
					//double w = u * (1.0 + tangentangle_current[m] * tangentangle_compare[n])*0.5;
					weight[segments_endpoints_idx[i][m]][segments_endpoints_idx[j][n]] = w;
					weight[segments_endpoints_idx[j][n]][segments_endpoints_idx[i][m]] = w;
				}
			}
		}
	}

	epsilon = 0;
	for (int i = 0; i < weight.size(); ++i){
		for (int j = 0; j < weight[i].size(); ++j){
			if (epsilon < weight[i][j])
				epsilon = weight[i][j];
		}
	}

	solvable = true;
	return solvable;
}
bool MyMTSP::Init(GraphArray *ngraph){
	vector<markerPath> paths;
	cout << "ngraph node 數量:" << ngraph->nodeNum;
	for (size_t i = 0; i <ngraph->array.size(); i++)
	{
		for each (pathNode* node in ngraph->array[i])
		{
			for (size_t n = 0; n < node->connects.size(); n++)
			{
				pathNode *next = node->connects[n];
				if (node->id < next->id){
					paths.push_back(markerPath(node, next, node->path[n]));
				}
			}
		}
	}
	cout << "總共" << paths.size() << "條path.";
	weight.clear();
	segments_endpoints_idx.clear();
	connectness.clear();
	solvable = false;
	for each (markerPath path in paths)
	{
		if (path.fullPath.size() < 2){
			printf("Not a line segment.\n");
			return solvable;
		}
	}
	num_node = 2 * paths.size() + 1;
	segments_endpoints_idx.resize(paths.size());
	for (int i = 0; i < paths.size(); ++i){
		segments_endpoints_idx[i].resize(2);
		segments_endpoints_idx[i][0] = 2 * i + 1;
		segments_endpoints_idx[i][1] = 2 * i + 2;
	}
	std::vector<bool> connectness_element(num_node, true);
	connectness.resize(num_node, connectness_element);
	for (int i = 0; i < paths.size(); ++i){
		int front_idx = paths[i].node1->id;
		int back_idx = paths[i].node2->id;
		int endpoint_idx_current[2] = { front_idx, back_idx };

		for (int j = i + 1; j < paths.size(); ++j){
			front_idx = paths[j].node1->id;
			back_idx = paths[j].node2->id;
			int endpoint_idx_compare[2] = { front_idx, back_idx };
			for (int m = 0; m < 2; ++m){
				for (int n = 0; n < 2; ++n){
					if (endpoint_idx_current[m] != endpoint_idx_compare[n]){
						connectness[segments_endpoints_idx[i][m]][segments_endpoints_idx[j][n]] = false;
						connectness[segments_endpoints_idx[j][n]][segments_endpoints_idx[i][m]] = false;
					}
				}
			}
		}
	}
	double u = 1;
	std::vector<double> weight_element(num_node, 0);
	weight.resize(num_node, weight_element);
	for (int i = 0; i < paths.size(); ++i){
		int n_point = paths[i].fullPath.size();
		OpenMesh::Vec3f point_front1 = paths[i].fullPath[0];
		OpenMesh::Vec3f point_front2 = paths[i].fullPath[1];
		OpenMesh::Vec3f point_back1 = paths[i].fullPath[n_point - 1];
		OpenMesh::Vec3f point_back2 = paths[i].fullPath[n_point - 2];
		OpenMesh::Vec3f endpoint_current[2] = { point_front1, point_back1 };
		OpenMesh::Vec3f tangentangle_current[2] = { (point_front1 - point_front2).normalize(), (point_back1 - point_back2).normalize() };
		for (int j = i + 1; j < paths.size(); ++j){
			n_point = paths[j].fullPath.size();
			point_front1 = paths[j].fullPath[0];
			point_front2 = paths[j].fullPath[1];
			point_back1 = paths[j].fullPath[n_point - 1];
			point_back2 = paths[j].fullPath[n_point - 2];
			OpenMesh::Vec3f endpoint_compare[2] = { point_front1, point_back1 };
			OpenMesh::Vec3f tangentangle_compare[2] = { (point_front1 - point_front2).normalize(), (point_back1 - point_back2).normalize() };
			for (int m = 0; m < 2; ++m){
				for (int n = 0; n < 2; ++n){
					double distance = (endpoint_current[m] - endpoint_compare[n]).length();
					double w = distance + u * (1.0 + OpenMesh::dot(tangentangle_current[m], tangentangle_compare[n]))*0.5;
					//double w = u * (1.0 + tangentangle_current[m] * tangentangle_compare[n])*0.5;
					weight[segments_endpoints_idx[i][m]][segments_endpoints_idx[j][n]] = w;
					weight[segments_endpoints_idx[j][n]][segments_endpoints_idx[i][m]] = w;
				}
			}
		}
	}

	epsilon = 0;
	for (int i = 0; i < weight.size(); ++i){
		for (int j = 0; j < weight[i].size(); ++j){
			if (epsilon < weight[i][j])
				epsilon = weight[i][j];
		}
	}

	solvable = true;
	return solvable;
}
bool MyMTSP::Solve(){
	if (!solvable){
		printf("Can not solve.\n");
		return false;
	}

	printf("Wire composition begin...\n");
	printf("Xi = %lf\n", epsilon);

	double clock_start, clock_end, cpu_time_used;
	clock_start = clock();

	IloEnv env;
	IloExpr expressions(env);
	IloRangeArray constraints(env);
	IloNumVarArray x(env);
	IloNumVarArray u(env);
	IloNumVar k;
	cout << "num_node:" << num_node;
	for (size_t i = 0; i < num_node*num_node; ++i){
		x.add(IloNumVar(env, 0, 1, ILOINT));
	}

	for (size_t i = 0; i < num_node; ++i){
		u.add(IloNumVar(env, 0, num_node, ILOINT));
	}

	k = IloNumVar(env, 1, num_node, ILOINT);

	for (int i = 0; i < num_node; ++i){
		for (int j = 0; j < num_node; ++j){
			if (i == j){
				constraints.add(x[i*num_node + j] == 0);
			}

			//object function
			expressions += weight[i][j] * x[i*num_node + j];

		}
	}
	expressions += epsilon * k;

	//constraint (5)
	for (int i = 1; i < num_node; ++i){
		IloExpr sum_in(env);
		IloExpr sum_out(env);
		for (int j = 0; j < num_node; ++j){
			if (i != j){
				sum_in += x[j*num_node + i];
				sum_out += x[i*num_node + j];
			}
		}
		constraints.add(sum_in == 1);
		constraints.add(sum_out == 1);
	}

	//constraint (6)
	IloExpr sum_in(env);
	IloExpr sum_out(env);
	for (int j = 1; j < num_node; ++j){
		sum_in += x[j*num_node + 0];
		sum_out += x[0 * num_node + j];
	}
	constraints.add((sum_in - k) == 0);
	constraints.add((sum_out - k) == 0);

	//constraint (7)
	constraints.add(u[0] == 0);
	for (int i = 1; i < num_node; ++i){
		constraints.add((u[i] + (num_node - 2)*x[i] - x[i*num_node]) - (num_node - 1) <= 0);
		constraints.add(u[i] + x[i] >= 2);
		for (int j = 1; j < num_node; ++j){
			if (i != j){
				constraints.add(u[i] - u[j] + num_node*x[i*num_node + j] + (num_node - 2)*x[j*num_node + i] - (num_node - 1) <= 0);
			}
		}
	}

	//constraints to ensure the endpoints of the same segment to be connected.
	for (int i = 0; i < segments_endpoints_idx.size(); ++i){
		int idx1 = segments_endpoints_idx[i][0];
		int idx2 = segments_endpoints_idx[i][1];
		constraints.add(x[idx1*num_node + idx2] + x[idx2*num_node + idx1] == 1);
	}

	//constraints to ensure the path to only go through the connected endpoints
	for (int i = 0; i < connectness.size(); ++i){
		for (int j = 0; j < connectness[i].size(); ++j){
			if (connectness[i][j] == false){
				constraints.add(x[i*num_node + j] + x[j*num_node + i] == 0);
			}
		}
	}
	//new constraint for my function 強制percent_y=0的點必須要是起點,強制percent_y =1的點必須是終點
	for each (int idx in dominant_startIdx)//強制dominant_startIdx內的node必定要作為虛擬節點0的下一點
	{
		constraints.add(x[idx] == 1);
	}
	for each(int idx in dominant_endIdx){//強制dominant_endIdx內的node必須要作為虛擬節點的上一點
		constraints.add(x[idx*num_node]==1);
	}
	//printf("Num of constraints: %d\n", constraints.getSize());
	IloModel model(env);
	model.add(IloMinimize(env, expressions));
	model.add(constraints);

	IloCplex cplex(model);
	cplex.setOut(env.getNullStream());
	//printf("Solving Begin...\n");
	//cout << "Solve start";
	if (!cplex.solve()){
		std::cout << "Fail to optimize the model." << std::endl;
		return false;
	}
	else{
		//std::cout << "Solved!" << std::endl;
	}

	results_x = IloNumArray(env);
	cplex.getValues(x, results_x);
	IloNumArray results_u(env);
	cplex.getValues(u, results_u);
	IloNum result_k = cplex.getValue(k);

	clock_end = clock();
	cpu_time_used = ((double)clock_end - clock_start) / CLOCKS_PER_SEC;
	printf("Done. Processing Time: %lf seconds.\n-----------------------\n", cpu_time_used);
	vector<pathNode*> startNode;
	return true;
}
vector<sequence> MyMTSP::getConnectSequence(){
	vector<sequence> results;
	vector<int> startIndex;
	for (int x = 1; x < num_node; x++)
	{
		if (results_x[x]==1)
		{
			startIndex.push_back(x);
		}
	}
	for each (int idx in startIndex)
	{
		int nextIdx = idx;
		sequence path;
		int count = 0;
		while (nextIdx!=0)
		{
			int pathsIdx = ((nextIdx - 1) / 2);
			bool positive = ((nextIdx-1)%2==0);
			if (count %2 == 0){
				path.push_back(indexAndDirect(pathsIdx, positive));
			}
			
			for (size_t x = 0; x < num_node; x++)
			{
				if (results_x[num_node*nextIdx + x])
				{
					nextIdx = x;
					break;
				}
			}
			count++;
		}
		results.push_back(path);
	}
	return results;
}