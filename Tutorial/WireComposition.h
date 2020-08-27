#pragma once
#define IL_STD
#include "ilcplex\ilocplex.h"
#include "ilconcert\iloexpression.h"
#include "meshCreater.h"
#include "imageConvert.h"
//用於使用最佳化來求解繞線順序
struct nodeMapping//已經被棄用
{
	vector<int> upNodeIndex;
	vector<int> lowNodeIndex;

	vector<vector<float>> weights;
	vector<IloNumVarArray*> optValues;
	pathNode *centerNode;
	static bool nodeIsVaild(pathNode *node){
		return node->connects.size() % 2 == 0;
	}
	void init(pathNode *node){
		centerNode = node;
		//cout << "初始化node" << node->id;
		int halfConnectNum = node->connects.size() / 2;
		vector<int> equalIndexs;
		for (size_t i = 0; i < node->connects.size(); i++)
		{
			pathNode *next = node->connects[i];
			if (next->y_2d < node->y_2d){
				upNodeIndex.push_back(i);
			}
			else if (next->y_2d > node->y_2d){
				lowNodeIndex.push_back(i);
			}
			else//相等的情況
			{
				equalIndexs.push_back(i);
			}
		}
		for each (int index in equalIndexs)//最後還是要處理y值相等的點,把索引值加入到比較小的那個組
		{
			if (upNodeIndex.size()<=lowNodeIndex.size())
			{
				upNodeIndex.push_back(index);
			}
			else
			{
				lowNodeIndex.push_back(index);
			}
		}
		weights.resize(upNodeIndex.size());
		//cout << "halfConnectNum:" << halfConnectNum << "upNodeIndex size:" << upNodeIndex.size() << "lowNodeIndex size:" << lowNodeIndex.size();
		for (int i = 0; i < upNodeIndex.size(); i++){
			int index1 = upNodeIndex[i];
			//cout << "index1:" << index1 << " x_2d:" << node->connects[index1]->x_2d<<"pos:("<<node->position<<")";
			OpenMesh::Vec3f vec1;
			if (node->path[index1].size()==0)
			{
				vec1 =(node->position - node->connects[index1]->position).normalize();
			}
			else{
				/*cout << "印出path[index1]:";
				for (size_t i = 0; i < node->path[index1].size(); i++)
				{
					cout << "(" << node->path[index1][i] << "),";
				}*/
				vec1 = (node->position - node->path[index1].front()).normalize();
			}
			cout << endl;
			for (int j = 0; j < lowNodeIndex.size(); j++){
				{
					int index2 = lowNodeIndex[j];
					OpenMesh::Vec3f vec2;
					if (node->path[index2].size()==0){
						vec2 = (node->connects[index2]->position- node->position).normalize();
					}
					else{
						vec2= (node->path[index2].front() - node->position).normalize();
					}
					//cout << "vec1:" << vec1 << " vec2:" << vec2 << "dot(vec1,vec2):" << dot(vec1, vec2);
					float dotResult = OpenMesh::dot(vec1, vec2);
					if (dotResult > 1){
						dotResult = 1;
					}
					else if (dotResult < -1){
						dotResult = -1;
					}
					float weight = acosf(dotResult)/M_PI;
					//cout << "index2:" << index2 << "x_2d:" << node->connects[index2]->x_2d << " weight:" << weight<<"angle:"<<acosf(weight)/M_PI<<"path[index2]:";
					/*for each (OpenMesh::Vec3f point in node->path[index2])
					{
						cout << "(" << point << "),  ";
					}
					cout <<"---"<<endl;*/
					weights[i].push_back(weight);
				}
			}
		}
		/*cout << "印出weight:";
		for (int i = 0; i < weights.size(); i++){
			cout << "{";
			for each (float w in weights[i])
			{
				cout << w << ",";
			}
			cout << "},";
		}
		cout << "node" << node->id<<"結束"<< endl;*/
	}
	void initOptValues(IloEnv env){
		cout << "initOptValues: node id:"<<centerNode->id;
		optValues.resize(weights.size());
		cout << "weight size:" << weights.size();
		for (size_t i = 0; i < weights.size(); i++)
		{
			optValues[i] = new IloNumVarArray(env);
			//vector<IloNumVarArray> numArray(1);
			
			cout << "i=" << i;
			cout << "optValues[i]->getSize" << optValues[i]->getSize();
			cout << " weights[i].size():" << weights[i].size()<<" ";
			for (size_t j = 0; j < weights[i].size(); j++)
			{
				cout << "j=" << j<<", ";
				optValues[i]->add(IloNumVar(env, 0, 1, ILOINT));
				cout << " add!";
			}
		}
	}
};
class WireComposition//已經被棄用
{
protected:
	GraphArray *graph;


public:
	//vector<vector<pathMapping>> mappings;
	//vector<vector<float>> NodeWeights;
	vector<nodeMapping> mappings;
	vector<int> solveFailIndexRecord;
	WireComposition(GraphArray *graph,vector<int> *startNodes,vector<int> *endNodes);
	~WireComposition();
	void test();
	bool Solve();
	void SolveOneByOne();
	vector<int> *startNodeids;
	vector<int> *endNodeids;
	vector<pair<pathNode*, pathNode*>> specialStartNodes;
};
struct markerPath{
	pathNode *node1;
	pathNode *node2;
	vector<OpenMesh::Vec3f> fullPath;
	markerPath(){

	}
	markerPath(pathNode *n1, pathNode *n2, vector<OpenMesh::Vec3f> path){
		node1 = n1;
		node2 = n2;
		fullPath.insert(fullPath.begin(), path.begin(), path.end());
		fullPath.insert(fullPath.begin(),node1->position);
		fullPath.push_back(node2->position);
	}
};
struct indexAndDirect{
	int index=-1;
	bool direct=true;//標記正反
	indexAndDirect(){
	}
	indexAndDirect(int idx,bool flag){
		index = idx;
		direct = flag;
	}
};
typedef vector<indexAndDirect> sequence;
class MyMTSP{
public:
	MyMTSP();
	bool Init(GraphArray *ngraph);//已被棄用
	bool Init(controlGraph pixelGraph);
	bool Solve();
	vector<sequence> getConnectSequence();
	//void SpecifySegmentSequence();
	//std::vector<std::vector<int>> GetSegmentSequence();
private:
	GraphArray *graph;
	bool solvable;
	int num_node;
	double epsilon;
	std::vector<std::vector<int>> segments_endpoints_idx;
	std::vector<std::vector<bool>> connectness;
	std::vector<std::vector<double>> weight;
	IloNumArray results_x;
	vector<int> dominant_startIdx;
	vector<int> dominant_endIdx;

};