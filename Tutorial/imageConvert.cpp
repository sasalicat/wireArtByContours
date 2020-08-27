#include "stdafx.h"
#include "imageConvert.h"
#include "WireComposition.h"
imageConvert::imageConvert(string filename)
{
	image = cv::imread(filename, cv::IMREAD_GRAYSCALE);
	if (!image.data)
	{
		cout << "讀取圖片" << filename << "失敗!"<<endl;
		system("pause");
	}
	//cout << "Mat type:" << image.type()<<"size:"<<image.size();
	if (useSkeleton){
		calSkeleton();
		cv::imwrite("../Resources/skel.png", skeletonImg);
	}
	else
	{
		cv::threshold(image, skeletonImg, skeletonThreshold, 255, cv::THRESH_BINARY);
		cv::imshow("after threshold", skeletonImg);
		cv::waitKey(0);
		cv::destroyAllWindows();
	}
	createData();

}
bool undirectContain(vector<pair<pixelPoint, pixelPoint>> edgeSet, pair<pixelPoint, pixelPoint> edge){
	for each (pair <pixelPoint,pixelPoint> e in edgeSet)
	{
		//cout << "e:(" << e.first.x << "," << e.first.y << ") >> (" << e.second.x << "," << e.second.y << "),    ";
		if (e.first.x == edge.first.x&&e.first.y == edge.first.y&&e.second.x == edge.second.x&&e.second.y == edge.second.y){
			//cout << endl;
			return true;
		}
		if (e.first.x == edge.second.x&&e.first.y == edge.second.y&&e.second.x == edge.first.x&&e.second.y == edge.first.y){
			//cout << endl;
			return true;
		}
	}
	//cout << endl;
	return false;
}
int indexof(vector<pixelPoint*> list, pixelPoint* elem){
	//cout << "list:";
	for (size_t i = 0; i < list.size(); i++)
	{
		//cout << "(" << list[i]->x << "," << list[i]->y << ")";

		if ((list[i]->x == elem->x) && (list[i]->y == elem->y)){
			//cout << "return" << i;
			return i;
		}
	}
	return -1;
}
int indexof(vector<int> list, int elem){
	for (size_t i = 0; i < list.size(); i++)
	{
		//cout << "(" << list[i]->x << "," << list[i]->y << ")";

		if (list[i] == elem){
			//cout << "return" << i;
			return i;
		}
	}
	return -1;
}
vector<pixelPoint*> imageConvert::findLastConnectPoint(int posx, int posy){
	vector<pixelPoint*> result;
	int x = posx - 1;
	col lastcol = datas[x];

	//cout<< "in findLastConnectPoint x:" << posx << "y:" << posy <<endl;
	if (x < 0){
		return result;
	}
	for (int y = posy - 1; y <= posy+1; y++){//對於點上一個col的上中下三點
		
		if (y >= 0 && y < skeletonImg.size().height){
			
			if (skeletonImg.at<unsigned char>(y, x)>0){//找到一個點
				int closePtIdx = -1;
				for (size_t i = 0; i < lastcol.height; i++)//對與那一列(col)所有點進行對比
				{
					//cout  <<"lastcol address:" << &lastcol<<"height:"<<lastcol.height<<endl;
					
					pixelPoint* pt = lastcol.points[i];
					//cout << "i="<<i<<" pt address" << pt<<", ";
					if (y-pt->y <=mergeYThreshold&&y-pt->y>=0){
						if (closePtIdx<0||abs(pt->y - y) < abs(lastcol.points[closePtIdx]->y-y)){//找y相差最小的pixelPoint
							closePtIdx = i;
						}
					}
				}
				if (indexof(result, lastcol.points[closePtIdx]) < 0){
					
					result.push_back(lastcol.points[closePtIdx]);
				}
			}

		}
	}
	//cout << endl;
	return result;
}
void imageConvert::createData(){
	
	width = skeletonImg.size().width;
	int height = skeletonImg.size().height;
	mergeYThreshold = 0.5*height;
	datas = new col[width];
	totalPixelPointNum=0;
	for ( int x = 0; x < width; x++)
	{
		//cout << "x=" << x << ":";
		vector<pixelPoint*> list;
		pixelPoint *firstMeetPoint = NULL;//此變數case2查看上一次case1遇到的是哪個pixelPoint,之後連接
		for (int y = 0; y < height; y++){
			//case1:從上往下第一次遇到有值像素
			if ((y == 0 && skeletonImg.at<unsigned char>(y, x)>0) || (skeletonImg.at<unsigned char>(y, x)>0 && skeletonImg.at<unsigned char>(y - 1, x) == 0)){
				pixelPoint *newPt = new pixelPoint(x, y, (float)x / (float)(width - 1), (float)y / (float)(height - 1), totalPixelPointNum++);
				//cout << "case1 newPt x:" << x << " y:" << y;
				flatArray.push_back(newPt);
				int suby = y;
				vector<pixelPoint*>last;
				while (suby<height&&skeletonImg.at<unsigned char>(suby, x)>0){//向下找所有同樣大於0的點,直到某個點能連接到其他pixelPoint或下方沒有任何>0的點了
					//cout << "suby:" << suby;
					last = findLastConnectPoint(x, suby);
					if (last.size() > 0){//找到能連接到的點
						//newPt->lastPoints = last;
						for each(pixelPoint *pt in last){
							bool found = false;
							for each (pixelPoint *ptHave in newPt->lastPoints)
							{
								if (pt == ptHave)
								{
									found = true;
									break;
								}
							}
							if (!found){
								newPt->lastPoints.push_back(pt);
							}
						}
					}
					suby++;
					if (suby - y > mergeYThreshold){//如果像素值依然是1但是超出了認為是同一個點的距離
						y = suby;	
						break;
					}
				}
				//cout << "停在y=" << y << "suby=" << suby<<" ";
				y = suby;
				//cout << "搜索后的lastPoint";
				for each(pixelPoint* pt in newPt->lastPoints){
					//cout << "(" << pt->x << "," << pt->y << ")";
					pt->nextPoints.push_back(newPt);
				}
				//cout<< endl;
				firstMeetPoint = newPt;
				list.push_back(newPt);
			}
			/*
			else if ((y == height - 1 && skeletonImg.at<unsigned char>(y, x)>0) || (skeletonImg.at<unsigned char>(y, x)>0 && skeletonImg.at<unsigned char>(y+1,x)==0))//理論上來說,這種情況只會出現在上一總情況之後,這裡算是撿漏,如果像素值>並且其下面的一個像素為0的話說明到達了>0像素行的尾端,則添加一個像素
			{
				pixelPoint *newPt = new pixelPoint(x, y, (float)x / (float)(width - 1), (float)y / (float)(height - 1), totalPixelPointNum++);
				//cout << "case2 newPt x:" << x << " y:" << y;
				flatArray.push_back(newPt);
				vector<pixelPoint*> last = findLastConnectPoint(x, y);
				if (last.size() > 0){//找到能連接到的點
					//newPt->lastPoints = last;
					for each(pixelPoint *pt in last){
						bool found = false;
						for each (pixelPoint *ptHave in newPt->lastPoints)
						{
							if (pt == ptHave)
							{
								found = true;
								break;
							}
						}
						if (!found){
							newPt->lastPoints.push_back(pt);
						}
					}
				}
				for each(pixelPoint* pt in newPt->lastPoints){
					
					pt->nextPoints.push_back(newPt);
				}
				if (firstMeetPoint != NULL){//理論上來說這裡一定會成立,因為case2一定是在case1之後發生
					firstMeetPoint->lastPoints.push_back(newPt);
					//newPt->nextPoint.push_back(firstMeetPoint);
				}
				list.push_back(newPt);
			}
			*/

		}
		
		col nowCol(list.data(), list.size());
		datas[x] = nowCol;
	}

	//修正y值.因為演算法是以從上到下第一個觸碰到的非0點為交點,所以當使用的elem線比一個像素更粗的時候,底端的交點就會出現y!=1的情況(因為碰到最下一列的像素不是作為交點坐標的像素)
	bool pixelMeetEnd = false;
	for (int x = 0; x < width; x++)//看看是不是有像素觸底
	{
		if (skeletonImg.at<unsigned char>(height - 1, width)>0){
			pixelMeetEnd = true;
			break;
		}
	}
	if (pixelMeetEnd){//開始修正,將y最大的像素的percent_y拉到1,其他像素根據自己的percent_y分攤這個變形量
		//cout << "修正datas:";
		float maxPercent_y=0;
		for (int x = 0; x < width;x++)
		{
			pixelPoint* colEnd = datas[x].points[datas[x].height - 1];
			if (maxPercent_y < colEnd->percentage_y){
				maxPercent_y = colEnd->percentage_y;
			}
		}
		for (int x = 0; x < width; x++){
			for (int y = 0; y < datas[x].height; y++)
			{
				datas[x].points[y]->percentage_y /= maxPercent_y;
			}
		}
	}

	//找交點,內部交點是lastPoints或nextPoints多於一個的點,潛在外部交點是percent_y=0或1且percent_x不等於0或1的點
	for (size_t x = 0; x < width; x++)
	{
		col colx = datas[x];
		for (size_t y = 0; y < colx.height; y++)
		{
			pixelPoint* pt = colx.points[y];
			//找內部交點
			//if (pt->lastPoints.size()>1 || pt->nextPoint.size()>1)
			if (pt->nextPoints.size()+pt->lastPoints.size()>2)
			{
				pt->intersection = true;
				intersectionPoint.push_back(pt);

			}
			//找潛在外部交點,為什麼潛在外部交點的定義是percent_y=0或1且percent_x不等於0或1的點
			//原因是潛在外部交點是為了在sample的時候額外切出來,這樣後面在創建pathGraph的時候不同distance的兩個ring之間的交點才會體現出來
		}
	}
}
void imageConvert::calSkeleton(){
	cv::Mat im;
	cv::threshold(image, im, skeletonThreshold, 255, cv::THRESH_BINARY);
	//cout << "skeletonImg type:" << skeletonImg.type() << "image type:" << image.type();
	cv::Mat morp_elem = getStructuringElement(cv::MORPH_CROSS, cv::Size(3, 3));

	cv::Mat skel = cv::Mat::zeros(im.size(), CV_8UC1);
	cv::Mat temp = cv::Mat::zeros(im.size(), CV_8UC1);
	//cout << "0.skel type:" << skel.type()<<"temp type:"<<temp.type() << endl;
	int i = 0;
	while (true){
		cout << "i=" << i<<" ";
		cv::morphologyEx(im, temp, cv::MORPH_OPEN, morp_elem);
		//cout << "1.skel type:" << skel.type() << "temp type:" << temp.type() << endl;
		cv::bitwise_not(temp, temp);
		//cout << "2.skel type:" << skel.type() << "temp type:" << temp.type() << endl;
		cv::bitwise_and(im, temp, temp);
		//cout << "3.skel type:" << skel.type() << "temp type:" << temp.type() << endl;
		//cv::imshow("第i次提取的輪廓", temp);
		//cv::waitKey(0);
		//cout << "skel size:" << skel.size()<<" type:"<<skel.type()<< " temp size:" << temp.size()<<"type:"<<temp.type()<<endl;
		cv::bitwise_or(skel, temp, skel);
		
		cv::erode(im, im, morp_elem);

		double max_val = -1;
		cv::minMaxLoc(im, NULL, &max_val, NULL, NULL);
		if (max_val == 0){
			break;
		}
		i++;
	}
	skeletonImg = skel;
	cv::imwrite("skel.jpg", skeletonImg);
	//cv::imshow("提取的骨骼:", skel);
	//cv::waitKey(0);
	//cv::destroyAllWindows();
}


vector<vector<pixelPoint>> imageConvert::sample(int Number){
	float interval = (float)width / (float)Number;
	vector<int> intervalHitX;
	for (int i = 1; i < Number; i++){//舉個例子如果width是5,interval是1,則Number是5,則需要記錄的HitX只有{1,2,3,4},處於x=0和x=5的點是sample的起點終點,一定會被記錄的
		intervalHitX.push_back(interval*i);
	}

	//vector<pixelPoint*> haveExtend;
	vector<pair<pixelPoint, pixelPoint>> extendRecord;//用於記錄已經採樣過了的連線,放置反向重複採樣
	vector<pair<pixelPoint*,pixelPoint*>> header;//採樣的起點,first是第一個pixelPoint,second則是header所來自的交點
	//vector<bool> headerMarker;//用於標記對應索引的header是向next還是last,向next時為true,反之為false
	vector<vector<pixelPoint>> result;
	for (int i = 0; i < datas[0].height; i++){
		//if (datas[0].points[i]->intersection){
			for each (pixelPoint* next in datas[0].points[i]->nextPoints)
			{
				//cout << "add next:(" << next->x << "," << next->y << ")";
				header.push_back(pair<pixelPoint*, pixelPoint*>(next, datas[0].points[i]));
			}
	}
	while (header.size() > 0){
		
		pixelPoint* now = header[0].first;
		pixelPoint* lineStartPoint = now;

		vector<pixelPoint> line;
		line.push_back(*header[0].second);
		//cout << "now->nextPoints.size:" << now->nextPoint.size()<<" intersection:"<<now->intersection;
		bool directToNext = header[0].first->x > header[0].second->x;//說明前進方向為正,即是朝node的next方向前進
		if (directToNext)
		{
			while (now->nextPoints.size() > 0 && !now->intersection){//如果當前是端點則nextpoint為0個,如果當前是交點則nextpoint多於1個

				if (indexof(intervalHitX, now->x)>=0){
					line.push_back(*now);
				}
				else if (now->percentage_y == 0||now->percentage_y==1)
				{
					line.push_back(*now);
				}
				now = now->nextPoints[0];//移動到下一個
			}
		}
		else
		{
			while (now->lastPoints.size() > 0 && !now->intersection){//如果當前是端點則nextpoint為0個,如果當前是交點則nextpoint多於1個
				if (indexof(intervalHitX, now->x)>=0){
					line.push_back(*now);
				}
				else if (now->percentage_y == 0 || now->percentage_y == 1)
				{
					line.push_back(*now);
				}
				now = now->lastPoints[0];//移動到下一個
			}
		}
		if (now != lineStartPoint){//有可能第一個點就直接是交點,這種情況不會進while,所以now會==lineStartPoint,這種情況不用在額外把自己加進去,因為一開始已經加過了
			line.push_back(*now);//因為交點不會進上面的while迴圈所以這裡要額外把交點加進去
		}
		
		header.erase(header.begin());
		//headerMarker.erase(headerMarker.begin());
		//cout << "終止于點(" << now->x << "," << now->y << ") intersection" << now->intersection << "nextPoint size:" << now->nextPoint.size();
		pair<pixelPoint, pixelPoint> edge(line.front(),line.back());
		if (!undirectContain(extendRecord, edge))
		{

			for each(pixelPoint* p in now->nextPoints)
			{
				header.push_back(pair<pixelPoint*, pixelPoint*>(p, now));//first是擴張路徑的起點,second是此起點的由來
				//headerMarker.push_back(true);
			}
			for each(pixelPoint* p in now->lastPoints)
			{
				header.push_back(pair<pixelPoint*, pixelPoint*>(p, now));//first是擴張路徑的起點,second是此起點的由來
				//headerMarker.push_back(false);
			}
			result.push_back(line);
			extendRecord.push_back(edge);
		}
		else{//如果這個終點已經被擴張過了,line被廢棄,因為肯定有一條line的反向存在

		}
		//cout << endl;
	}
	return result;
}
void imageConvert::debugDatas(){
	for (int y = 0; y < width; y++){
		cout << y << ":";
		for (int x = 0; x < datas[y].height; x++){
			pixelPoint* pt = datas[y].points[x];
			cout << "(" << pt->x << "," << pt->y << ");";
		}
		cout << "     ";
	}
	cout << "debug datas 結束";
}
imageConvert::~imageConvert()
{
}

bool leftKeyDown = false;
int activePointIndex =-1;
controlGraph *mainGraph = NULL;
const int TRAGET_INDEX=2;
void mouseEventInConnect(int event, int x, int y, int flags, void *param)
{
	switch (event)
	{
		case CV_EVENT_LBUTTONDOWN:{
			leftKeyDown = true;
			break;
		}
		case CV_EVENT_LBUTTONUP:{
			leftKeyDown = false;
			break;
		}

	}
	if (leftKeyDown){
		int width = (float)mainGraph->Image->size().width / 3;
		float percent_x = (float)(x-width) / width;
		percentPoint tragetPos(percent_x, (float)y / (float)mainGraph->Image->size().height);
		//mainGraph->contorlPoints[TRAGET_INDEX]->line->deform(mainGraph->contorlPoints[TRAGET_INDEX]->point, tragetPos);
		imageConvert::drawBaseGraph(mainGraph->Image,mainGraph->mainGraphs);
		imageConvert::drawGraphPoints(mainGraph->Image,mainGraph->mainGraphs);
		imageConvert::drawControlPoints(mainGraph->Image, mainGraph->contorlPoints);
		//cv::circle(*mainGraph->Image, cv::Point2f(x, y), 4, cv::Scalar(0, 255, 0), 2);
		cv::imshow("graphShow", *mainGraph->Image);
	}
}
void imageConvert::drawBaseGraph(cv::Mat *image, vector<graphLine*> graph){
	(*image) = cv::Scalar(0, 0, 0);//清屏

	int width = image->size().width / 3; int height = image->size().height;
	for (size_t t = 0; t < 3; t++)
	{
		for each (graphLine *line in graph)
		{
			int begin_x = t*width;
			for (size_t i = 0; i < line->points->size() - 1; i++)
			{
				cv::Point2i p1(begin_x + width*line->points->operator[](i).percentage_x, height*line->points->operator[](i).percentage_y);
				cv::Point2i p2(begin_x + width*line->points->operator[](i + 1).percentage_x, height*line->points->operator[](i+1).percentage_y);
				cv::line(*image, p1, p2, cv::Scalar(225, 255, 255), 1);
			}
		}
	}
}
void imageConvert::drawBaseGraph(cv::Mat *image, vector<graphLine*> graph, int width, int height){
	(*image) = cv::Scalar(0, 0, 0);//清屏
	for (size_t t = 0; t < 3; t++)
	{
		for each (graphLine *line in graph)
		{
			int begin_x = t*width;
			for (size_t i = 0; i < line->points->size() - 1; i++)
			{
				cv::Point2i p1(begin_x + width*line->points->operator[](i).percentage_x, height*line->points->operator[](i).percentage_y);
				cv::Point2i p2(begin_x + width*line->points->operator[](i + 1).percentage_x, height*line->points->operator[](i+1).percentage_y);
				cv::line(*image, p1, p2, cv::Scalar(225, 255, 255), 1);
			}
		}
	}
}
void  imageConvert::drawControlPoints(cv::Mat *image, vector<controlPoint*> cpts){
	int width = image->size().width / 3; int height = image->size().height;
	for each (controlPoint* cpt in cpts)
	{
		if (cpt->matchPoint!=NULL)//已經有配對的點
			cv::circle(*image, cv::Point2i(width + cpt->points[0]->percentage_x*width, cpt->points[0]->percentage_y*height),3,cv::Scalar(0,255,0),2);
		else
			cv::circle(*image, cv::Point2i(width + cpt->points[0]->percentage_x*width, cpt->points[0]->percentage_y*height), 3, cv::Scalar(0, 0, 255), 2);
	}
}
void imageConvert::drawGraphPoints(cv::Mat *image, vector<graphLine*>graph){
	int width = image->size().width / 3; int height = image->size().height;
	for each (graphLine *line in graph)
	{
		int begin_x = width;

		for (size_t i = 0; i < line->points->size(); i++)
		{
			cv::Point2i pos(begin_x + width*line->points->operator[](i).percentage_x, height*line->points->operator[](i).percentage_y);
			cv::circle(*image,pos ,2,cv::Scalar(255,0,0),2);

		}
	}
}
void imageConvert::drawNextGraph(cv::Mat *image, vector<graphLine*> graph, float start_percentX, int width, int height){
	for each (graphLine *line in graph)
	{
		int begin_x = start_percentX*width;
		for (size_t i = 0; i < line->points->size() - 1; i++)
		{
			cv::Point2i p1(begin_x + width*line->points->operator[](i).percentage_x, height*line->points->operator[](i).percentage_y + height);
			cv::Point2i p2(begin_x + width*line->points->operator[](i + 1).percentage_x, height*line->points->operator[](i+1).percentage_y + height);
			cv::line(*image, p1, p2, cv::Scalar(225, 255, 255), 1);
		}
	}
}
controlPointGroup imageConvert::horizontalMatch(vector<controlPoint*> controls, float margeThreshold){
	controlPointGroup result;
	//cout << "controls:";
		
	for each (controlPoint* cpt in controls)
	{
		//cout << "(" << cpt->points[0]->percentage_x << "," << cpt->points[0]->percentage_y << ")";
		if (cpt->points[0]->percentage_x == 0.0f){
			result.group1.push_back(cpt);
		}
		else if (cpt->points[0]->percentage_x == 1.0f)
		{
			result.group2.push_back(cpt);
		}
	}

	for (int i = 0; i < result.group1.size(); i++){
		//cout << "group1 第" << i << "點:(" << result.group1[i]->points[0]->percentage_x << "," << result.group1[i]->points[0]->percentage_y << ")";
		for each (controlPoint* testPoint in result.group2)//
		{
			if (abs(testPoint->points[0]->percentage_y - result.group1[i]->points[0]->percentage_y) < margeThreshold){
				if (abs(testPoint->points[0]->percentage_y - 0.5f) > abs(result.group1[i]->points[0]->percentage_y - 0.5f)){//testPoint比較接近圖片的角落
					percentPoint pos(result.group1[i]->points[0]->percentage_x, testPoint->points[0]->percentage_y);
					result.group1[i]->deform(pos);
					//cout << "變形到:(" << pos.x << "," << pos.y << ");";
				}
				else if (abs(testPoint->points[0]->percentage_y - 0.5f) < abs(result.group1[i]->points[0]->percentage_y - 0.5f)){//當前點比較接近圖片角落
					percentPoint pos(testPoint->points[0]->percentage_x, result.group1[i]->points[0]->percentage_y);
					testPoint->deform(pos); //line->deform(testPoint->point, pos);
				}
				result.group1[i]->matchPoint = testPoint;
				testPoint->matchPoint = result.group1[i];
			}
		}
	}
	return result;
}
bool imageConvert::horizonalIsOk(controlPointGroup match){
	//cout << "horizonalIsOk開始" << endl;
	int connectCount = 0;
	for each (controlPoint* point in match.group1)
	{
		if (point->matchPoint!=NULL)
		{
			connectCount++;
		}
	}
	if (connectCount>0)//橫向連接成立的條件是至少有一個連接點是有連接的
	{
		return true;
	}
	else
	{
		return false;
	}
}
controlPointGroup imageConvert::verticalMatch(vector<controlPoint*> controls, float margeThreshold,bool& success){
	controlPointGroup result;
	//cout << "controls:";
	for each (controlPoint* cpt in controls)
	{
		//cout << "(" << cpt->points[0]->percentage_x << "," << cpt->points[0]->percentage_y << ")";
		if (cpt->points[0]->percentage_y == 0.0f){
			int insertIdx=0;
			for (; insertIdx <result.group1.size(); insertIdx++)
			{
				if (result.group1[insertIdx]->points[0]->percentage_x >= cpt->points[0]->percentage_x){
					break;
				}

			}
			//cout << "group1 insert at:" << insertIdx;
			result.group1.insert(result.group1.begin()+insertIdx,cpt);
			//result.group1.push_back(cpt);
		}
		else if (cpt->points[0]->percentage_y == 1.0f)
		{
			int insertIdx = 0;
			for (; insertIdx <result.group2.size(); insertIdx++)
			{
				if (result.group2[insertIdx]->points[0]->percentage_x >= cpt->points[0]->percentage_x){
					break;
				}
			}
			//cout << "group2 insert at:" << insertIdx;
			result.group2.insert(result.group2.begin() + insertIdx, cpt);
			//result.group2.push_back(cpt);
		}
	}
	percentPoint upLeft(0, 0);
	percentPoint upRight(1, 0);
	percentPoint btmLeft(0, 1);
	percentPoint btmRight(1, 1);
	if (result.group1.size() == result.group2.size() || abs((int)result.group1.size() - (int)result.group2.size()) == 1){
		success = false;
		if (result.group1.size() == result.group2.size()){
			success = true;
		}
		else{
			if (result.group1.size() > (int)result.group2.size()){//如果圖片上方的點數量比下方多1個說明理論上有一個重合點位於y=0 和y=1,只有這樣每一層交點數量才會一致
				if (result.group1.front()->points[0]->percentage_x <= margeThreshold){//由於誤差,有可能交點並不是剛剛好位於0,如果小於修正範圍,則把該點y值修正到0
					result.group1.front()->deform(upLeft);
				}
				if (1 - result.group1.back()->points[0]->percentage_x <= margeThreshold){
					result.group1.back()->deform(upRight);
				}
				if (result.group1.front()->points[0]->percentage_x == 0.0f&&result.group1.back()->points[0]->percentage_x == 1.0f){
				//if (result.group1.front()->points[0]->percentage_x == result.group1.back()->points[0]->percentage_y){//判斷是不是完全修正成功了,因為有兩個點如果有其中一個沒法修正的話也不能算是完全成功了
					success = true;
				}
			}
			else//如果圖片下方的點數量比上方多1
			{
				if (result.group2.front()->points[0]->percentage_x <= margeThreshold){
					result.group2.front()->deform(btmLeft);
				}
				if (1 - result.group2.back()->points[0]->percentage_x <= margeThreshold){
					result.group2.back()->deform(btmRight);
				}
				if (result.group2.front()->points[0]->percentage_x == 0.0f&&result.group2.back()->points[0]->percentage_x == 1.0f){
					success = true;
				}
			}
		}
		if (success){
			bool margeEndPoint = result.group2.front()->points[0]->percentage_x == 0.0f&&result.group2.back()->points[0]->percentage_x == 1.0f;//result.group2.front()->points[0]->percentage_x == result.group2.back()->points[0]->percentage_y;
			vector<controlPoint*> splicingControlPoints;//拼接兩個
			splicingControlPoints.insert(splicingControlPoints.end(), result.group2.begin(), result.group2.end());
			if (margeEndPoint)//跳過一個y=0的點來實現合併相同點的效果
				splicingControlPoints.insert(splicingControlPoints.end(), result.group2.begin()+1, result.group2.end());
			else
				splicingControlPoints.insert(splicingControlPoints.end(), result.group2.begin(), result.group2.end());
			int pointCanChoose = splicingControlPoints.size() - result.group1.size() + 1;
			if (result.group2.size()==1)//如果下面的點只有一個的話對應哪個都好
			{
				result.group2[0]->matchPoint = result.group1[0];
				for each (controlPoint* cpt in result.group1)
				{
					cpt->matchPoint = result.group2[0];
				}
			}

			for (size_t i = 0; i < pointCanChoose; i++)
			{
				float error = 0;
				for (size_t j = 0; j < result.group1.size()-1; j++)
				{
					
					float splicing2 = splicingControlPoints[i + j + 1]->points[0]->percentage_x;
					if (i + j + 1 >= result.group2.size()){
						splicing2 += 1;
					}
					float splicing1 = splicingControlPoints[i + j]->points[0]->percentage_x;
					if (i + j >= result.group2.size()){
						splicing1 += 1;
					}

					error += abs((splicing2 - splicing1) - (result.group1[j + 1]->points[0]->percentage_x - result.group1[j]->points[0]->percentage_x));

				}
				if (error == 0){
					for (size_t j = 0; j < result.group1.size(); j++){
						controlPoint* matchPoint = NULL;
						if (margeEndPoint){
							matchPoint = result.group2[i + j%(result.group2.size()-1)];//這會導致圖片下方的控制點的最後一點一定沒有matchpoint,不過在這個情況下最後一點就是最前一點
						}
						else
						{
							matchPoint = result.group2[i + j%result.group2.size()];
						}
						result.group1[j]->matchPoint = matchPoint;
						matchPoint->matchPoint = result.group1[j];
					}
					break;
				}
				else if (error<mergeVerticalThreshold)//通過改動group1的percent_x來強制match
				{
					for (size_t j = 0; j < result.group1.size() - 1; j++){
						float splicing2 = splicingControlPoints[i + j + 1]->points[0]->percentage_x;
						if (i + j + 1 >= result.group2.size()){
							splicing2 += 1;
						}
						float splicing1 = splicingControlPoints[i + j]->points[0]->percentage_x;
						if (i + j >= result.group2.size()){
							splicing1 += 1;
						}
						float correctInterval = splicing2 - splicing1;
						result.group1[j + 1]->deform(percentPoint(result.group1[j]->points[0]->percentage_x + correctInterval,0));
					}
					for (size_t j = 0; j < result.group1.size(); j++){
						controlPoint* matchPoint = NULL;
						
						if (margeEndPoint){
							matchPoint = result.group2[(i + j) % (result.group2.size() - 1)];//這會導致圖片下方的控制點的最後一點一定沒有matchpoint,不過在這個情況下最後一點就是最前一點
						}
						else
						{
							//cout << "i+j%result.group2.size()=" << i + j%result.group2.size();
							matchPoint = result.group2[(i + j)%result.group2.size()];
						}
						//cout << "點j=" << j << "matchpoint address:" << matchPoint << endl;
						result.group1[j]->matchPoint = matchPoint;
						matchPoint->matchPoint = result.group1[j];
					}
					break;

				}
			}

		}
		else
		{
			return result;
		}
	}
	return result;
}
bool imageConvert::inExtraStartDirect(int id1, int id2){
	for each (pair<int,int> direct in extraStartDirect)
	{
		if (direct.first==id1&&direct.second ==id2)
		{
			return true;
		}
	}
	return false;
}
void imageConvert::showGraphicConnect(vector<vector<pixelPoint>> &graph){
	int ori_w = skeletonImg.size().width;
	int ori_h = skeletonImg.size().height;
	cout << "ori_w:" << ori_w << ",ori_h:" << ori_h;
	float ratio = (float)ori_w / (float)ori_h;
	int height = 300;//200;
	int width = ratio*height;
	cv::Mat  newone(height, width*3, CV_8UC3, cv::Scalar(0,0,0));
	int count = 0;

	//切割輸入的graph,用lineGraphs的keyPoints如果不切割的話後面的nodeGraph無法找到可以連接的node
	int origrn_graph_size = graph.size();
	for (size_t i = 0; i <origrn_graph_size; i++)
	{
		vector<pixelPoint*> interKpts;
		for (size_t j = 0; j < graph[i].size(); j++)
		{
			pixelPoint *kpt = &graph[i][j];
			if (kpt->percentage_y==0)
			{
				//cout << "kpt:(" << kpt->percentage_x << "," << kpt->percentage_y << "):last size:" << kpt->lastPoints.size() << " next size:" << kpt->nextPoint.size() << ";";
			}
			//位於圖片底端或者頂端但又不是位於最左或最右的中間交點需要額外切割出來
			if ((kpt->percentage_y == 0 || kpt->percentage_y == 1) && (kpt->percentage_x != 0 && kpt->percentage_x != 1) && (kpt->lastPoints.size()>0 && kpt->nextPoints.size()>0))
			{
				interKpts.push_back(kpt);
			}
		}
		if (interKpts.size() > 0){//需要切割
			vector<vector<pixelPoint>> segmPart;
			int beginSegmIdx = 0;
			for each (pixelPoint* ikpt in interKpts)
			{
				int endSegmIdx = -1;
				for (int j = 0; j < graph[i].size(); j++){
					if (graph[i][j].percentage_x == ikpt->percentage_x&&graph[i][j].percentage_y == ikpt->percentage_y){
						endSegmIdx = j;
						break;
					}
				}
				if (endSegmIdx < 0){
					cout << "錯誤的切割!endSegmIdx:" << endSegmIdx;
					system("pause");
				}
				vector<pixelPoint> newPart;
				for (size_t j = beginSegmIdx; j <= endSegmIdx; j++)
				{
					newPart.push_back(graph[i][j]);
				}
				segmPart.push_back(newPart);
				beginSegmIdx = endSegmIdx;
			}
			vector<pixelPoint> newPart;
			for (size_t j = beginSegmIdx; j < graph[i].size(); j++)
			{
				newPart.push_back(graph[i][j]);
			}
			segmPart.push_back(newPart);
			//切割完成但是要放回去并沒有那麼簡單,隨便放會改變graph的結構,舉個例子 graph有3個line 編號0,1,2,如果切割了0,生成3,4,如果直接移除0,把3,4push back
			//graph的結構會變為:1,2,3,4,看起來沒問題,但是如果發現lineGraphs[1] interKpts>0需要切割,此時會回去找graph[1],但是此時索引值1指向的是原來的2!
			graph[i] = segmPart[0];
			for (int j = 1; j < segmPart.size(); j++)
			{
				graph.push_back(segmPart[j]);
			}
			//同樣是上面case 0,1,2切割后產生 3,4,以這個順序插入后結果會是3,1,2,4,不會影響下一次切割
		}
	}

	vector<graphLine*> lineGraphs(graph.size());
	//cout << "lineGraphs總共有" << graph.size() << "條line";
	for (size_t i = 0; i < graph.size(); i++)
	{
		lineGraphs[i] = new graphLine(&graph[i],ori_w,ori_h);
		//cout << "lineGraph address:" << lineGraphs[i] << endl;
		//cout <<"i="<<i <<"points:";
		for each (pixelPoint point in graph[i])
		{
			//cout << "(" << point.x << "," << point.y << ")";
		}
	}
	//cout << endl;
	for (size_t i = 0; i < lineGraphs.size(); i++)
	{
		for (size_t j = 0; j < lineGraphs.size(); j++)
		{
			if (i != j){
				lineGraphs[i]->tryConnect(lineGraphs[j]);
			}
		}
		lineGraphs[i]->finshConnect();
	}
	vector<pixelPoint*> keyPoints;
	vector<controlPoint*> controlPoints;

	for each (graphLine *lineGraph in lineGraphs)
	{
		//keyPoints.insert(keyPoints.end(), lineGraphs->keyPoints.begin(), lineGraphs->keyPoints.end());
		for each (pixelPoint* kpt in lineGraph->keyPoints)
		{
			bool haveControlPoint = false;
			for each (controlPoint* cpt in controlPoints)
			{
				float cpt_x = cpt->points[0]->percentage_x;
				float cpt_y = cpt->points[0]->percentage_y;
				if (cpt_x == kpt->percentage_x&&cpt_y == kpt->percentage_y){
					cpt->points.push_back(kpt);
					cpt->lines.push_back(lineGraph);
					haveControlPoint = true;
				}
			}
			if (!haveControlPoint){//如果沒有找到和自身坐標相同的控制點
				controlPoint *cpt = new controlPoint();
				cpt->points.push_back(kpt);
				cpt->lines.push_back(lineGraph);
				controlPoints.push_back(cpt);
			}
		}

	}
	/*cout << "印出control Points:";
	for each (controlPoint* cpt in controlPoints)
	{
		cout << "(" << cpt->points[0]->percentage_x << "," << cpt->points[0]->percentage_y << ") point num:" << cpt->points.size() <<"match:"<<cpt->matchPoint<< endl;
	}*/
	
	controlPointGroup matchResult_horizontal = horizontalMatch(controlPoints, 0.02f);

	bool horizontal_ok = horizonalIsOk(matchResult_horizontal);
	mainGraph = new controlGraph();
	mainGraph->mainGraphs = lineGraphs;
	mainGraph->Image = &newone;
	mainGraph->contorlPoints = controlPoints;
	imageConvert::drawBaseGraph(mainGraph->Image, lineGraphs);
	imageConvert::drawGraphPoints(mainGraph->Image, lineGraphs);
	imageConvert::drawControlPoints(mainGraph->Image, controlPoints);

	
	//cout << "horizontal_ok:" << horizontal_ok;
	if (!horizontal_ok){
		cv::putText(newone,"error! none horizontal connect",cv::Point2i(10,height/2),cv::FONT_HERSHEY_SIMPLEX,0.7,cv::Scalar(75,75,255));
	}
	cv::imshow("graphShow", newone);
	//cv::setMouseCallback("graphShow",mouseEventInConnect);
	cv::waitKey(0);
	bool vertical_ok = false;
	controlPointGroup matchResult_vertical = verticalMatch(controlPoints, 0.02f, vertical_ok);
	cv::Mat image2(height*2, width * 3, CV_8UC3, cv::Scalar(0, 0, 0));
	imageConvert::drawBaseGraph(&image2, lineGraphs,width,height);
	float start_x = matchResult_vertical.group1[0]->matchPoint->points[0]->percentage_x - matchResult_vertical.group1[0]->points[0]->percentage_x + 1;//如果在這裡報Access violation reading location 0x0000000000000000.一般是縱向連接找不到匹配的點 matchResult_vertical.group1[0]->matchPoint ==NULL導致
	imageConvert::drawNextGraph(&image2, lineGraphs,start_x,width,height);
	float firstPointOffset = matchResult_vertical.group1[0]->points[0]->percentage_x;
	//cout << "畫點在(" << start_x << "," << height << ")";
	if (!vertical_ok){
		cv::putText(image2, "error! none  vertical connect", cv::Point2i(10, height / 2), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(75, 75, 255));
		cv::circle(image2, cv::Point2i((start_x+firstPointOffset)*width, height), 3, cv::Scalar(0, 0, 255), 2);
	}
	else{
		cv::circle(image2, cv::Point2i((start_x+firstPointOffset)*width, height), 3, cv::Scalar(0, 255, 0), 2);
	}
	cv::imshow("graphShow", image2);
	//cv::setMouseCallback("graphShow",mouseEventInConnect);
	cv::waitKey(0);
	if (horizontal_ok&&vertical_ok){
		this->lineGraph = lineGraphs;
		this->Shift_percent = matchResult_vertical.group1[0]->matchPoint->points[0]->percentage_x - matchResult_vertical.group1[0]->points[0]->percentage_x;//解釋下Shift_percent,這個值的意義是每個level起始點要向內shift多少grid長度的百分比,
		if (this->Shift_percent < 0){
			this->Shift_percent += 1.0f;
		}
		//cout << "shift_percent:" << this->Shift_percent;
		//這個值是有上一牌控制點的第一個(result2.group1[0])所匹配的點的percentage_x決定的
	}
	cv::destroyAllWindows();
	MyMTSP mtsp;
	mtsp.Init(*mainGraph);
	if (mtsp.Solve()){
		//cout << "mtsp solve success! great!";
		vector<sequence> connectSequence = mtsp.getConnectSequence();
		cv::Mat  demoImg(height, width, CV_8UC3, cv::Scalar(0, 0, 0));

		for each (sequence seq1 in connectSequence)
		{
			int startId;
			int nextId;
			if (seq1[0].direct)
			{
				startId = graph[seq1[0].index].front().id;
				nextId = graph[seq1[0].index].back().id;
			}
			else
			{
				startId = graph[seq1[0].index].back().id;
				nextId = graph[seq1[0].index].front().id;
			}
			int endId;
			if (seq1.back().direct){
				endId= graph[seq1.back().index].back().id;
			}
			else{
				endId= graph[seq1.back().index].front().id;
			}


			connectDirect[startId][-1] = nextId;
			if (flatArray[startId]->percentage_y != 0||flatArray[endId]->percentage_y!=1){
				extraStartDirect.push_back(pair<int, int>(startId,nextId));
			}
			for (size_t i = 0; i <seq1.size()-1; i++)
			{
				int nowId;
				int lastId;
				int nextId;
				if (seq1[i].direct){
					nowId = graph[seq1[i].index].back().id;
					lastId = graph[seq1[i].index].front().id;
				}
				else
				{
					nowId = graph[seq1[i].index].front().id;
					lastId = graph[seq1[i].index].back().id;
				}
				if (seq1[i+1].direct){
					if (nowId != graph[seq1[i+1].index].front().id)
					{
						cout << "連續的序列交點不是同一個點!";
						system("pause");
					}
					nextId = graph[seq1[i+1].index].back().id;
				}
				else
				{
					if (nowId != graph[seq1[i+1].index].back().id)
					{
						cout << "連續的序列交點不是同一個點!";
						system("pause");
					}
					nextId = graph[seq1[i+1].index].front().id;
				}
				connectDirect[nowId][lastId] = nextId;
			}
		}
		cv::destroyAllWindows();

	}
	else
	{
		cout << "solve fail......";
	}
	
}
