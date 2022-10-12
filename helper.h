#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "metis.h"
#include "epanet2_2.h"
#include "json.hpp"
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <numeric>

using namespace std;

map<int, string> nodeIndexToID;
map<string, int> linkIDToIndex;
map<int, string> linkIndexToID;
map<string, int> nodeIDToIndex;
vector<double> nodeAllLinkLength;
vector<vector<int>> initSolution;
int nodeCount, tankCount, linkCount, patCount, curveCount, controlCount, ruleCount;
string inpFile, cutEdgeFile, partitionFile;
int nParts;
set<string> skipPipeIDs;
vector<double> initPressures;
map<int, int> nodePartition;		//节点分区，map[节点index]所属分区
map<int, int> edgeIndex;			//割边和解的映射关系，map[解的下标]割边在管网中的index
vector<double> linkInitStatus;
EN_Project ph;
ofstream output("resource/test.txt");

int pressureNoDiffCount;

void initEpanetData();
void initCutEdge();
void initPartition();
void initPressure();

void init() {
	// 从文件读取配置
	ifstream fin("resource/config.json");   // 注意此处是相对路径
	nlohmann::json j;
	fin >> j;
	fin.close();
	inpFile = j["inpFile"];
	cutEdgeFile = j["cutEdgeFile"];
	partitionFile = j["partitionFile"];
	nParts = j["nParts"];
	vector<string> skipPipeIDsVec = j["skipPipeIDs"].get<vector<string>>();
	initSolution = j["initSolution"].get<vector<vector<int>>>();
	skipPipeIDs = set<string>(skipPipeIDsVec.begin(), skipPipeIDsVec.end());
	skipPipeIDsVec.assign(skipPipeIDs.begin(), skipPipeIDs.end());

	if (!output) {
		cout << "打开文件失败2！" << endl;
		exit(1);
	}

	initEpanetData();
	initCutEdge();
	initPartition();
	//initPressure();
}

void initEpanetData() {
	// 从inp文件创建工程
	//EN_Project ph;
	int errcode;
	EN_createproject(&ph);
	errcode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
	if (errcode != 0)
	{
		// Call functions that perform desired analysis
		cout << "open inp file err, code =" << errcode << endl;
	}

	errcode = EN_setdemandmodel(ph, EN_PDA, 10, 40, 0.5); //0、修改 最小服务水头 （modena 14 40 exnet：18 50） 303行
	if (errcode != 0) {
		cout << "EN_setdemandmodel, code=" << errcode << endl;
	}

	// 从inp文件读取count，并输出到控制台
	EN_getcount(ph, EN_NODECOUNT, &nodeCount);
	EN_getcount(ph, EN_TANKCOUNT, &tankCount);
	EN_getcount(ph, EN_LINKCOUNT, &linkCount);
	EN_getcount(ph, EN_PATCOUNT, &patCount);
	EN_getcount(ph, EN_CURVECOUNT, &curveCount);
	EN_getcount(ph, EN_CONTROLCOUNT, &controlCount);
	EN_getcount(ph, EN_RULECOUNT, &ruleCount);
	cout << "node count = " << nodeCount << endl;
	cout << "tank count = " << tankCount << endl;
	cout << "link count = " << linkCount << endl;
	cout << "pat count = " << patCount << endl;
	cout << "curve count = " << curveCount << endl;
	cout << "control count = " << controlCount << endl;
	cout << "rule count = " << ruleCount << endl << endl;

	for (int i = 1; i <= nodeCount; i++) {
		char id[50];
		EN_getnodeid(ph, i, id);
		nodeIndexToID[i] = id;
		nodeIDToIndex[id] = i;
		nodeAllLinkLength.push_back(0);
	}
	vector<pair<int, int>> links(0);
	for (int i = 1; i <= linkCount; i++) {
		char id[50];
		EN_getlinkid(ph, i, id);
		linkIDToIndex[id] = i;
		linkIndexToID[i] = id;
		double status;
		EN_getlinkvalue(ph, i, EN_INITSTATUS, &status);
		linkInitStatus.push_back(status);
		int start, end;
		EN_getlinknodes(ph, i, &start, &end);
		links.push_back(make_pair(start, end));
		double length;
		EN_getlinkvalue(ph, i, EN_LENGTH, &length);
		nodeAllLinkLength[start - 1] += length;
		nodeAllLinkLength[end - 1] += length;
	}

	long t;
	EN_openH(ph);
	// initialize hydraulics; don't save them to file
	EN_initH(ph, EN_NOSAVE);
	// solve hydraulics
	EN_runH(ph, &t);
	// user-supplied function to process results
	for (int i = 1; i <= nodeCount; i++) {
		double value;
		EN_getnodevalue(ph, i, EN_PRESSURE, &value);
		initPressures.push_back(value);
	}
	EN_closeH(ph);
	//// 关掉工程
	//EN_deleteproject(ph);
}

void initCutEdge() {
	string line;
	ifstream cutEdge(cutEdgeFile);
	if (!cutEdge) {
		cout << "打开文件失败！" << cutEdgeFile << endl;
		exit(1);//失败退回操作系统    
	}
	getline(cutEdge, line);
	istringstream tmp(line);
	int edgeNum;
	tmp >> edgeNum;
	for (int i = 0; i < edgeNum; i++) {
		getline(cutEdge, line);
		istringstream tmp(line);
		int vecIndex;
		string linkID;
		tmp >> vecIndex >> linkID;
		edgeIndex[vecIndex] = linkIDToIndex[linkID];
	}
	cutEdge.close();
	cout << "解的下标：管段的index" << endl;
	for (auto iter = edgeIndex.begin(); iter != edgeIndex.end(); iter++) {
		cout << iter->first << " : " << iter->second << endl;
	}
}

void closeEpanet() {
	// 关掉工程
	EN_deleteproject(ph);
	output.close();
}

void initPartition() {
	ifstream partition(partitionFile);
	if (!partition) {
		cout << "打开文件失败！" << partitionFile << endl;
		exit(1);//失败退回操作系统    
	}
	string line;
	getline(partition, line);
	istringstream tmp(line);
	int num;
	tmp >> num;
	for (int i = 0; i < num; i++) {
		getline(partition, line);
		istringstream tmp(line);
		int nodeIndex, partitionNum;
		tmp >> nodeIndex >> partitionNum;
		nodePartition[nodeIndex] = partitionNum;
	}
	partition.close();
}

//获取所有节点的压力
void initPressure() {
	//EN_Project ph;
	//int errcode;
	//EN_createproject(&ph);
	//errcode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
	//if (errcode != 0)
	//{
	//	// Call functions that perform desired analysis
	//	cout << "open inp file err, code =" << errcode << endl;
	//}
	int  i;
	long t;
	EN_openH(ph);
	// initialize hydraulics; don't save them to file
	EN_initH(ph, EN_NOSAVE);
	// solve hydraulics
	EN_runH(ph, &t);
	// user-supplied function to process results
	for (int i = 1; i <= nodeCount; i++) {
		double value;
		EN_getnodevalue(ph, i, EN_PRESSURE, &value);
		initPressures.push_back(value);
	}
	EN_closeH(ph);
	//// 关掉工程
	//EN_deleteproject(ph);
	for (int i = 0; i < initPressures.size(); i++) {
		if (isnan(initPressures[i])) {
			cout << "initPressures[" << i << "] is nan" << endl;
		}
	}
}

double variance(vector<double> list) {
	double sum = 0;
	for (int i = 0; i < list.size(); i++) {
		sum += list[i];
	}
	double mean = sum / list.size(); //均值
	double accum = 0.0;
	for (int i = 0; i < list.size(); i++) {
		accum += (list[i] - mean) * (list[i] - mean) / list.size();
	}
	double stdev = accum; //方差
	//if (isnan(accum / (list.size() - 1))) {
	//	cout << "accum / (list.size() - 1) nan," << accum / (list.size() - 1) << ", size=" << list.size() << ", mean" << mean << ", sum" << sum << endl;
	//	for (int i = 0; i < list.size(); i++) {
	//		cout << list[i] << ", ";
	//	}
	//	cout << endl;
	//}
	//if (isnan(stdev)) {
	//	cout << "stdev nan," << stdev << endl;
	//}
	return stdev;
}

void PrintNoDiffSolution(vector<double> pressures, vector<int> solution) {
	if (pressures.size() != initPressures.size()) {
		cout << "pressure size not equal to initPressure size" << endl;
		return;
	}
	bool isEqual = true;
	for (int i = 0; i < pressures.size(); i++) {
		if (isnan(pressures[i]) || isnan(initPressures[i])) {
			continue;
		}
		if (pressures[i] != initPressures[i]) {
			isEqual = false;
			break;
		}
	}
	if (isEqual) {
		/*for (int i = 0; i < solution.size(); i++) {
			cout << solution[i] << ",";
		}
		cout << endl;*/
	}
}

// 输入：一个解向量（1表示流量计，0表示阀门），分区时的切割边
// hasPressureBelowThreshold 表示是否有压力小于阈值
// 输出：目标函数值，长度为3
// 目标1：流量计个数 （越小越好）
// 目标2：管网漏损量降低比率,(A-B)/A，A表示放置前管网总漏损量，B表示放置solution后管网总漏损量 （越大越好）
// 目标3：压力均衡性（各分区的方差的平方和的平均值+（管网所有点的压力-固定值）的平方和的平均） （越小越好）
vector<double> calculateObjectives(vector<int> solution, bool *hasPressureBelowThreshold) {
	//EN_Project ph;
	//int errcode;
	//EN_createproject(&ph);
	//errcode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
	//if (errcode != 0)
	//{
	//	// Call functions that perform desired analysis
	//	cout << "271 open inp file err, code =" << errcode << endl;
	//}
	int errCode;
	long t;
	errCode = EN_openH(ph);
	if (errCode != 0) {
		cout << "EN_openH, code=" << errCode << endl;
	}
	errCode = EN_setdemandmodel(ph, EN_PDA, 10, 40, 0.5); //1、修改 最小服务水头 （modena 14 40 exnet：18 50） step2在374行
	if (errCode != 0) {
		cout << "EN_setdemandmodel, code=" << errCode << endl;
	}
	for (int i = 0; i < solution.size(); i++) {
		errCode = EN_setlinkvalue(ph, edgeIndex[i], EN_INITSTATUS, EN_OPEN); //默认管段的初始状态是打开状态
		if (errCode != 0) {
			cout << "292 EN_setlinkvalue, code=" << errCode << "" << ", edgeIndex=" << edgeIndex[i] << ", edgeID=" << linkIndexToID[edgeIndex[i]] << endl;
		}
	}
	// user-supplied function to set parameters
	for (int i = 0; i < solution.size(); i++) {
		if (solution[i] == 0) {
			// 将管段cutLinkIndex[i]设置为阀门
			errCode = EN_setlinkvalue(ph, edgeIndex[i], EN_INITSTATUS, EN_CLOSED);
			if (errCode != 0) {
				cout << "301 EN_setlinkvalue, code=" << errCode << "" << ", edgeIndex=" << edgeIndex[i] << ", edgeID=" << linkIndexToID[edgeIndex[i]] << endl;
			}
		}
		else {
			// 将管段cutLinkIndex[i]设置为流量计，无需处理
		}
	}
	// initialize hydraulics; don't save them to file
	errCode = EN_initH(ph, EN_NOSAVE);
	if (errCode != 0) {
		cout << "EN_initH, code=" << errCode << endl;
	}
	
	// solve hydraulics
	errCode = EN_runH(ph, &t);
	if (errCode != 0) {
		vector<double> rst;
		rst.push_back(DBL_MAX);
		rst.push_back(DBL_MAX);
		rst.push_back(DBL_MAX);
		return rst;
	}
	//cout << "En_runH ok" << endl;
	// user-supplied function to process results
	vector<double> pressures;
	for (int i = 1; i <= nodeCount; i++) {
		double value;
		EN_getnodevalue(ph, i, EN_PRESSURE, &value);
		pressures.push_back(value);
	}
	EN_closeH(ph);
	vector<pair<int, int>> links(0);
	vector<double> newNodeAllLinkLength;
	for (int i = 1; i <= nodeCount; i++) {
		newNodeAllLinkLength.push_back(0);
	}
	for (int i = 1; i <= linkCount; i++) {
		int start, end;
		EN_getlinknodes(ph, i, &start, &end);
		links.push_back(make_pair(start, end));
		double length;
		EN_getlinkvalue(ph, i, EN_LENGTH, &length);
		newNodeAllLinkLength[start - 1] += length;
		newNodeAllLinkLength[end - 1] += length;
	}
	//// 关掉工程
	//EN_deleteproject(ph);

	//PrintNoDiffSolution(pressures, solution);
	// 
	int count1 = 0, count2 = 0;
	for (int i = 0; i < pressures.size(); i++) {
		if (skipPipeIDs.find(linkIndexToID[i+1]) != skipPipeIDs.end()) { 
			continue; //在config.json中跳过水箱水库等节点（只保留正常供水节点）
		}
		if (pressures[i] < 0) {//2、修改 修改369和380行（modena 0 14 1 1 ；exnet 0 18 1 8） 下一步在446行
			count1++;
		}
		if (pressures[i] < 10) { //modena是14，exnet是10
			count2++;
		}
	}
	if (count2 > 8) {
		*hasPressureBelowThreshold = true;
	}
	if (count1 > 5) { //容忍压力<0的个数为1，也就是压力小于0的个数最多为1 对应上面的程序
		vector<double> rst;
		rst.push_back(DBL_MAX);//最大值
		rst.push_back(DBL_MAX);
		rst.push_back(DBL_MAX);//如果这里不满足，在log中输出特别大的值
		return rst;
	}
	if (count2 > 15) { //容忍压力<14最小服务水头的个数为8，也就是压力小于14的个数最多为8//14是最小服务水头，低于14，需水量就为0 
		vector<double> rst;
		rst.push_back(DBL_MAX);
		rst.push_back(DBL_MAX);
		rst.push_back(DBL_MAX);
		return rst;
	}

	// 目标函数1
	int obj1 = 0;
	for (int i = 0; i < solution.size(); i++) {
		if (solution[i] == 1) {
			obj1++;
		}
	}
	// 目标函数2
	double obj2 = 0;
	double A = 0, B = 0;
	for (int i = 0; i < nodeCount; i++) {
		// 有负压时会变成无穷
		if (initPressures[i] > 0) {  //计算目标函数把负压忽略掉
			A += 0.5 * pow(10, -8) * pow(initPressures[i], 1.18) * nodeAllLinkLength[i];
		}
	}
	for (int i = 0; i < nodeCount; i++) {
		if (pressures[i] > 0) {
			B += 0.5 * pow(10, -8) * pow(pressures[i], 1.18) * newNodeAllLinkLength[i];
		}
	}
	obj2 = (B - A) / A * 100;
	if (isnan(obj2)) {
		cout << "obj2 is nan" << endl;
	}
	if (obj2 > 0) { //目标函数2是负数，将这个设定为小于负10的舍掉 3、修改 modena 是0 Exnet是-2）
		vector<double> rst;
		rst.push_back(DBL_MAX);
		rst.push_back(DBL_MAX);
		rst.push_back(DBL_MAX);
		return rst;
	}
	// 目标函数3
	double obj3 = 0;
	vector<vector<double>> partPressures(nParts);
	for (int i = 1; i <= nodeCount; i++) {
		if (pressures[i - 1] > 0) {
			partPressures[nodePartition[i]].push_back(pressures[i - 1]);
		}
	}
	double totalVar1 = 0;
	for (int i = 0; i < nParts; i++) {
		totalVar1 += variance(partPressures[i]) / nParts;
	}
	double totalVar2 = 0;
	for (int i = 0; i < pressures.size(); i++) {
		if (pressures[i] > 0) {
			totalVar2 += ((pressures[i] - 14) * (pressures[i] - 14)) / pressures.size(); //3、修改 改成对应管网的最小服务水头
		}
	}
	obj3 = (totalVar1 * 10 + totalVar2); //obj3 = (totalVar1*10 + totalVar2);

	// 所有目标函数都是越小越好
	vector<double> objs;
	objs.push_back(obj1);
	objs.push_back(obj2);
	objs.push_back(obj3); //目标函数不运行哪个就将那个注释掉
	return objs;
}