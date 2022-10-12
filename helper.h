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
map<int, int> nodePartition;		//�ڵ������map[�ڵ�index]��������
map<int, int> edgeIndex;			//��ߺͽ��ӳ���ϵ��map[����±�]����ڹ����е�index
vector<double> linkInitStatus;
EN_Project ph;
ofstream output("resource/test.txt");

int pressureNoDiffCount;

void initEpanetData();
void initCutEdge();
void initPartition();
void initPressure();

void init() {
	// ���ļ���ȡ����
	ifstream fin("resource/config.json");   // ע��˴������·��
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
		cout << "���ļ�ʧ��2��" << endl;
		exit(1);
	}

	initEpanetData();
	initCutEdge();
	initPartition();
	//initPressure();
}

void initEpanetData() {
	// ��inp�ļ���������
	//EN_Project ph;
	int errcode;
	EN_createproject(&ph);
	errcode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
	if (errcode != 0)
	{
		// Call functions that perform desired analysis
		cout << "open inp file err, code =" << errcode << endl;
	}

	errcode = EN_setdemandmodel(ph, EN_PDA, 10, 40, 0.5); //0���޸� ��С����ˮͷ ��modena 14 40 exnet��18 50�� 303��
	if (errcode != 0) {
		cout << "EN_setdemandmodel, code=" << errcode << endl;
	}

	// ��inp�ļ���ȡcount�������������̨
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
	//// �ص�����
	//EN_deleteproject(ph);
}

void initCutEdge() {
	string line;
	ifstream cutEdge(cutEdgeFile);
	if (!cutEdge) {
		cout << "���ļ�ʧ�ܣ�" << cutEdgeFile << endl;
		exit(1);//ʧ���˻ز���ϵͳ    
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
	cout << "����±꣺�ܶε�index" << endl;
	for (auto iter = edgeIndex.begin(); iter != edgeIndex.end(); iter++) {
		cout << iter->first << " : " << iter->second << endl;
	}
}

void closeEpanet() {
	// �ص�����
	EN_deleteproject(ph);
	output.close();
}

void initPartition() {
	ifstream partition(partitionFile);
	if (!partition) {
		cout << "���ļ�ʧ�ܣ�" << partitionFile << endl;
		exit(1);//ʧ���˻ز���ϵͳ    
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

//��ȡ���нڵ��ѹ��
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
	//// �ص�����
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
	double mean = sum / list.size(); //��ֵ
	double accum = 0.0;
	for (int i = 0; i < list.size(); i++) {
		accum += (list[i] - mean) * (list[i] - mean) / list.size();
	}
	double stdev = accum; //����
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

// ���룺һ����������1��ʾ�����ƣ�0��ʾ���ţ�������ʱ���и��
// hasPressureBelowThreshold ��ʾ�Ƿ���ѹ��С����ֵ
// �����Ŀ�꺯��ֵ������Ϊ3
// Ŀ��1�������Ƹ��� ��ԽСԽ�ã�
// Ŀ��2������©�������ͱ���,(A-B)/A��A��ʾ����ǰ������©������B��ʾ����solution�������©���� ��Խ��Խ�ã�
// Ŀ��3��ѹ�������ԣ��������ķ����ƽ���͵�ƽ��ֵ+���������е��ѹ��-�̶�ֵ����ƽ���͵�ƽ���� ��ԽСԽ�ã�
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
	errCode = EN_setdemandmodel(ph, EN_PDA, 10, 40, 0.5); //1���޸� ��С����ˮͷ ��modena 14 40 exnet��18 50�� step2��374��
	if (errCode != 0) {
		cout << "EN_setdemandmodel, code=" << errCode << endl;
	}
	for (int i = 0; i < solution.size(); i++) {
		errCode = EN_setlinkvalue(ph, edgeIndex[i], EN_INITSTATUS, EN_OPEN); //Ĭ�Ϲܶεĳ�ʼ״̬�Ǵ�״̬
		if (errCode != 0) {
			cout << "292 EN_setlinkvalue, code=" << errCode << "" << ", edgeIndex=" << edgeIndex[i] << ", edgeID=" << linkIndexToID[edgeIndex[i]] << endl;
		}
	}
	// user-supplied function to set parameters
	for (int i = 0; i < solution.size(); i++) {
		if (solution[i] == 0) {
			// ���ܶ�cutLinkIndex[i]����Ϊ����
			errCode = EN_setlinkvalue(ph, edgeIndex[i], EN_INITSTATUS, EN_CLOSED);
			if (errCode != 0) {
				cout << "301 EN_setlinkvalue, code=" << errCode << "" << ", edgeIndex=" << edgeIndex[i] << ", edgeID=" << linkIndexToID[edgeIndex[i]] << endl;
			}
		}
		else {
			// ���ܶ�cutLinkIndex[i]����Ϊ�����ƣ����账��
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
	//// �ص�����
	//EN_deleteproject(ph);

	//PrintNoDiffSolution(pressures, solution);
	// 
	int count1 = 0, count2 = 0;
	for (int i = 0; i < pressures.size(); i++) {
		if (skipPipeIDs.find(linkIndexToID[i+1]) != skipPipeIDs.end()) { 
			continue; //��config.json������ˮ��ˮ��Ƚڵ㣨ֻ����������ˮ�ڵ㣩
		}
		if (pressures[i] < 0) {//2���޸� �޸�369��380�У�modena 0 14 1 1 ��exnet 0 18 1 8�� ��һ����446��
			count1++;
		}
		if (pressures[i] < 10) { //modena��14��exnet��10
			count2++;
		}
	}
	if (count2 > 8) {
		*hasPressureBelowThreshold = true;
	}
	if (count1 > 5) { //����ѹ��<0�ĸ���Ϊ1��Ҳ����ѹ��С��0�ĸ������Ϊ1 ��Ӧ����ĳ���
		vector<double> rst;
		rst.push_back(DBL_MAX);//���ֵ
		rst.push_back(DBL_MAX);
		rst.push_back(DBL_MAX);//������ﲻ���㣬��log������ر���ֵ
		return rst;
	}
	if (count2 > 15) { //����ѹ��<14��С����ˮͷ�ĸ���Ϊ8��Ҳ����ѹ��С��14�ĸ������Ϊ8//14����С����ˮͷ������14����ˮ����Ϊ0 
		vector<double> rst;
		rst.push_back(DBL_MAX);
		rst.push_back(DBL_MAX);
		rst.push_back(DBL_MAX);
		return rst;
	}

	// Ŀ�꺯��1
	int obj1 = 0;
	for (int i = 0; i < solution.size(); i++) {
		if (solution[i] == 1) {
			obj1++;
		}
	}
	// Ŀ�꺯��2
	double obj2 = 0;
	double A = 0, B = 0;
	for (int i = 0; i < nodeCount; i++) {
		// �и�ѹʱ��������
		if (initPressures[i] > 0) {  //����Ŀ�꺯���Ѹ�ѹ���Ե�
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
	if (obj2 > 0) { //Ŀ�꺯��2�Ǹ�����������趨ΪС�ڸ�10����� 3���޸� modena ��0 Exnet��-2��
		vector<double> rst;
		rst.push_back(DBL_MAX);
		rst.push_back(DBL_MAX);
		rst.push_back(DBL_MAX);
		return rst;
	}
	// Ŀ�꺯��3
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
			totalVar2 += ((pressures[i] - 14) * (pressures[i] - 14)) / pressures.size(); //3���޸� �ĳɶ�Ӧ��������С����ˮͷ
		}
	}
	obj3 = (totalVar1 * 10 + totalVar2); //obj3 = (totalVar1*10 + totalVar2);

	// ����Ŀ�꺯������ԽСԽ��
	vector<double> objs;
	objs.push_back(obj1);
	objs.push_back(obj2);
	objs.push_back(obj3); //Ŀ�꺯���������ĸ��ͽ��Ǹ�ע�͵�
	return objs;
}