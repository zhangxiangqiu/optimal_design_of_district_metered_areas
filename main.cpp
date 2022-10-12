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
#include "mopso.h"

#define SIGNAL_NUM -1000

using namespace std;

struct Link
{
	int index;
	string id;
	int node1;
	int node2;
	string nodeID1;
	string nodeID2;
	double length;
	double diameter;
};
struct CutEdge 
{
	int index;
	string id;
	double diameter; 
	int start;
	int end;
	string nodeID1; 
	string nodeID2; 
	double length;
	int startPart;
	int endPart; 
};

bool judge(Link a, Link b) {
	return a.node1 < b.node1;
}

vector<idx_t> part;
vector<Link> allLinks;
map<int, vector<Link>> nodeLinksMap;
vector<double> nodeBaseDemand(0);
vector<double> nodeTotalLength;

void genGraphFromInpFile(string filename) {

	EN_Project ph;
	int errcode;
	EN_createproject(&ph);
	errcode = EN_open(ph, (char*)filename.c_str(), "resource/report.rpt", "");
	if (errcode != 0)
	{
		
		cout << "open inp file err, code =" << errcode << endl;
	}

	
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

	int nodeWeightNum = 2, edgeWeightNum = 0; 
	set<int> nodes;
	vector<pair<int, int>> links(0);
	vector<int> edgeWeights(0);
	vector<int> nodeWeights(nodeCount * nodeWeightNum);
	nodeTotalLength = vector<double>(nodeCount);
	for (int i = 1; i <= nodeCount; i++) {
		char id[50];
		EN_getnodeid(ph, i, id); 
		nodeIndexToID[i] = id;
		double pattern;
		if (i < 600) {
			pattern = 1.20399;
		}
		else if (i < 1200) {  
			pattern = 1.25399;
		}
		else {
			pattern = 1.15399;
		}
		double baseDemand;
		EN_getnodevalue(ph, i, EN_BASEDEMAND, &baseDemand);
		nodeBaseDemand.push_back(baseDemand * pattern);
	}
	for (int i = 1; i <= linkCount; i++) { 
		int start, end;
		EN_getlinknodes(ph, i, &start, &end);
		links.push_back(make_pair(start, end));
		double length, diameter;
		EN_getlinkvalue(ph, i, EN_LENGTH, &length);
		EN_getlinkvalue(ph, i, EN_DIAMETER, &diameter);
		edgeWeights.push_back(diameter);
		nodeTotalLength[start] += length;
		nodes.insert(start);
		nodes.insert(end);
		char linkID[50], nodeID1[50], nodeID2[50];
		EN_getlinkid(ph, i, linkID);
		EN_getnodeid(ph, start, nodeID1);
		EN_getnodeid(ph, end, nodeID2); 
		Link link = {
			i,
			linkID,
			start,
			end,
			nodeID1,
			nodeID2,
			length,
			diameter,
		};
		allLinks.push_back(link);
		
	}
	
	EN_deleteproject(ph);

	for (int i = 0; i < allLinks.size(); i++) {
		for (int j = i + 1; j < allLinks.size(); j++) {
			if ((allLinks[i].node1 == allLinks[j].node1) && (allLinks[i].node2 == allLinks[j].node2)) {
				cout << "注意这里，有重复管段。管段ID 《" << allLinks[i].id << "》和管段ID《" << allLinks[j].id << "》的起始点重复，都是<" << allLinks[i].nodeID1 << "," << allLinks[i].nodeID2 << ">" << endl;
			}
		}
	}

	double totalBaseDemand = 0.0;
	double totalNodeLength = 0.0;
	for (int i = 0; i < nodeCount; i++) {
		totalBaseDemand += nodeBaseDemand[i];
		totalNodeLength += nodeTotalLength[i];
	}

	int index = 0; 

	for (int i = 0; i < nodeCount; i++) {
		
		nodeWeights[index++] = nodeBaseDemand[i]; 
		nodeWeights[index++] = nodeTotalLength[i]; 
		
	}

	for (int i = 0; i < allLinks.size(); i++) {
		nodeLinksMap[allLinks[i].node1].push_back(allLinks[i]);
		nodeLinksMap[allLinks[i].node2].push_back(allLinks[i]);
	}

	// 写入graph.txt（带权重）
	ofstream outgraph("resource/graph.txt");
	if (!outgraph) {
		cout << "打开文件失败3！" << endl;
		exit(1);
	}
	outgraph << nodes.size() << " " << links.size() << " " << nodeWeightNum << " " << edgeWeightNum << endl;
	for (int i = 0; i < nodes.size(); i++) {
		for (int j = 0; j <= links.size(); j++) {
			int first = links[j].first - 1;
			int second = links[j].second - 1;
			if (first == i) {
				outgraph << second << " " << edgeWeights[j] << " ";
			}
			else if (second == i) {
				outgraph << first << " " << edgeWeights[j] << " ";
			}
		}
		outgraph << SIGNAL_NUM;
		for (int j = 0; j < nodeWeightNum; j++) {
			outgraph << " " << nodeWeights[i * nodeWeightNum + j];
		}
		if (i != nodes.size() - 1) {
			outgraph << endl;
		}
	}
	outgraph.close();
}

vector<idx_t> splitGraph(vector<idx_t> xadj, vector<idx_t> adjncy, vector<idx_t> adjwgt, vector<idx_t> vwgt, int n_Parts, bool useKWay, bool useEdgeWeight) {
	idx_t nVertices = xadj.size() - 1; 
	idx_t nEdges = adjncy.size() / 2;    
	idx_t nWeights = vwgt.size() / nVertices;	
	idx_t nParts = n_Parts;   
	idx_t objval;
	std::vector<idx_t> part(nVertices, 0);

	idx_t options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_UFACTOR] = 500; 
	

	cout << "nWeights=" << nWeights << endl;

	idx_t* adjwgtData = adjwgt.data();
	if (!useEdgeWeight) {
		adjwgtData = NULL;
	}

	int ret = 0;
	if (useKWay) {
		ret = METIS_PartGraphKway(&nVertices, &nWeights, xadj.data(), adjncy.data(), vwgt.data(), NULL, adjwgtData, &nParts, NULL, NULL, options, &objval, part.data()); //倒数第四个是看是否是多权重，如果是多权重应该设置多个数值
		
	} 
	else {
		ret = METIS_PartGraphRecursive(&nVertices, &nWeights, xadj.data(), adjncy.data(), vwgt.data(), NULL, adjwgtData, &nParts, NULL, NULL, options, &objval, part.data());
		
	}
	if (ret == METIS_OK) {
		cout << "graph split success, code=" << ret << endl << endl;
	}
	else {
		cout << "graph split fail, code=" << ret << endl << endl;
	}
	cout << "graph split objval " << objval << endl;

	
	
	return part;
}

void genNewInpFromResult(string filename) {
	// 从inp文件创建工程
	EN_Project ph;
	int errcode;
	EN_createproject(&ph);
	errcode = EN_open(ph, (char*)filename.c_str(), "resource/report.rpt", "");
	if (errcode != 0)
	{
		
		cout << "open inp file err, code =" << errcode << endl;
	}
	for (int i = 1; i <= nodeCount; i++) {
		EN_setnodevalue(ph, i, EN_ELEVATION, (part[i-1] % 5) * 25 + 10);
	}
	string s = "resource/inp/Modena111.inp";
	EN_saveinpfile(ph, (char*)s.c_str());

	EN_deleteproject(ph);
}
void generateObjectiveFile(string filename) { 
	ofstream outobjective(filename);
	if (!outobjective) {
		cout << "打开文件失败3！" << endl;
		exit(1);
	}
	vector<double> baseDemandVal(nParts);
	vector<double> lengthVal(nParts);
	double t = 0;
	for (int i = 0; i < part.size(); i++) {
		baseDemandVal[part[i]] += nodeBaseDemand[i];
		lengthVal[part[i]] += nodeTotalLength[i];
		t += nodeBaseDemand[i];
	}
	for (int i = 0; i < nParts; i++) {
		outobjective << i;
		outobjective << "," << baseDemandVal[i];
		outobjective << "," << lengthVal[i] << endl;
	}
	outobjective.close();
}
void generateCutEdgeFile(string filename, string filename1) {
	vector<CutEdge> cutEdges;
	for (int i = 0; i < allLinks.size(); i++) {
		Link link = allLinks[i];
		int start = link.node1;
		int end = link.node2;
		int part1 = part[start - 1];
		int part2 = part[end - 1];
		if (part1 != part2) {
			CutEdge cutEdge = {
			 link.index,
			 link.id,
			 link.diameter,
			 start,
			 end,
			 link.nodeID1,
			 link.nodeID2,
			 link.length,
			 part1,
			 part2,
			};
			cutEdges.push_back(cutEdge);
		}
	}
	ofstream out(filename);
	if (!out) {
		cout << "打开文件失败4！" << endl;
		exit(1);
	}
	ofstream out1(filename1);
	if (!out1) {
		cout << "打开文件失败5！" << endl;
		exit(1);
	}
	out << cutEdges.size() << endl;
	out << setw(5) << "index" << setw(10) << "id" << setw(10) << "diameter" << setw(10) << "start" << setw(10) << "end" << setw(10) << "startID" << setw(10) << "endID" << setw(10) << "length" << setw(10) << "startPart" << setw(10) << "endPart" << endl;
	for (int i = 0; i < cutEdges.size(); i++) {
		CutEdge ce = cutEdges[i];
		out << setw(5) << ce.index << setw(10) << ce.id << setw(10) << ce.diameter << setw(10) << ce.start << setw(10) << ce.end << setw(10) << ce.nodeID1 << setw(10) << ce.nodeID2 << setw(10) << ce.length << setw(10) << ce.startPart << setw(10) << ce.endPart << endl;
	}
	out1 << cutEdges.size() << endl;
	for (int i = 0; i < cutEdges.size(); i++) {
		CutEdge ce = cutEdges[i];
		out1 << i << "," << ce.id << endl;
	}
	out.close();
	out1.close();
}

void run_Split() {
	
	std::ifstream fin("resource/config.json");   
	nlohmann::json j;
	fin >> j;
	fin.close();
	string filename = j["inpFile"];
	bool useKWay = j["useKWay"];
	string cutEdgeFile = j["cutEdgeFile"];
	string partitionFile = j["partitionFile"];
	nParts = j["nParts"];

	
	genGraphFromInpFile(filename);

	
	ifstream ingraph("resource/graph.txt");
	if (!ingraph) {
		cout << "打开文件失败！" << endl;
		exit(1);
	}

	
	int vexnum, edgenum, vexWeightNum, edgeWeightFlag;
	string line;
	getline(ingraph, line);
	istringstream tmp(line);
	tmp >> vexnum >> edgenum >> vexWeightNum >> edgeWeightFlag;
	vector<idx_t> xadj(0);
	vector<idx_t> adjncy(0); 
	vector<idx_t> adjwgt(0); 
	vector<idx_t> vwgt(0);	

	idx_t aa, ww;
	for (int i = 0; i < vexnum; i++) {
		xadj.push_back(adjncy.size());
		getline(ingraph, line);
		istringstream tmp(line);
		while (tmp >> aa >> ww) {
			if (aa == SIGNAL_NUM) {
				vwgt.push_back(ww);
				for (int j = 1; j < vexWeightNum; j++) {
					tmp >> ww;
					vwgt.push_back(ww);
				}
			}
			else {
				adjncy.push_back(aa);
				adjwgt.push_back(ww);
			}
		}
	}
	xadj.push_back(adjncy.size());
	ingraph.close();

	bool useEdgeWeight = true;
	if (edgeWeightFlag == 0) {
		useEdgeWeight = false;
	}


	part = splitGraph(xadj, adjncy, adjwgt, vwgt, nParts, useKWay, useEdgeWeight);


	ofstream outpartition(partitionFile);
	if (!outpartition) {
		cout << "打开文件失败2！" << endl;
		exit(1);
	}
	ofstream outobjective("resource/objective.csv");
	if (!outpartition) {
		cout << "打开文件失败3！" << endl;
		exit(1);
	}
	vector<int> objVal(nParts * vexWeightNum);
	outpartition << part.size() << endl;
	for (int i = 0; i < part.size(); i++) {
		outpartition << i + 1 << " " << part[i] << endl;
		for (int j = 0; j < vexWeightNum; j++) {//节点权重个数
			objVal[part[i] * vexWeightNum] += vwgt[i * vexWeightNum];
			objVal[part[i] * vexWeightNum + 1] += vwgt[i * vexWeightNum + 1];
		}
		//cout << nodeIndexToID[i + 1] << " " << part[i] << endl;
	}
	for (int i = 0; i < objVal.size(); i += vexWeightNum) {
		outobjective << i;
		for (int j = 0; j < vexWeightNum; j++) {
			outobjective << "," << objVal[i + j];
		}
		outobjective << endl;
	}
	outpartition.close();
	outobjective.close();
	generateObjectiveFile("resource/objective.csv");

	genNewInpFromResult(filename);
	generateCutEdgeFile("resource/allCutEdge.txt","resource/allCutEdge1.csv");
}

void run_Mopso() { // 多目标优化
	srand((unsigned)time(NULL));

	Multi_Object_Particle_Swarm_Optimization mospo;
	mospo.MOPSO();
}

void testEpanet() {
	string inpFile = "resource/inp/exnet2modified.inp"; // 可以修改inp
	vector<double> initPressures;
	vector<double> newPressures;
	
	
	for (int i = 0; i < 1; i++) {
		EN_Project ph;
		int errcode;
		EN_createproject(&ph);
		errcode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
		if (errcode != 0)
		{
			// Call functions that perform desired analysis
			cout << "open inp file err, code =" << errcode << endl;
		}
		EN_getcount(ph, EN_NODECOUNT, &nodeCount);
		EN_getcount(ph, EN_LINKCOUNT, &linkCount);
		cout << "count1 = " << linkCount;
		
		for (int i = nodeCount; i >= nodeCount - 5; i--) {
			char value[50];
			EN_getnodeid(ph, i, value);
			cout << value << endl;
		}
		long t = 0;
		EN_openH(ph);
		// initialize hydraulics; don't save them to file
		EN_initH(ph, EN_NOSAVE);
		// solve hydraulics
		int errCode = EN_runH(ph, &t);
		// user-supplied function to process results
		for (int i = nodeCount; i >= nodeCount - 5; i--) {
			double value;
			EN_getnodevalue(ph, i, EN_PRESSURE, &value);
			cout << setiosflags(ios::fixed) << setprecision(6) << value << endl;
		}
		
		EN_closeH(ph);
	
		EN_deleteproject(ph);
	}


	vector<string> solution1 = { "179" };
	EN_createproject(&ph);
	int errcode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
	if (errcode != 0)
	{
		// Call functions that perform desired analysis
		cout << "open inp file err, code =" << errcode << endl;
	}
	EN_getcount(ph, EN_LINKCOUNT, &linkCount);
	vector<int> solution;
	for (int i = 0; i < solution1.size(); i++) {
		int index;
		EN_getnodeindex(ph, (char*)solution1[i].c_str(), &index);
		solution.push_back(index);
	}
	long t = 0;
	EN_openH(ph);

	// user-supplied function to set parameters
	for (int i = 0; i < solution.size(); i++) {
		
		EN_setlinkvalue(ph, solution[i], EN_INITSTATUS, EN_CLOSED);
	}
	for (int i = 1; i <= linkCount; i++) {
		double value;
		EN_getlinkvalue(ph, i, EN_INITSTATUS, &value);
		cout << value << ",";
	}
	cout << endl;
	for (int i = 0; i < solution.size(); i++) {
		cout << solution[i] << ",";
	}
	cout << endl;
	// initialize hydraulics; don't save them to file
	EN_initH(ph, EN_NOSAVE);

	// solve hydraulics
	errcode = EN_runH(ph, &t);
	if (errcode != 0) {
		cout << "EN_runH err, code=" << errcode << endl;
	}
	// user-supplied function to process results
	for (int i = 1; i <= nodeCount; i++) {
		double value;
		EN_getnodevalue(ph, i, EN_PRESSURE, &value);
		newPressures.push_back(value);
	}
	double flow;
	EN_getlinkvalue(ph, 151, EN_FLOW, &flow);
	cout << flow << endl;
	EN_closeH(ph);

	EN_deleteproject(ph);

	
}

void run_Pressure() {
	string inpFile = "";
	string partFile = "";
	string rstFile = "";
	vector<string> solution = { "2173","2407","2443","2538","2586","2717","2899","2945","2979","3030","3037","3181","3479","3549","3730","4894","5053","5241","5165","2757","2619","3047","2365","2656","3166","2904","4041","3172","2139","2679" };//输入阀门的id(0对应的）还要修改skippipe modena "269", "270", "271", "272" Exnet "3001", "3002", "9003", "9004", "9005", "9006", "9007"
	EN_Project ph;//上面的是关闭的阀门的id
	int errCode;
	long t;
	vector<int> nodeIndex;
	vector<int> nodePartition;
	vector<double> initPressures;
	vector<double> newPressures;
	int nodeCount, linkCount;

	

	EN_createproject(&ph);
	errCode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
	if (errCode != 0)
	{
	
		cout << "open inp file err, code =" << errCode << endl;
		return;
	}
	EN_getcount(ph, EN_NODECOUNT, &nodeCount);
	EN_getcount(ph, EN_LINKCOUNT, &linkCount);

	errCode = EN_openH(ph);
	if (errCode != 0) {
		cout << "EN_openH, code=" << errCode << endl;
		return;
	}
	errCode = EN_setdemandmodel(ph, EN_PDA, 10, 40, 0.5);// 4、修改 modena 14 40 exnet 18 50 
	if (errCode != 0) {
		cout << "EN_setdemandmodel, code=" << errCode << endl;
		return;
	}

	errCode = EN_initH(ph, EN_NOSAVE);
	if (errCode != 0) {
		cout << "EN_initH, code=" << errCode << endl;
		return;
	}

	
	errCode = EN_runH(ph, &t);
	if (errCode != 0) {
		cout << "EN_runH, code=" << errCode << endl;
		return;
	}
	for (int i = 1; i <= nodeCount; i++) {//所有节点
		double value;
		EN_getnodevalue(ph, i, EN_PRESSURE, &value);
		initPressures.push_back(value);
		char id[50];
		EN_getnodeid(ph, i, id);
		if (value < 14) {
			cout << "1111 you xiaoyu 14de yali  " << value << " " << id << endl;
		}
	}
	errCode = EN_closeH(ph);
	if (errCode != 0 && errCode != 6) {
		cout << "EN_closeH, code=" << errCode << endl;
		return;
	}
	EN_deleteproject(ph);


	EN_createproject(&ph);
	errCode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
	if (errCode != 0)
	{
		// Call functions that perform desired analysis
		cout << "open inp file err, code =" << errCode << endl;
		return;
	}
	EN_getcount(ph, EN_NODECOUNT, &nodeCount);
	EN_getcount(ph, EN_LINKCOUNT, &linkCount);

	errCode = EN_openH(ph);
	if (errCode != 0) {
		cout << "EN_openH, code=" << errCode << endl;
		return;
	}
	errCode = EN_setdemandmodel(ph, EN_PDA, 10, 40, 0.5); 
	if (errCode != 0) {
		cout << "EN_setdemandmodel, code=" << errCode << endl;
		return;
	}
	for (int i = 0; i < solution.size(); i++) {
		int index;
		char* id = (char*)solution[i].c_str();
		EN_getlinkindex(ph, id, &index);
		
		errCode = EN_setlinkvalue(ph, index, EN_INITSTATUS, EN_CLOSED);
		if (errCode != 0) {
			cout << "641 EN_setlinkvalue, code=" << errCode << "" << ", edgeIndex=" << index << ", edgeID=" << solution[i] << endl;
		}
	}

	errCode = EN_initH(ph, EN_NOSAVE);
	if (errCode != 0) {
		cout << "EN_initH, code=" << errCode << endl;
		return;
	}


	errCode = EN_runH(ph, &t);

	for (int i = 1; i <= nodeCount; i++) {
		double value;
		EN_getnodevalue(ph, i, EN_PRESSURE, &value);
		newPressures.push_back(value);
		char id[50];
		EN_getnodeid(ph, i, id);
		if (value < 14) {
			cout << "2222 you xiaoyu 14de yali  " << value << " " << id << endl;
		}
	}
	errCode = EN_closeH(ph);
	if (errCode != 0) {
		cout << "EN_closeH, code=" << errCode << endl;
		return;
	}
	EN_deleteproject(ph);


	ifstream inPartFile(partFile);
	if (!inPartFile) {
		cout << "打开文件失败2！" << endl;
		exit(1);  
	}
	int count;
	string line;
	getline(inPartFile, line);
	istringstream tmp(line);
	tmp >> count;
	for (int i = 0; i < count; i++) {
		getline(inPartFile, line);
		istringstream tmp(line);
		int index, partition;
		tmp >> index >> partition;
		nodeIndex.push_back(index);
		nodePartition.push_back(partition);
	}

	ofstream outFile(rstFile);
	if (!outFile) {
		cout << "打开文件失败3！" << endl;
		exit(1);
	}
	for (int i = 0; i < nodeIndex.size(); i++) {
		outFile << nodeIndex[i] << "," << nodePartition[i] << "," << initPressures[i] << "," << newPressures[i] << endl;
	}
	outFile.close();

	
}

void run_Test() {
	string inpFile = "E:\\博四下\\metis软件使用\\test\\resource\\inp\\exnet2modified-mopso2.inp"; 
	vector<string> solution = { "04" };
	int nodeCount, linkCount;

	for (int i = 0; i < solution.size(); i++) {
			int errCode;
			long t;

			EN_createproject(&ph);
			errCode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
			if (errCode != 0)
			{
				
				cout << "open inp file err, code =" << errCode << endl;
				return;
			}
			EN_getcount(ph, EN_NODECOUNT, &nodeCount);
			EN_getcount(ph, EN_LINKCOUNT, &linkCount);

			int index1, index2;
			string id1 = "3001", id2 = "3002";
			EN_getnodeindex(ph, (char*)id1.c_str(), &index1);
			EN_getnodeindex(ph, (char*)id2.c_str(), &index2);

			errCode = EN_openH(ph);
			if (errCode != 0) {
				cout << "EN_openH, code=" << errCode << endl;
				return;
			}
			errCode = EN_setdemandmodel(ph, EN_PDA, 18, 50, 0.5);
			if (errCode != 0) {
				cout << "EN_setdemandmodel, code=" << errCode << endl;
				return;
			}
			int index;
			char* id = (char*)solution[i].c_str();
			EN_getlinkindex(ph, id, &index);
			
			errCode = EN_setlinkvalue(ph, index, EN_INITSTATUS, EN_CLOSED);
			if (errCode != 0) {
				cout << "EN_setlinkvalue, code=" << errCode << "" << ", edgeIndex=" << index << ", edgeID=" << solution[i] << endl;
				return;
			}
			
			errCode = EN_initH(ph, EN_NOSAVE);
			if (errCode != 0) {
				cout << "EN_initH, code=" << errCode << endl;
				return;
			}

			
			errCode = EN_runH(ph, &t);
			if (errCode != 0) {
				cout << "EN_runH, code=" << errCode << endl;
				return;
			}
			for (int j = 1; j <= nodeCount; j++) {
				double value;
				EN_getnodevalue(ph, j, EN_PRESSURE, &value);
				if (j == index1 || j == index2) {
					continue;
				}
				if (value < 10) {
					char id[50];
					EN_getnodeid(ph, j, id);
					cout << solution[i] << " " << id << endl;
				}
			}
			errCode = EN_closeH(ph);
			if (errCode != 0) {
				cout << "EN_closeH, code=" << errCode << endl;
				return;
			}

			EN_deleteproject(ph);
		}
	
}

void run_Flowrate() {
	string inpFile = "E:\\学位论文\\DMA结果部分\\test\\resource\\inp\\Modena.inp";
	vector<string> cutedgeID = { "4","27","33","56","82","103","150","154","178","192","205","212","223","242","254","271" };//为了获取这些管段id的flowrateA 
	vector<double> flows;
	EN_Project ph;
	int errCode;
	long t;
	vector<int> linkIndex;
	int nodeCount, linkCount;


	EN_createproject(&ph);
	errCode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
	if (errCode != 0)
	{
		
		cout << "open inp file err, code =" << errCode << endl;
		return;
	}
	EN_getcount(ph, EN_NODECOUNT, &nodeCount);
	EN_getcount(ph, EN_LINKCOUNT, &linkCount);

	errCode = EN_openH(ph);
	if (errCode != 0) {
		cout << "EN_openH, code=" << errCode << endl;
		return;
	}
	errCode = EN_setdemandmodel(ph, EN_PDA, 14, 40, 0.5);
	if (errCode != 0) {
		cout << "EN_setdemandmodel, code=" << errCode << endl;
		return;
	}
	// initialize hydraulics; don't save them to file
	errCode = EN_initH(ph, EN_NOSAVE);
	if (errCode != 0) {
		cout << "EN_initH, code=" << errCode << endl;
		return;
	}

	// solve hydraulics(2.0)
	errCode = EN_runH(ph, &t);
	if (errCode != 0) {
		cout << "EN_runH, code=" << errCode << endl;
		return;
	}
	for (int i = 0; i < cutedgeID.size(); i++) {
		int linkIndex;
		EN_getlinkindex(ph, (char*)cutedgeID[i].c_str(), &linkIndex);
		double value;
		EN_getlinkvalue(ph, linkIndex, EN_FLOW, &value);
		flows.push_back(value);
		cout << "cutedgeID =" << cutedgeID[i] << ",Flowrate=" << value << endl;
	}
	errCode = EN_closeH(ph);
	if (errCode != 0 && errCode != 6) {
		cout << "EN_closeH, code=" << errCode << endl;
		return;
	}
	EN_deleteproject(ph);
}

void run_initialsolution() {
	string inpFile = "";
	vector<string> closedcutedge = { "3636","2443","3730","2218","3744","2173","2716","3479","3489","3564","3253","3166","2365","4894","2819","5271","5241","2534","3309","2141","3549","2757","2700","3181","4919","4153","2904","2586","2717","3037","2734","2899" };//输入阀门的id(0对应的）还要修改skippipe modena "269", "270", "271", "272" Exnet "3001", "3002", "9003", "9004", "9005", "9006", "9007"
	EN_Project ph;
	int errCode;
	long t;
	vector<int> nodeIndex;
	vector<double> newPressures;
	int nodeCount, linkCount;

	
	EN_createproject(&ph);
	errCode = EN_open(ph, (char*)inpFile.c_str(), "resource/report.rpt", "");
	if (errCode != 0)
	{
		// Call functions that perform desired analysis
		cout << "open inp file err, code =" << errCode << endl;
		return;
	}
	EN_getcount(ph, EN_NODECOUNT, &nodeCount);
	EN_getcount(ph, EN_LINKCOUNT, &linkCount);

	errCode = EN_openH(ph);
	if (errCode != 0) {
		cout << "EN_openH, code=" << errCode << endl;
		return;
	}
	errCode = EN_setdemandmodel(ph, EN_PDA, 10, 40, 0.5); 
	if (errCode != 0) {
		cout << "EN_setdemandmodel, code=" << errCode << endl;
		return;
	}
	for (int i = 0; i < closedcutedge.size(); i++) {
		int index;
		char* id = (char*)closedcutedge[i].c_str();
		EN_getlinkindex(ph, id, &index);
	
		errCode = EN_setlinkvalue(ph, index, EN_INITSTATUS, EN_CLOSED);
		if (errCode != 0) {
			cout << "942 EN_setlinkvalue, code=" << errCode << "" << ", edgeIndex=" << index << ", edgeID=" << closedcutedge[i] << endl;
		}
	}
	
	errCode = EN_initH(ph, EN_NOSAVE);
	if (errCode != 0) {
		cout << "EN_initH, code=" << errCode << endl;
		return;
	}

	// solve hydraulics
	errCode = EN_runH(ph, &t);
	/*if (errCode != 0) {
		cout << "EN_runH, code=" << errCode << endl;
		return;
	}*/
	for (int i = 1; i <= nodeCount; i++) {
		double value;
		EN_getnodevalue(ph, i, EN_PRESSURE, &value);
		newPressures.push_back(value);
		char id[50];
		EN_getnodeid(ph, i, id);
		if (value < 10) {
			
			cout << "nodeID=" << id << ",pressure=" << value << endl;
		}
	}
	errCode = EN_closeH(ph);
	if (errCode != 0) {
		cout << "EN_closeH, code=" << errCode << endl;
		return;
	}
	EN_deleteproject(ph);
}

int main() {
	_setmaxstdio(8096);
	run_Split(); 
	
}
