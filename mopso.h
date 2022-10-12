#include <iostream>
#include <iomanip>
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
#include "helper.h"
#include <random>

using namespace std;


#define N 500					
#define T 3					
#define I 1500                  
#define OBJ_NUM 3				
#define MAX_ARCHIVE_NUM 50	
#define MESH_DIV 10				
#define w 0.2					
#define c1 2					
#define c2 2				
#define Lower 0.35				
#define Upper 0.45				
		
#define minV -50
#define maxV 50

vector<int> pareto(double fitness[][OBJ_NUM], int iNum);
bool dominate(double fitness1[OBJ_NUM], double fitness2[OBJ_NUM]);
double getRand(double min, double max);
double sigmoid(double x);

class Multi_Object_Particle_Swarm_Optimization 
{
private:
	int		i, j;
	int		it;										
	int		Population[N][T];						
	double	Population_Fitness[N][OBJ_NUM];			
	double	Population_V[N][T];						
	double	Sigmoid_V[N][T];						

	int		Population_Copy[N][T];					
	double	Population_Copy_Fitness[N][OBJ_NUM];	

	int     Archive[10000][T];					
	double  Archive_Fitness[10000][OBJ_NUM];		
	int		CurArchiveNum;							
	int		PressureBelowThresholdCount;			

	int		Pbest_Population[N][T];					
	double	Pbest_Fitness[N][OBJ_NUM];				
	int		Gbest_individual[T];					
	double	Gbest_Fitness[OBJ_NUM];					

	bool	Converged;							

	ofstream outStream;

public:
	/*������Init��������ʼ������λ�ã������ٶȣ���Ӧ�ȣ��������ţ��ⲿ�浵��ȫ�����š�*/ //��������λ��ȫ����
	void Init();
	/*������MOPSO����������Ⱥ�㷨���塷*/
	void MOPSO();
	/*������Fitness�������Դ���������Ӧ��ֵ��*/
	vector<double> Fitness(int* input_solution, bool* hasPressureBelowThreshold); //*��ָ���������ŵĵ�ַ
	/*������Update���������������ٶȡ�λ�á���Ӧ�ȡ��������š��ⲿ�浵��ȫ�����š�*/
	void Update();
	
	vector<double> CalculateProbability();
	/*������GetGbest����������ȫ�����š�*/
	void GetGbest();
	
	void Output(int generation);
	
	void WriteResult();
};


void Multi_Object_Particle_Swarm_Optimization::Init()
{
	//��ʼ��epanet�������
	init();
	
	for (i = 0; i < initSolution.size(); i++) {
		for (j = 0; j < T; j++) {
			Population[i][j] = initSolution[i][j];
		}
	}
	for (i = initSolution.size(); i < N; i++) {	//N��Ⱥ��С�����Ӹ�����
		for (j = 0; j < T; j++) { //T ����������
			/*if (skipPipeIDs.find(linkIndexToID[edgeIndex[j]]) != skipPipeIDs.end()) {
				Population[i][j] = 1;
				continue;
			}*/
			/*if (i < T) {
				if (j == i % T) {
					Population[i][j] = 0;
				}
				else {
					Population[i][j] = 1;
				}
				continue;
			}*/
			if (i == 0) {
				Population[i][j] = 1; 
				continue;
			}
			int x = rand() % 2;
			if (x == 1) {
				Population[i][j] = 0;
			}
			else {
				Population[i][j] = 1;//��ά����
			}
		
		}
	}
	
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < T; j++) {
			Population_V[i][j] = minV + rand() / double(RAND_MAX / (maxV - minV));//i�����Ӹ�����j�ǽ���������
		}
	}
	
	PressureBelowThresholdCount = 0;
	for (i = 0; i < N; i++)
	{
		bool hasPressureBelowThreshold = false;
		vector<double> objs = Fitness(Population[i], &hasPressureBelowThreshold);
		for (j = 0; j < OBJ_NUM; j++) {//Ŀ�꺯���ĸ���
			Population_Fitness[i][j] = objs[j];
		}
		if (hasPressureBelowThreshold) {
			PressureBelowThresholdCount++;
		}
	}
	
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < T; j++) {
			Pbest_Population[i][j] = Population[i][j];
		}
	}
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < OBJ_NUM; j++) {
			Pbest_Fitness[i][j] = Population_Fitness[i][j];
		}
	}
	
	vector<int> archiveIndex = pareto(Population_Fitness, N);
	for (i = 0; i < archiveIndex.size(); i++) {
		for (j = 0; j < T; j++) {
			Archive[i][j] = Population[archiveIndex[i]][j];
		}
		for (j = 0; j < OBJ_NUM; j++) {
			Archive_Fitness[i][j] = Population_Fitness[archiveIndex[i]][j];
		} 
	}
	CurArchiveNum = archiveIndex.size();
	
	GetGbest();
	
	Converged = false;
	
	outStream.open("resource/log.txt", ios::out | ios::trunc);
}


void Multi_Object_Particle_Swarm_Optimization::MOPSO()
{
	
	Init();//��ʼ��
	Output(0);
	/*������ʼ*/
	for (it = 1; it <= I; it++)
	{
		cout << "mopso ��������" << it << "�֣��� " << CurArchiveNum << " ���ⲿ����" << endl;
		if (Converged) {
			break;
		}
		/*�ڡ����ú���Update�����������ٶȡ�λ�á���Ӧ�ȡ��������š��ⲿ�浵��ȫ�����š�*/
		Update();//����
		Output(it);
	}
	WriteResult();
	outStream.close();
	closeEpanet();
}

/*������Update���壺���������ٶȡ�λ�á���Ӧ�ȡ��������š��ⲿ�浵��ȫ�����š�*/
void Multi_Object_Particle_Swarm_Optimization::Update()
{
	
	for (i = 0; i < N; i++) {
		for (j = 0; j < T; j++) { 
			Population_V[i][j] = w * Population_V[i][j] + c1 * getRand(Lower, Upper) * (double(Pbest_Population[i][j]) - Population[i][j]) + c2 * getRand(Lower, Upper) * (double(Gbest_individual[j]) - Population[i][j]);
			if (Population_V[i][j] > maxV) {
				Population_V[i][j] = maxV;
			}
			if (Population_V[i][j] < minV) {
				Population_V[i][j] = minV;
			}
			Sigmoid_V[i][j] = sigmoid(Population_V[i][j]);
		}
	}

	
	for (i = 0; i < N; i++) {
		for (j = 0; j < T; j++) {
			
			double r = rand() / (RAND_MAX * 1.0);
			if (Sigmoid_V[i][j] <= r) {    
			Population_Copy[i][j] = 1;
			}
			else {
				Population_Copy[i][j] = 0;
			}
			
		}
	}

	//����������Ӧ�ȣ�Population_Copy_Fitness
	PressureBelowThresholdCount = 0;
	for (i = 0; i < N; i++)
	{
		bool hasPressureBelowThreshold = false;
		vector<double> objs = Fitness(Population_Copy[i], &hasPressureBelowThreshold);
		for (j = 0; j < OBJ_NUM; j++) {
			Population_Copy_Fitness[i][j] = objs[j];
		}
		if (hasPressureBelowThreshold) {
			PressureBelowThresholdCount++;
		}
	}
	
	for (i = 0; i < N; i++) {//����Ⱥ����
		
		bool oldDominateNew = dominate(Population_Fitness[i], Population_Copy_Fitness[j]);
		
		bool newDominateOld = dominate(Population_Fitness[j], Population_Copy_Fitness[i]);
	
		if (oldDominateNew) {
			for (j = 0; j < T; j++) {
				Pbest_Population[i][j] = Population[i][j];
			}
		}
	
		if (newDominateOld) {
			for (j = 0; j < T; j++) {
				Pbest_Population[i][j] = Population_Copy[i][j];  
			}
		}
		
		if (!oldDominateNew && !newDominateOld) {
			double r = rand() / (RAND_MAX * 1.0);
			if (r >= 0.5) {
				for (j = 0; j < T; j++) {
					Pbest_Population[i][j] = Population_Copy[i][j]; //pbest�Ǿֲ�����λ��
				}
			}
			else {
				for (j = 0; j < T; j++) {
					Pbest_Population[i][j] = Population[i][j];
				}
			}
		}
	}
	for (i = 0; i < N; i++) {
		for (j = 0; j < T; j++) {
			Population[i][j] = Population_Copy[i][j];
		}
		for (j = 0; j < OBJ_NUM; j++) {
			Population_Fitness[i][j] = Population_Copy_Fitness[i][j];
		}
	}

	
	vector<int> dominateIndex = pareto(Population_Fitness, N);
	for (i = 0; i < dominateIndex.size(); i++) {
		for (j = 0; j < T; j++) {
			Archive[CurArchiveNum + i][j] = Population[dominateIndex[i]][j];
		}
		for (j = 0; j < OBJ_NUM; j++) {
			Archive_Fitness[CurArchiveNum + i][j] = Population_Fitness[dominateIndex[i]][j];
		}
	}
	CurArchiveNum += dominateIndex.size();
	
	dominateIndex = pareto(Archive_Fitness, CurArchiveNum);
	for (i = 0; i < dominateIndex.size(); i++) {
		for (j = 0; j < T; j++) {
			Archive[i][j] = Archive[dominateIndex[i]][j];
		}
		for (j = 0; j < OBJ_NUM; j++) {
			Archive_Fitness[i][j] = Archive_Fitness[dominateIndex[i]][j];
		}
	}
	CurArchiveNum = dominateIndex.size();
	
	if (CurArchiveNum > MAX_ARCHIVE_NUM) {
		vector<double> Probability = CalculateProbability();
		vector<int> indexes(Probability.size(), 0);
		for (int i = 0; i != indexes.size(); i++) {
			indexes[i] = i;
		}
		sort(indexes.begin(), indexes.end(),
			[&](const int& a, const int& b) {
			return (Probability[a] > Probability[b]);
		});
		int Archive_Copy[MAX_ARCHIVE_NUM][T];
		double Archive_Copy_Fitness[MAX_ARCHIVE_NUM][OBJ_NUM];
		for (i = 0; i < MAX_ARCHIVE_NUM; i++) {
			int index = indexes[i];
			for (j = 0; j < T; j++) {
				Archive_Copy[i][j] = Archive[index][j];
			}
			for (j = 0; j < OBJ_NUM; j++) {
				Archive_Copy_Fitness[i][j] = Archive_Fitness[index][j];
			}
		}
		for (i = 0; i < MAX_ARCHIVE_NUM; i++) {
			for (j = 0; j < T; j++) {
				Archive[i][j] = Archive_Copy[i][j];
			}
			for (j = 0; j < OBJ_NUM; j++) {
				Archive_Fitness[i][j] = Archive_Copy_Fitness[i][j];
			}
		}
		CurArchiveNum = MAX_ARCHIVE_NUM;
	}

	
	GetGbest();
	
	Converged = false;
}

vector<double> Multi_Object_Particle_Swarm_Optimization::Fitness(int* input_solution, bool* hasPressureBelowThreshold)
{
	vector<int> solution;
	for (int j = 0; j < T; j++) {//T�ǽ���������
		solution.push_back(input_solution[j]);
	}
	return calculateObjectives(solution, hasPressureBelowThreshold);
}

/*��CalculateProbability���壺�����ⲿ�������ĸ��ʣ�ӵ����ԽС������Խ�󣩡�*/
vector<double> Multi_Object_Particle_Swarm_Optimization::CalculateProbability() {
	double Max[OBJ_NUM] = { DBL_MIN };
	double Min[OBJ_NUM] = { DBL_MAX };
	for (i = 0; i < CurArchiveNum; i++) {
		for (j = 0; j < OBJ_NUM; j++) {
			Max[j] = max(Max[j], Archive_Fitness[i][j]);
			Min[j] = max(Min[j], Archive_Fitness[i][j]);
		}
	}
	double Mod[OBJ_NUM];
	for (i = 0; i < OBJ_NUM; i++) {
		Mod[i] = (Max[i] - Min[i]) / MESH_DIV;
	}
	map<string, vector<int>> pos2IndexMap;
	for (i = 0; i < CurArchiveNum; i++) {
		string key = "";
		for (j = 0; j < OBJ_NUM; j++) {
			if (j != 0) {
				key += "#";
			}
			key += to_string(floor((Archive_Fitness[i][j] - Min[j]) / Mod[j]));
		}
		pos2IndexMap[key].push_back(i);
	}
	vector<int> count(CurArchiveNum);
	for (auto iter = pos2IndexMap.begin(); iter != pos2IndexMap.end(); iter++) {
		int c = iter->second.size();
		for (i = 0; i < c; i++) {
			count[iter->second[i]] = c;
		}
	}
	vector<double> Probability(CurArchiveNum);
	double ProbabilitySum = 0;
	for (i = 0; i < CurArchiveNum; i++) {
		Probability[i] = 1 / double(pow(count[i], 3));
		ProbabilitySum += Probability[i];
	}
	for (i = 0; i < CurArchiveNum; i++) {
		Probability[i] /= ProbabilitySum;
	}
	return Probability;
}

/*��GetGbest�������壺����ȫ�����š�*/
void Multi_Object_Particle_Swarm_Optimization::GetGbest()
{
	vector<double> Probability = CalculateProbability();
	double r = rand() / (RAND_MAX * 1.0);
	double Probability_Sum = Probability[0];
	int selected = 0;
	for (i = 1; i < CurArchiveNum; i++) {
		if (Probability_Sum <= r && r < Probability_Sum + Probability[i]) {
			selected = i - 1;
		}
	}
	for (i = 0; i < T; i++) {
		Gbest_individual[i] = Archive[selected][i];
	}
	for (i = 0; i < OBJ_NUM; i++) {
		Gbest_Fitness[i] = Archive_Fitness[selected][i];
	}
}

void Multi_Object_Particle_Swarm_Optimization::Output(int generation) {
	if (generation == I) {
		cout << "���һ�֣�CurArchiveNum=" << CurArchiveNum << endl;
	}
	if (generation == 0) {
		outStream << "��ʼ�⣨��" << CurArchiveNum << "���ⲿ������" << "����" << PressureBelowThresholdCount << "�������ѹ��С����ֵ����" << endl;
	}
	else {
		outStream << "��" << generation << "���⣨��" << CurArchiveNum << "���ⲿ������" << "����" << PressureBelowThresholdCount << "�������ѹ��С����ֵ����" << endl;
	}
	for (i = 0; i < CurArchiveNum; i++) {
		outStream << "����" << "��[";
		for (j = 0; j < T; j++) {
			outStream << Archive[i][j];
			if (j == T - 1) {
				outStream << "]";
			}
			else {
				outStream << ",";
			}
		}
		outStream << "����Ӧ�ȣ�[";
		for (j = 0; j < OBJ_NUM; j++) {
			outStream << setiosflags(ios::fixed) << setprecision(6) << Archive_Fitness[i][j];
			if (j == OBJ_NUM - 1) {
				outStream << "]";
			}
			else {
				outStream << ",";
			}
		}
		outStream << endl;
	}
	outStream << endl;
}

void Multi_Object_Particle_Swarm_Optimization::WriteResult() {
	ofstream resStream;
	resStream.open("resource/result.csv", ios::out | ios::trunc);
	for (i = 0; i < CurArchiveNum; i++) {
		for (j = 0; j < OBJ_NUM; j++) {
			double data = Archive_Fitness[i][j];
			if (j == 1) {
				data = -data;
			}
			resStream << setiosflags(ios::fixed) << setprecision(6) << data;
			if (j != OBJ_NUM - 1) {
				resStream << ",";
			}
		}
		if (i != CurArchiveNum - 1) {
			resStream << endl;
		}
	}
	resStream.close();
}

vector<int> pareto(double fitness[][OBJ_NUM], int iNum) {//pareto�߽� iNum������
	vector<int> rst;
	for (int i = 0; i < iNum; i++) {
		bool dominated;//bool����true or false
		for (int j = 0; j < iNum; j++) {
			if (i == j) {
				continue;
			}
			dominated = dominate(fitness[j], fitness[i]);
			// ����true��ʾ��j֧���i
			if (dominated) {
				break;
			}
		}
		if (!dominated) {
			rst.push_back(i);
		}
	}
	return rst;
}

// �ж�fitness1�Ƿ�֧��fitness2
bool dominate(double fitness1[OBJ_NUM], double fitness2[OBJ_NUM]) {
	int count1 = 0, count2 = 0, count3 = 0;
	for (int k = 0; k < OBJ_NUM; k++) {
		if (fitness1[k] < fitness2[k]) {
			count1++;
		}
		else if (fitness1[k] == fitness2[k]) {
			count2++;
		}
		else {
			count3++;
		}
	}
	if (count1 + count2 == OBJ_NUM && count1 > 0) {
		return true;
	}
	return false;
}

double getRand(double min, double max)
{
	return min + rand() / double(RAND_MAX / (max - min));
}

double sigmoid(double x)
{
	return (1 / (1 + exp(-x)));
}