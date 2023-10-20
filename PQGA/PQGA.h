/*
* PQGA
* �����沢�������Ŵ��㷨
*/
#pragma once
#define _USE_MATH_DEFINES 
#define NORMAL 0	  //��ͨ����
#define POP_BEST 1	  //��Ⱥ�����Ÿ���
#define POP_WORST 2	  //��Ⱥ��������
#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include <string>
#include <float.h>
#include <thread>
#include <deque>
using namespace std;

const int PACK_MAX = 3820;					//�����������
const int MAX_GEN = 20;						//����������
const int TABU_GEN = 50;							//������������������
const int POP_SIZE = 200;					//��Ⱥ��С
const int GENE_NUM = 100;					//����������������ĸ���
const int GENE_LEN = 1;					    //ÿ������ı��볤�ȣ���ÿ�������Ķ����Ʊ��볤��
const int CHROM_LEN = GENE_NUM * GENE_LEN;	//����Ķ����Ʊ��볤��
const double MIGRATE_RATE = 0.1;			//�������
const double INIT_AMPLITUDE = 1 / sqrt(2);	//���Ŷ���֮һ���������ڳ�ʼ����Ⱥ
const double PI = M_PI;						//Pi����
const double K1 = 0.001 * PI;				//��С��ת��
const double K2 = 0.05 * PI;				//�����ת��
const double EPSLION = 0.1;					//H-ep��

static vector<double> PACK_WEIGHT;
static vector<double> PACK_VAL;

//���ӱ���
struct qubit {
	double alpha;
	double beta;
	qubit() : alpha(INIT_AMPLITUDE), beta(INIT_AMPLITUDE) {}
	qubit(double a, double b) : alpha(a), beta(b) {}
	string toString() {
		string qubitStr;
		qubitStr += "qubit alpha :" + to_string(alpha) + "\n";
		qubitStr += "qubit beta :" + to_string(beta) + "\n";
		return qubitStr;
	}
};

//�������Χ
struct range {
	double floor;
	double ceil;
	range() : floor(0.0), ceil(0.0) {}
	range(double x, double y) : floor(x), ceil(y) {}
	string toString() {
		return "floor = " + to_string(floor)
			+ "\nceil = " + to_string(ceil)
			+ "\n";
	}
};

class Individual {
private:
	vector<qubit> mChrom;	  //��qubit�����Ⱦɫ��
	double mFitness;		  //��Ӧֵ
	string mBinary;			  //�����Ʊ���
	vector<double> mGenesDec; //ÿ�����򣨱�������ʮ���Ʊ�ʾ
	int mSpecFlag;			  //��¼�������ı�־��0����ͨ���� 1�����Ÿ��� 2�������壩
	bool isPrintQubit = false;
	bool isPrintGeneDec = true;
public:
	Individual() {
		//�޲ι��캯����������ʼ����Ⱥ
		mChrom.resize(CHROM_LEN);
		for (int i = 0; i < CHROM_LEN; i++) {
			qubit initQ = qubit();
			mChrom[i] = initQ;
		}
		this->mFitness = 0.0;
		this->mBinary = "";
		this->mGenesDec.resize(GENE_NUM, 0);
		this->mSpecFlag = 0;
	};
	Individual(vector<qubit> chrom, double fitness, string binary) {
		this->mChrom = chrom;
		this->mFitness = fitness;
		this->mBinary = binary;
		this->mGenesDec.resize(GENE_NUM, 0);
		this->mSpecFlag = 0;
	}
	vector<qubit> getChrom() {
		return this->mChrom;
	}
	void setQubitByPos(int index, double a, double b) {
		qubit q = qubit(a, b);
		this->mChrom[index] = q;
	}
	double getFitness() {
		return this->mFitness;
	}
	void setFitness(double fitness) {
		this->mFitness = fitness;
	}
	vector<double> getGeneDec() {
		return this->mGenesDec;
	}
	void setGeneDec(vector<double> geneDec) {
		this->mGenesDec = geneDec;
	}
	int getSpecFlag() {
		return this->mSpecFlag;
	}
	void setSpecFlag(int specFlag) {
		this->mSpecFlag = specFlag;
	}
	string getBinary() {
		return this->mBinary;
	}
	void setBinary(const string binary) {
		this->mBinary = binary;
	}
	string toString() {
		string ans = "";
		//��ӡqubit
		if (isPrintQubit) {
			for (auto qbit : mChrom) {
				ans += qbit.toString();
			}
		}
		//��ӡʮ��������
		if (isPrintGeneDec) {
			for (int i = 0; i < mGenesDec.size(); i++) {
				ans += "x" + to_string(i) + " = " + to_string(mGenesDec[i]) + "\n";
			}
		}
		//��ӡ��Ӧֵ
		ans += "mFitness = " + to_string(mFitness);
		ans += "\n";
		//��ӡ�����Ʊ���
		ans += "mBinary = " + mBinary + "\n";
		return ans;
	}
};

//���ɱ�item
struct tabuItem {
	pair<int, int> vertexs;	//�Ի�����������λ��
	double fitDiff;			//�Ի������ԭ������������Ӧ��
	Individual newIndv;
	tabuItem() : vertexs(pair<int, int>()), fitDiff(-1), newIndv(Individual()) {}
	tabuItem(pair<int, int> v, double diff, Individual indv) : vertexs(v), fitDiff(diff), newIndv(indv) {}

	bool operator==(const tabuItem& rhs)const {
		if (vertexs.first == rhs.vertexs.first && vertexs.second == rhs.vertexs.second ||
			vertexs.first == rhs.vertexs.second && vertexs.second == rhs.vertexs.first) {
			return true;
		}
		else {
			return false;
		}
	}

	string toString() {
		string ans;
		ans += "i = ";
		ans += to_string(vertexs.first / GENE_LEN);
		ans += '\n';
		ans += "j = ";
		ans += to_string(vertexs.second / GENE_LEN);
		ans += '\n';
		ans += "fit diff = ";
		ans += to_string(fitDiff);
		ans += '\n';
		return ans;
	}
};

double srand();

void initPop(int M, int N, vector<Individual>& population, int strategyFlag);

void initUniverse(vector<Individual>& u1, vector<Individual>& u2, vector<Individual>& u3, vector<Individual>& u4);

void collapse(vector<Individual>& population);

void qGateRAS_1(vector<Individual>& population, Individual& best);

void qGateRAS_2(vector<Individual>& population, Individual& best);

void qGateAdaptive(vector<Individual>& population, Individual& best, double f_max, double f_min);

vector<double> decodeBinary(string binary, vector<range> bound);

double packFunc(Individual& indv, vector<double> weight, vector<double> value);

double objFunc(Individual& indv);

double objFuncShaffer(Individual& indv);

void calFitness(vector<Individual>& population, Individual& best, double f_max, double f_min);

void tabu(vector<Individual>& population, Individual& best);

void migration(vector<Individual>& u1, vector<Individual>& u2, vector<Individual>& u3, vector<Individual>& u4);

void printPopulation(vector<Individual>& pop);

void maQuantumAlgorithm();