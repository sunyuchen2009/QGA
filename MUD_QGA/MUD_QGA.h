/*
* MUD_QGA
* 多用户检测（MultipleUserDetection）场景下的量子遗传算法
*/
#pragma once
#define _USE_MATH_DEFINES 
#define NORMAL 0	    //普通个体
#define POP_BEST 1	    //种群中最优个体
#define POP_WORST 2	    //种群中最差个体
#define GA_FIT 0		//遗传算法适应度计算
#define TS_FIT 1		//禁忌搜索适应度计算
#define UNITARY_FUNC 0
#define SHAFFER_FUNC 1
#define SHUBERT_FUNC 2
#define HUMPBACK_FUNC 3
#define DEJONES_FUNC 4

#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include <string>
#include <float.h>
#include <set>
#include <queue>
#include <deque>
#include <fstream>
#include <unordered_set>
using namespace std;

const int MAX_GEN = 500;							//最大迭代次数
const int TABU_GEN = 50;							//禁忌搜索最大迭代次数
const int POP_SIZE = 500;							//种群大小
const int GENE_NUM = 10;							//基因个数，即变量的个数
const int GENE_LEN = 4;								//每个基因的编码长度，即每个变量的二进制编码长度
const int CHROM_LEN = GENE_NUM * GENE_LEN;			//个体的二进制编码长度
const double MIGRATE_RATE = 0.1;					//移民比率
const double INIT_AMPLITUDE = 1 / sqrt(2);			//根号二分之一常量，用于初始化种群
const double PI = M_PI;								//Pi常量
const double K1 = 0.001 * PI;						//最小旋转角
const double K2 = 0.005 * PI;						//最大旋转角
const double EPSLION = 0.1;							//H-ep门

static vector<double> DISTANCE;						//节点到基站距离
static vector<double> GAINS;						//每个节点到基站的信道增益

//量子比特
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
//个体基因范围
struct range {
	double floor;
	double ceil;
	range(double x, double y) : floor(x), ceil(y) {}
};


class Individual {
private:
	vector<qubit> mChrom;	  //用qubit编码的染色体
	double mFitness;		  //适应值
	string mBinary;			  //二进制编码
	vector<double> mGenesDec; //每个基因（变量）的十进制表示
	int mSpecFlag;			  //记录特殊个体的标志（0：普通个体 1：最优个体 2：最差个体）
	bool isPrintQubit = false;
	bool isPrintGeneDec = false;
public:
	Individual() {
		//无参构造函数，创建初始化种群
		mChrom.resize(CHROM_LEN);
		for (int i = 0; i < CHROM_LEN; i++) {
			qubit initQ = qubit();
			mChrom[i] = initQ;
		}
		this->mFitness = INT_MIN;
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
	//获得当前个体与tarIndv个体之间的海明距离
	int getHammingDis(Individual& tarIndv) {
		string curBin = this->mBinary;
		string tarBin = tarIndv.getBinary();
		int ans = 0;
		for (int i = 0; i < CHROM_LEN; i++) {
			if (curBin[i] != tarBin[i]) {
				ans++;
			}
		}
		return ans;
	}

	virtual ~Individual() {}

	string toString() {
		string ans = "";
		
		//打印qubit
		if (isPrintQubit) {
			for (auto qbit : mChrom) {
				ans += qbit.toString();
			}
		}
		//打印十进制描述
		if (isPrintGeneDec) {
			for (int i = 0; i < mGenesDec.size(); i++) {
				ans += "x" + to_string(i) + " = " + to_string(mGenesDec[i]) + "\n";
			}
		}
		//打印适应值
		ans += "mFitness = " + to_string(mFitness);
		ans += "\n";
		//打印二进制编码
		ans += "mBinary = " + mBinary + "\n";
		return ans;
	}
};

//禁忌表item
struct tabuItem {
	pair<int, int> vertexs;	//对换的两个顶点位置
	double fitDiff;			//对换后相较原个体提升的适应度
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