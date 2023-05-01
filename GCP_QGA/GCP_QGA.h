#pragma once
#define _USE_MATH_DEFINES 
#define NORMAL 0	  //��ͨ����
#define POP_BEST 1	  //��Ⱥ�����Ÿ���
#define POP_WORST 2	  //��Ⱥ��������
#define GA_FIT 0
#define TS_FIT 1

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

const int MAX_GEN = 5000;							//����������
const int TABU_GEN = 10;							//������������������
const int POP_SIZE = 200;							//��Ⱥ��С
const int GENE_NUM = 2;							//����������������ĸ���
const int GENE_LEN = 25;								//ÿ������ı��볤�ȣ���ÿ�������Ķ����Ʊ��볤��
const int CHROM_LEN = GENE_NUM * GENE_LEN;			//����Ķ����Ʊ��볤��
const double MIGRATE_RATE = 0.1;					//�������
const double INIT_AMPLITUDE = 1 / sqrt(2);			//���Ŷ���֮һ���������ڳ�ʼ����Ⱥ
const double PI = M_PI;								//Pi����
const double K1 = 0.001 * PI;						//��С��ת��
const double K2 = 0.008 * PI;						//�����ת��
const double EPSLION = 0.1;							//H-ep��

static vector<pair<int, int>> GCP_EDGE;				//����ͼ�߼���

//���ӱ���
struct qubit {
	double alpha;
	double beta;
	qubit() : alpha(INIT_AMPLITUDE), beta(INIT_AMPLITUDE) {}
	qubit(double a, double b) : alpha(a), beta(b) {}
};
//�������Χ
struct range {
	double floor;
	double ceil;
	range(double x, double y) : floor(x), ceil(y) {}
};


class Individual {
private:
	vector<qubit> mChrom;	  //��qubit�����Ⱦɫ��
	double mFitness;		  //��Ӧֵ
	string mBinary;			  //�����Ʊ���
	vector<double> mGenesDec; //ÿ�����򣨱�������ʮ���Ʊ�ʾ
	int mSpecFlag;			  //��¼�������ı�־��0����ͨ���� 1�����Ÿ��� 2�������壩
	int mSameCnt;			  //��ǰ����������ڶ�������ͬ��ɫ�Ķ�����Ϊ0���ǿ��ý�
	int mColorNum;			  //��ǰ�������ʹ����ɫ������ԽСԽ��
public:
	Individual() {
		//�޲ι��캯����������ʼ����Ⱥ
		mChrom.resize(CHROM_LEN);
		for (int i = 0; i < CHROM_LEN; i++) {
			qubit initQ = qubit();
			mChrom[i] = initQ;
		}
		this->mFitness = INT_MIN;
		this->mBinary = "";
		this->mGenesDec.resize(GENE_NUM, 0);
		this->mSpecFlag = 0;
		this->mSameCnt = -1;
		this->mColorNum = -1;
	};
	Individual(vector<qubit> chrom, double fitness, string binary) {
		this->mChrom = chrom;
		this->mFitness = fitness;
		this->mBinary = binary;
		this->mGenesDec.resize(GENE_NUM, 0);
		this->mSpecFlag = 0;
		this->mSameCnt = -1;
		this->mColorNum = -1;
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
	int getSameCnt() {
		return this->mSameCnt;
	}
	void setSameCnt(int sameCnt) {
		this->mSameCnt = sameCnt;
	}
	int getColorNum() {
		return this->mColorNum;
	}
	void setColorNum(int colorNum) {
		this->mColorNum = colorNum;
	}
	string getBinary() {
		return this->mBinary;
	}
	void setBinary(const string binary) {
		this->mBinary = binary;
	}
	//��õ�ǰ������tarIndv����֮��ĺ�������
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
	string toString() {
		string ans = "";
		////��ӡ��
		//ans += "qubit alpha :";
		//for (qubit q : mChrom) {
		//	ans += to_string(q.alpha) + " ";
		//}
		//ans += "\n";
		////��ӡ��
		//ans += "qubit beta :";
		//for (qubit q : mChrom) {
		//	ans += to_string(q.beta) + " ";
		//}
		ans += '\n';
		for (int i = 0; i < mGenesDec.size(); i++) {
			ans += "x" + to_string(i) + " = " + to_string(mGenesDec[i]) + "\n";
		}
		//��ӡ��Ӧֵ
		ans += "mFitness = " + to_string(mFitness);
		ans += "\nmSameCnt = " + to_string(mSameCnt);
		ans += "\nmColorNum = " + to_string(mColorNum);
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