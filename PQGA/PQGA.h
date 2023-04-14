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
using namespace std;

const int MAX_GEN = 50;						//����������
const int POP_SIZE = 100;					//��Ⱥ��С
const int GENE_NUM = 2;						//����������������ĸ���
const int GENE_LEN = 21;					//ÿ������ı��볤�ȣ���ÿ�������Ķ����Ʊ��볤��
const int CHROM_LEN = GENE_NUM * GENE_LEN;	//����Ķ����Ʊ��볤��
const double MIGRATE_RATE = 0.1;			//�������
const double INIT_AMPLITUDE = 1 / sqrt(2);	//���Ŷ���֮һ���������ڳ�ʼ����Ⱥ
const double PI = M_PI;						//Pi����
const double K1 = 0.001 * PI;				//��С��ת��
const double K2 = 0.05 * PI;				//�����ת��
const double EPSLION = 0.1;					//H-ep��

struct qubit {
	double alpha;
	double beta;
	qubit() : alpha(INIT_AMPLITUDE), beta(INIT_AMPLITUDE) {}
	qubit(double a, double b) : alpha(a), beta(b) {}
};
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
		this->mGenesDec = { 0, 0 };
		this->mSpecFlag = 0;
	};
	Individual(vector<qubit> chrom, double fitness, string binary) {
		this->mChrom = chrom;
		this->mFitness = fitness;
		this->mBinary = binary;
		this->mGenesDec = { 0, 0 };
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
		////��ӡ��
		//ans += "qubit alpha :";
		//for (qubit q : mChrom) {
		//	ans += to_string(q.alpha);
		//}
		//ans += "\n";
		////��ӡ��
		//ans += "qubit beta :";
		//for (qubit q : mChrom) {
		//	ans += to_string(q.beta);
		//}
		for (int i = 0; i < mGenesDec.size(); i++) {
			ans += "x" + to_string(i) + " = " + to_string(mGenesDec[i]) + "\n";
		}
		//��ӡ��Ӧֵ
		ans += "mFitness = " + to_string(mFitness);
		ans += "\n";
		//��ӡ�����Ʊ���
		ans += "mBinary = " + mBinary + "\n";
		return ans;
	}
};