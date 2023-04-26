#include "GCP_QGA.h"
using namespace std;

//生成[0, 1]之间的随机数
double srand() {
    int N = rand() % 999;
    double random = static_cast<double>(N) / 1000.0;;//随机产生0到1的小数
    //cout << random << endl;
    return random;
}

/*
* 初始化种群，使用量子比特编码
* @param M 种群大小
* @param N 每个基因的量子比特编码长度
* @param strategyFlag 初始化策略：0（默认两个量子态等概率）
                                  1（小生境初始化策略） (α, β) = (sqrt(j / n), sqrt(1 - j / n))
*/
void initPop(int M, int N, vector<Individual>& population, int strategyFlag) {
    switch (strategyFlag)
    {
    case 0:
        //将种群中所有的染色体初始化为根号二分之一
        population.resize(M);
        for (int i = 0; i < M; i++) {
            population[i] = Individual();
        }
        break;
    case 1: //小生境
        population.resize(M);
        for (int i = 0; i < M; i++) {
            vector<qubit> chrom(CHROM_LEN);
            for (int j = 0; j < CHROM_LEN; j++) {
                qubit initQ = qubit(sqrt(j / CHROM_LEN), sqrt(1 - j / CHROM_LEN));
                chrom[j] = initQ;
            }
            population[i] = Individual(chrom, 0, "");
        }
        break;
    default:
        break;
    }
}

/*
* 塌缩函数，完成对种群的一次测量
*/
void collapse(vector<Individual>& population) {
    for (int i = 0; i < POP_SIZE; i++) {
        string binTemp = "";
        for (int j = 0; j < CHROM_LEN; j++) {
            double pick = srand();
            double alpha = population[i].getChrom()[j].alpha;
            //cout << "pick = " << pick << '\n';
            //cout << "alpha^2 = " << alpha * alpha << '\n';

            if (pick > alpha * alpha) {
                binTemp += '1';
            }
            else {
                binTemp += '0';
            }
        }
        population[i].setBinary(binTemp);
        //同时初始化specFlag;
        population[i].setSpecFlag(0);
        //population[i].setBinary("1111011111110100101111110100101101100110");
    }
}

/*
* 量子旋转门QGate RAS_1
*   x	best 	f(x)>=f(best)	delta	    ab>0	ab<0	a=0	    b=0
*   0	0		false			0		    0		0		0	    0
*   0	0		true			0		    0		0		0	    0
*   0	1		false			0   	    0		0		0	    0
*   0	1		true			0.01pi	    -1		1		+-1	    0
*   1	0		false			0.01pi	    -1		1		+-1	    0
*   1	0		true			0.01pi	    1		-1		0	    +-1
*   1	1		false			0.01pi	    1		-1		0	    +-1
*   1	1		true			0.01pi	    1		-1		0	    +-1
*/
void qGateRAS_1(vector<Individual>& population, Individual& best) {
    double bestFit = best.getFitness();                         //最优个体适应度
    for (int i = 0; i < POP_SIZE; i++) {
        for (int j = 0; j < CHROM_LEN; j++) {
            double alpha = population[i].getChrom()[j].alpha;   //α
            double beta = population[i].getChrom()[j].beta;     //β
            char x = population[i].getBinary()[j];              //当前个体第j个编码
            char b = best.getBinary()[j];                       //最优个体第j个编码
            double delta = 0.0;                                 //旋转角大小
            int s = 0;                                          //旋转角方向，即正负号
            double curFit = population[i].getFitness();         //当前个体适应度
            //当前个体与最优个体相同，不旋转
            if ((x == '0' && b == '0')) {
                delta = 0;
                s = 0;
            }
            else if ((x == '0' && b == '1') && (curFit < bestFit)) {
                delta = 0;
                s = 0;
            }
            else if (x == '0' && b == '1' && curFit >= bestFit) {
                delta = 0.01 * PI;
                if (alpha * beta > 0) {
                    s = -1;
                }
                else if (alpha * beta < 0) {
                    s = 1;
                }
                else if (alpha == 0) {
                    s = (srand() > 0.5) ? 1 : -1;
                }
                else if (beta == 0) {
                    s = 0;
                }
            }
            else if (x == '1' && b == '0' && curFit < bestFit) {
                delta = 0.01 * PI;
                if (alpha * beta > 0) {
                    s = -1;
                }
                else if (alpha * beta < 0) {
                    s = 1;
                }
                else if (alpha == 0) {
                    s = (srand() > 0.5) ? 1 : -1;
                }
                else if (beta == 0) {
                    s = 0;
                }
            }
            else if (x == '1' && b == '0' && curFit >= bestFit) {
                delta = 0.01 * PI;
                if (alpha * beta > 0) {
                    s = 1;
                }
                else if (alpha * beta < 0) {
                    s = -1;
                }
                else if (alpha == 0) {
                    s = 0;
                }
                else if (beta == 0) {
                    s = (srand() > 0.5) ? 1 : -1;
                }
            }
            else if (x == '1' && b == '1' && curFit < bestFit) {
                delta = 0.01 * PI;
                if (alpha * beta > 0) {
                    s = 1;
                }
                else if (alpha * beta < 0) {
                    s = -1;
                }
                else if (alpha == 0) {
                    s = 0;
                }
                else if (beta == 0) {
                    s = (srand() > 0.5) ? 1 : -1;
                }
            }
            else if (x == '1' && b == '1' && curFit >= bestFit) {
                delta = 0.01 * PI;
                if (alpha * beta > 0) {
                    s = 1;
                }
                else if (alpha * beta < 0) {
                    s = -1;
                }
                else if (alpha == 0) {
                    s = 0;
                }
                else if (beta == 0) {
                    s = (srand() > 0.5) ? 1 : -1;
                }
            }
            double e = s * delta;       //旋转角
            double newAlpha = alpha * cos(e) - beta * sin(e);
            double newBeta = alpha * sin(e) + beta * cos(e);
            population[i].setQubitByPos(j, newAlpha, newBeta);
        }
    }
}

/*
* 量子旋转门QGate RAS_2
*   x	best 	f(x)>f(best)	delta	    ab>0	ab<0	a=0	    b=0
*   0	0		false			0		    0		0		0	    0
*   0	0		true			0		    0		0		0	    0
*   0	1		false			0.01pi	    1		-1		0	    +-1
*   0	1		true			0.01pi	    -1		1		+-1	    0
*   1	0		false			0.01pi	    -1		1		+-1	    0
*   1	0		true			0.01pi	    1		-1		0	    +-1
*   1	1		false			0		    0		0		0	    0
*   1	1		true			0		    0		0		0	    0
*/
void qGateRAS_2(vector<Individual>& population, Individual& best) {
    double bestFit = best.getFitness();                         //最优个体适应度
    for (int i = 0; i < POP_SIZE; i++) {
        for (int j = 0; j < CHROM_LEN; j++) {
            double alpha = population[i].getChrom()[j].alpha;   //α
            double beta = population[i].getChrom()[j].beta;     //β
            char x = population[i].getBinary()[j];              //当前个体第j个编码
            char b = best.getBinary()[j];                       //最优个体第j个编码
            double delta = 0.0;                                 //旋转角大小
            int s = 0;                                          //旋转角方向，即正负号
            double curFit = population[i].getFitness();         //当前个体适应度
            //当前个体与最优个体相同，不旋转
            if ((x == '0' && b == '0') || (x == '1' && b == '1')) {
                delta = 0;
                s = 0;
            }
            else if ((x == '0' && b == '1') && (curFit < bestFit)) {
                delta = 0.01 * PI;
                if (alpha * beta > 0) {
                    s = 1;
                }
                else if (alpha * beta < 0) {
                    s = -1;
                }
                else if (alpha == 0) {
                    s = 0;
                }
                else if (beta == 0) {
                    s = (srand() > 0.5) ? 1 : -1;
                }
            }
            else if (x == '0' && b == '1' && curFit >= bestFit) {
                delta = 0.01 * PI;
                if (alpha * beta > 0) {
                    s = -1;
                }
                else if (alpha * beta < 0) {
                    s = 1;
                }
                else if (alpha == 0) {
                    s = (srand() > 0.5) ? 1 : -1;
                }
                else if (beta == 0) {
                    s = 0;
                }
            }
            else if (x == '1' && b == '0' && curFit < bestFit) {
                delta = 0.01 * PI;
                if (alpha * beta > 0) {
                    s = -1;
                }
                else if (alpha * beta < 0) {
                    s = 1;
                }
                else if (alpha == 0) {
                    s = (srand() > 0.5) ? 1 : -1;
                }
                else if (beta == 0) {
                    s = 0;
                }
            }
            else if (x == '1' && b == '0' && curFit >= bestFit) {
                delta = 0.01 * PI;
                if (alpha * beta > 0) {
                    s = 1;
                }
                else if (alpha * beta < 0) {
                    s = -1;
                }
                else if (alpha == 0) {
                    s = 0;
                }
                else if (beta == 0) {
                    s = (srand() > 0.5) ? 1 : -1;
                }
            }
            double e = s * delta;       //旋转角
            double newAlpha = alpha * cos(e) - beta * sin(e);
            double newBeta = alpha * sin(e) + beta * cos(e);
            population[i].setQubitByPos(j, newAlpha, newBeta);
        }
    }
}

/*
* 量子旋转门QGate Adaptive
*   x	best 	f(x)>=f(best)	delta	    ab>0	ab<0	a=0	    b=0
*   0	0		false			0		    0		0		0	    0
*   0	0		true			0		    0		0		0	    0
*   0	1		false			Θ   	    1		-1		0	    +-1
*   0	1		true			Θ   	    -1		1		+-1	    0
*   1	0		false			Θ   	    -1		1		+-1	    0
*   1	0		true			Θ   	    1		-1		0	    +-1
*   1	1		false			0		    0		0		0	    0
*   1	1		true			0		    0		0		0	    0
* 其中，旋转角Θ的计算公式如下:
* Θ = K1 + (K2 - K1) * (f_i - f_min) / (f_max - f_min)  [f_min != f_max]
* Θ = K1                                                [f_min == f_max]
* f_i、f_max、f_min分别戴代表当前个体的适应度、当前种群最优个体适应度以及最差适应度
* K1、K2为两个正常数且K1 < K2，用于控制收敛速度
*/
void qGateAdaptive(vector<Individual>& population, Individual& best, double& f_max, double& f_min) {
    const double CONST_ARG = (f_max > f_min) ? (K2 - K1) / (f_max - f_min) : 0;

    double bestFit = best.getFitness();                         //最优个体适应度
    for (int i = 0; i < POP_SIZE; i++) {
        //根据当前个体计算旋转角
        double theta = K1 + (population[i].getFitness() - f_min) * CONST_ARG;
        //个体变异概率，默认为80%
        double mutationPick = srand();
        if (mutationPick > 0.8) {
            continue;
        }
        for (int j = 0; j < CHROM_LEN; j++) {
            double alpha = population[i].getChrom()[j].alpha;   //α
            double beta = population[i].getChrom()[j].beta;     //β
            double delta = 0.0;                                 //旋转角大小
            char x = population[i].getBinary()[j];              //当前个体第j个编码
            char b = best.getBinary()[j];                       //最优个体第j个编码
            int s = 0;                                          //旋转角方向，即正负号
            double curFit = population[i].getFitness();         //当前个体适应度
            //当前个体与最优个体相同，不旋转
            if ((x == '0' && b == '0') || (x == '1' && b == '1')) {
                delta = 0;
                s = 0;
            }
            else if ((x == '0' && b == '1') && (curFit < bestFit)) {
                delta = theta;
                if (alpha * beta > 0) {
                    s = 1;
                }
                else if (alpha * beta < 0) {
                    s = -1;
                }
                else if (alpha == 0) {
                    s = 0;
                }
                else if (beta == 0) {
                    s = (srand() > 0.5) ? 1 : -1;
                }
            }
            else if (x == '0' && b == '1' && curFit >= bestFit) {
                delta = theta;
                if (alpha * beta > 0) {
                    s = -1;
                }
                else if (alpha * beta < 0) {
                    s = 1;
                }
                else if (alpha == 0) {
                    s = (srand() > 0.5) ? 1 : -1;
                }
                else if (beta == 0) {
                    s = 0;
                }
            }
            else if (x == '1' && b == '0' && curFit < bestFit) {
                delta = theta;
                if (alpha * beta > 0) {
                    s = -1;
                }
                else if (alpha * beta < 0) {
                    s = 1;
                }
                else if (alpha == 0) {
                    s = (srand() > 0.5) ? 1 : -1;
                }
                else if (beta == 0) {
                    s = 0;
                }
            }
            else if (x == '1' && b == '0' && curFit >= bestFit) {
                delta = theta;
                if (alpha * beta > 0) {
                    s = 1;
                }
                else if (alpha * beta < 0) {
                    s = -1;
                }
                else if (alpha == 0) {
                    s = 0;
                }
                else if (beta == 0) {
                    s = (srand() > 0.5) ? 1 : -1;
                }
            }
            double e = s * delta;       //旋转角
            /*cout << "-----------------------\n";
            cout << "old alpha = " << alpha << '\n';
            cout << "old beta = " << beta << '\n';
            cout << "cur fit = " << curFit << '\n';
            cout << "best fit = " << bestFit << '\n';
            cout << "x = " << x << '\n';
            cout << "b = " << b << '\n';
            cout << "s = " << s << '\n';
            cout << "a * b = " << alpha * beta << '\n';*/
            
            if (e == 0) {
                //H-ep门
                if (alpha * alpha <= EPSLION && beta * beta >= 1 - EPSLION) {
                    population[i].setQubitByPos(j, sqrt(EPSLION), sqrt(1 - EPSLION));
                }
                else if (beta * beta <= EPSLION && alpha * alpha >= 1 - EPSLION) {
                    population[i].setQubitByPos(j, sqrt(1 - EPSLION), sqrt(EPSLION));
                }
            }
            else {
                double newAlpha = alpha * cos(e) - beta * sin(e);
                double newBeta = alpha * sin(e) + beta * cos(e);
                population[i].setQubitByPos(j, newAlpha, newBeta);
            }
            /*cout << "new alpha = " << population[i].getChrom()[j].alpha << '\n';
            cout << "new beta = " << population[i].getChrom()[j].beta << '\n';
            cout << "-----------------------\n";*/
        }
    }

}

/*
* 解码二进制编码，返回十进制，并转移到变量所限制的区间
* @param binary: 个体的二进制编码，包括多个基因
* @param bound: 变量的上下限
* @return: 各个基因（变量）的数组
*/
vector<double> decodeBinary(string binary, vector<range> bound) {
    vector<double> x(GENE_NUM);
    //binary = "1111011111110100101111110100101101100110";
    for (int i = 0; i < binary.size(); i += GENE_LEN) {
        int index = i / GENE_LEN;
        string curGene = binary.substr(i, GENE_LEN);
        double var = stoi(curGene, nullptr, 2);
        //将变量值映射到对应的range区间
        var = bound[index].floor + var / (pow(2, GENE_LEN) - 1) * (bound[index].ceil - bound[index].floor);
        x[index] = static_cast<double>(var);
        //cout << "x" << index << " = " << x[index] << endl;
    }
    return x;
}

/*
* 待优化的目标函数如下：
* f(x, y) = x * sin(4 * Pi * x) + y * sin(20 * Pi * y)
* x:[-3.0, 12.1]
* y:[4.1, 5.8]
* @param binary：个体的二进制编码
*/
double objFunc(Individual& indv) {
    string binary = indv.getBinary();
    vector<double> x(GENE_NUM); //各个变量的数组
    vector<range> bound = { range(-3.0, 12.1), range(4.1, 5.8) };
    x = decodeBinary(binary, bound);
    double fitness = x[0] * sin(4 * PI * x[0]) + x[1] * sin(20 * PI * x[1]);
    indv.setGeneDec(x);
    return fitness;
}

/*
* Shaffer's Func:
* 此函数有无限多个局部极大值点，其中只有一个(0, 0)为全局最大
* x: [-10, 10]
* y: [-10, 10]
*/
double objFuncShaffer(Individual& indv) {
    string binary = indv.getBinary();
    vector<double> x(GENE_NUM); //各个变量的数组
    vector<range> bound = { range(-10, 10), range(-10, 10) };
    x = decodeBinary(binary, bound);
    double fitness = 0.5 - (pow(sin(sqrt(x[0] * x[0] + x[1] * x[1])), 2) - 0.5) / (pow(1 + 0.001 * (x[0] * x[0] + x[1] * x[1]), 2));
    indv.setGeneDec(x);
    return fitness;
}

/*
* Gragh Coloring Problem Fitness:
* 每个基因代表了当前顶点的颜色编号（从0  - 10）
*/
double gcpFunc(Individual& indv) {
    string binary = indv.getBinary();
    vector<double> x(GENE_NUM, 0);
    vector<range> bound = { GENE_NUM, range(0, 5)};
    set<double> colors;
    //vector<pair<int, int>> edge = { {0, 1}, {0, 4}, {0, 5}, {1, 2},{1, 0}, {1, 6}, {2, 1}, {2, 3}, {2, 7}, {3, 2}, {3, 4}, {3, 8}, {4, 3}, {4, 8}, {4, 9}, {5, 0}, {5, 7}, {5, 8}, {6, 1}, {6, 8}, {6, 9}, {7, 2}, {7, 5}, {7, 9}, {8, 3}, {8, 5}, {8, 6}, {9, 4}, {9, 6}, {9, 7} };
    /*vector<vector<int>> edge = { {0, 1, 4, 5},
                                {1, 0, 2, 6},
                                {2, 1, 3, 7},
                                {3, 2, 4, 8},
                                {4, 0, 3, 9},
                                {5, 0, 7, 8},
                                {6, 1, 8, 9},
                                {7, 2, 5, 9},
                                {8, 3, 5, 6},
                                {9, 4, 6, 7}, };*/
    vector<vector<int>> edge = { {0, 1, 2, 3, 4},
                                {1, 0, 3},
                                {2, 0, 3},
                                {3, 0, 1, 2},
                                {4, 0, 3} };
    for (int i = 0; i < binary.size(); i += GENE_LEN) {
        int index = i / GENE_LEN;
        string curGene = binary.substr(i, GENE_LEN);
        int var = stoi(curGene, nullptr, 2);
        //将变量值映射到对应的range区间
        var = bound[index].floor + var / (pow(2, GENE_LEN) - 1) * (bound[index].ceil - bound[index].floor);
        x[index] = static_cast<double>(var);
        colors.emplace(x[index]);
        //cout << "x" << index << " = " << x[index] << endl;
    }
    //cout << "----------------------\n";
    indv.setGeneDec(x);
    //相邻顶点颜色相同的顶点个数
    int sameCnt = 0;
    //使用的颜色种数
    int colorNum = colors.size();
    for (auto& e : edge) {
        set<double> temp;
        for (auto& v : e) {
            //使用集合来验证重复元素
            if (temp.find(x[v]) != temp.end()) {
                sameCnt++;
            }
            else {
                temp.emplace(x[v]);
            }
        }
    }
    double fitness = -(sameCnt + colorNum);
    indv.setSameCnt(sameCnt);
    indv.setColorNum(colorNum);
    return fitness;
}

/*
* Gragh Coloring Problem Fitness:
* queen5_5: 25个顶点，320条边
*/
double gcpQueen5_5Func(Individual& indv) {
    string binary = indv.getBinary();
    vector<double> x(GENE_NUM, 0);
    vector<range> bound = { GENE_NUM, range(0, 20) };
    set<double> colors;
    //vector<pair<int, int>> edge = { {0, 1}, {0, 4}, {0, 5}, {1, 2},{1, 0}, {1, 6}, {2, 1}, {2, 3}, {2, 7}, {3, 2}, {3, 4}, {3, 8}, {4, 3}, {4, 8}, {4, 9}, {5, 0}, {5, 7}, {5, 8}, {6, 1}, {6, 8}, {6, 9}, {7, 2}, {7, 5}, {7, 9}, {8, 3}, {8, 5}, {8, 6}, {9, 4}, {9, 6}, {9, 7} };
    vector<vector<int>> edge = { {0, 6, 12, 18, 24, 1, 2, 3, 4, 5, 10, 15, 20},
                                {1, 7, 13, 19, 5, 2, 3, 4, 6, 11, 16, 21, 0},
                                {2, 8, 14, 6, 10, 3, 4, 7, 12, 17, 22, 1, 0},
                                {3, 9, 7, 11, 15, 4, 8, 13, 18, 23, 2, 1, 0},
                                {4, 8, 12, 16, 20, 9, 14, 19, 24, 3, 2, 1, 0},
                                {5, 11, 17, 23, 6, 7, 8, 9, 10, 15, 20, 1, 0},
                                {6, 12, 18, 24, 10, 7, 8, 9, 11, 16, 21, 5, 2, 1, 0},
                                {7, 13, 19, 11, 15, 8, 9, 12, 17, 22, 6, 5, 3, 2, 1},
                                {8, 14, 12, 16, 20, 9, 13, 18, 23, 7, 6, 5, 4, 3, 2},
                                {9, 13, 17, 21, 14, 19, 24, 8, 7, 6, 5, 4, 3},
                                {10, 16, 22, 11, 12, 13, 14, 15, 20, 6, 5, 2, 0},
                                {11, 17, 23, 15, 12, 13, 14, 16, 21, 10, 7, 6, 5, 3, 1},
                                {12, 18, 24, 16, 20, 13, 14, 17, 22, 11, 10, 8, 7, 6, 4, 2, 0},
                                {13, 19, 17, 21, 14, 18, 23, 12, 11, 10, 9, 8, 7, 3, 1},
                                {14, 18, 22, 19, 24, 13, 12, 11, 10, 9, 8, 4, 2},
                                {15, 21, 16, 17, 18, 19, 20, 11, 10, 7, 5, 3, 0},
                                {16, 22, 20, 17, 18, 19, 21, 15, 12, 11, 10, 8, 6, 4, 1},
                                {17, 23, 21, 18, 19, 22, 16, 15, 13, 12, 11, 9, 7, 5, 2},
                                {18, 24, 22, 19, 23, 17, 16, 15, 14, 13, 12, 8, 6, 3, 0},
                                {19, 23, 24, 18, 17, 16, 15, 14, 13, 9, 7, 4, 1},
                                {20, 21, 22, 23, 24, 16, 15, 12, 10, 8, 5, 4, 0},
                                {21, 22, 23, 24, 20, 17, 16, 15, 13, 11, 9, 6, 1},
                                {22, 23, 24, 21, 20, 18, 17, 16, 14, 12, 10, 7, 2},
                                {23, 24, 22, 21, 20, 19, 18, 17, 13, 11, 8, 5, 3},
                                {24, 23, 22, 21, 20, 19, 18, 14, 12, 9, 6, 4, 0} };
    for (int i = 0; i < binary.size(); i += GENE_LEN) {
        int index = i / GENE_LEN;
        string curGene = binary.substr(i, GENE_LEN);
        int var = stoi(curGene, nullptr, 2);
        //将变量值映射到对应的range区间
        var = bound[index].floor + var / (pow(2, GENE_LEN) - 1) * (bound[index].ceil - bound[index].floor);
        x[index] = static_cast<double>(var);
        colors.emplace(x[index]);
        //cout << "x" << index << " = " << x[index] << endl;
    }
    indv.setGeneDec(x);
    //相邻顶点颜色相同的顶点个数
    int sameCnt = 0;
    //使用的颜色种数
    int colorNum = colors.size();
    for (auto& e : edge) {
        set<double> temp;
        for (auto& v : e) {
            //使用集合来验证重复元素
            if (temp.find(x[v]) != temp.end()) {
                sameCnt++;
            }
            else {
                temp.emplace(x[v]);
            }
        }
    }
    //double fitness = 1.0 / (sameCnt * 50000 + colorNum);
    double fitness = 1.0/sameCnt;
    indv.setSameCnt(sameCnt);
    indv.setColorNum(colorNum);
    return fitness;
}

/*
* 适应度计算，计算种群中所有个体的适应度值
* @param population: 种群
*/
void calFitness(vector<Individual>& population, Individual& best, double& f_max, double& f_min) {
    double bestFit = -DBL_MAX;
    double maxFit = -DBL_MAX;
    double minFit = DBL_MAX;
    if (best.getFitness() != 0) {
        bestFit = best.getFitness();
    }
    int bestIdx = 0; //最优个体位置
    int maxIdx = 0;  //f_max位置
    int minIdx = 0;  //f_min位置
    for (int i = 0; i < population.size(); i++) {
        double fitness = gcpFunc(population[i]);
        //记录最优个体
        bestIdx = (fitness > bestFit) ? i : bestIdx;
        bestFit = max(bestFit, fitness);
        population[i].setFitness(fitness);
        //记录本代中最优个体
        maxIdx = (fitness > maxFit) ? i : maxIdx;
        maxFit = max(maxFit, fitness);
        //记录本代中最差个体
        minIdx = (fitness < minFit) ? i : minIdx;
        minFit = min(minFit, fitness);
    }
    //如果本代种群中有比最优个体适应度更高的个体，更新best
    //cout << "best indv = " << best.getFitness() << endl;
    //cout << "bestFit = " << bestFit << endl;
    if (bestFit != best.getFitness()) {
        best = population[bestIdx];
    }
    //根据适应度，从大到小排序种群，并根据移民比率修改flag
    sort(population.begin(), population.end(), [](Individual& a, Individual& b) {return a.getFitness() > b.getFitness(); });
    //for (int i = 0; i <= POP_SIZE * 10 / 100; i++) {
    //    population[i].setSpecFlag(1);
          //cout << population[i].getFitness() << endl;
    //}
    
    //记录最优个体
    f_max = maxFit;
    //population[maxIdx].setSpecFlag(1);
    //记录最差个体
    f_min = minFit;
    //population[minIdx].setSpecFlag(2);
}

/*
* 变异操作，量子非门
*/
void mutation(vector<Individual>& population) {
    for (int i = POP_SIZE - 1; i >= POP_SIZE - POP_SIZE / 10; i--) {
        population[i] = Individual();
        /*for (int j = 0; j < CHROM_LEN; j++) {
            swap(population[i].getChrom()[j].alpha, population[i].getChrom()[j].beta);
            cout << i << " a = " << population[i].getChrom()[j].alpha << endl;
            cout << i << " b = " << population[i].getChrom()[j].beta << endl;
        }*/
        
    }
}

void printPopulation(vector<Individual>& pop) {
    for (int i = 0; i < pop.size(); i++) {
        cout << pop[i].toString() << endl;
    }
}

/*
* QGA主函数
*/
void quantumAlgorithm() {
    Individual best;
    vector<Individual> population;

    //初始化种群
    initPop(POP_SIZE, CHROM_LEN, population, 0);
    //对种群进行一次测量，得到二进制编码
    collapse(population);
    //计算适应度，找出最优个体后再迭代进化
    double f_max = 0;
    double f_min = 0;
    calFitness(population, best, f_max, f_min);
    //进化迭代
    for (int gen = 0; gen < MAX_GEN; gen++) {
        cout << "当前进化代数： " << gen << endl;
        //测量种群
        collapse(population);
        //计算适应度
        double f_max = 0;
        double f_min = 0;
        calFitness(population, best, f_max, f_min);
        //量子旋转门
        qGateAdaptive(population, best, f_max, f_min);
        //qGateRAS_1(population, best);
        //qGateRAS_2(population, best);
        mutation(population);

        //printPopulation(population);
        cout << "##############################\n";
        cout << "best chrom:\n" << best.toString() << endl;
        cout << "##############################\n";
    }
    cout << "--------------------------\n";
    printPopulation(population);
}

int main()
{
    quantumAlgorithm();

    std::cout << "Hello World!\n";
}
