#include "PACK_QGA.h"
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
    case 2: //随机初始化，针对TSP问题
        population.resize(M);
        for (int i = 0; i < M; i++) {
            vector<qubit> chrom(CHROM_LEN);
            for (int j = 0; j < CHROM_LEN; j++) {
                double seed = srand();
                qubit initQ = qubit(sqrt(seed), sqrt(1 - seed));
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
        //population[i].setBinary("00010010011010010000");
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
    if (f_max == f_min) return;
    const double CONST_ARG = (f_max > f_min) ? (K2 - K1) / (f_max - f_min) : 0;

    double bestFit = best.getFitness();                         //最优个体适应度
    for (int i = 0; i < POP_SIZE; i++) {
        //根据当前个体计算旋转角
        double theta = K1 + (population[i].getFitness() - f_min) * CONST_ARG;
        //个体变异概率，默认为80%
        double mutationPick = srand();
        if (mutationPick > 1) {
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
* 0-1背包问题
* @param indv: 个体解
* @param weight: 物品重量
* @param value: 物品价值
*/
double packFunc(Individual& indv, vector<int> weight, vector<int> value) {
    string binary = indv.getBinary();
    vector<double> x(GENE_NUM, 0);

    for (int i = 0; i < binary.size(); i += GENE_LEN) {
        int index = i / GENE_LEN;
        string curGene = binary.substr(i, GENE_LEN);
        int var = stoi(curGene, nullptr, 2);

        //1代表选，0代表不选
        x[index] = static_cast<double>(var);
        //cout << "x" << index << " = " << x[index] << endl;
    }
    indv.setGeneDec(x);

    double fitness = 0.0;
    double weightTemp = 0.0;
    for (int i = 0; i < weight.size(); i++) {
        weightTemp += x[i] * weight[i];
    }
    //cout << "wei = " << weightTemp << endl;

    indv.setColorNum(weightTemp);
    for (int i = 0; i < value.size(); i++) {
        fitness += x[i] * value[i];
    }
    indv.setSameCnt(fitness);
    if (weightTemp > 282) {
        fitness -= 100 * (weightTemp - 282);
    }

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
        double fitness = packFunc(population[i], PACK_WEIGHT, PACK_VAL);
        //记录最优个体
        bestIdx = (fitness >= bestFit) ? i : bestIdx;
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
* TO-DO：下标问题
*/
void mutation(vector<Individual>& population) {
    //for (int i = POP_SIZE - 1; i >= POP_SIZE - POP_SIZE / 2; i--) {
    for (int i = POP_SIZE - 1; i >= 0; i--) {
        population[i] = Individual();
        /*for (int j = 0; j < CHROM_LEN; j++) {
            swap(population[i].getChrom()[j].alpha, population[i].getChrom()[j].beta);
            cout << i << " a = " << population[i].getChrom()[j].alpha << endl;
            cout << i << " b = " << population[i].getChrom()[j].beta << endl;
        }*/

    }
}

/*
* 禁忌搜索算子，用于局部搜索
*/
void tabu(vector<Individual>& population, Individual& best) {
    //评价当前种群的与最优个体间的平均海明距离，最大距离为个体编码长度
    int hammingSum = 0;
    for (auto& x : population) {
        hammingSum += x.getHammingDis(best);
    }
    double avgHamming = static_cast<double>(hammingSum) / POP_SIZE;
    cout << "haming = " << avgHamming << '\n';
    //当平均海明距离小于20% * 最大海明距离时，进行禁忌搜索
    if (avgHamming > 0.2 * CHROM_LEN) {
        return;
    }
    //当平均海明距离小于10%时，认为算法收敛，进行量子灾变
    if (avgHamming <= 0.1 * CHROM_LEN) {
        mutation(population);
        return;
    }
    // return;
    //建立禁忌表
    deque<tabuItem> tabuList(10);

    //将best作为初始解，迭代局部最优
    Individual tempChrom = Individual();
    vector<tabuItem> tabuTemp;  //这里用于存储best所有对换位置结果

    for (int m = 0; m < TABU_GEN; m++) {
        cout << "Tabu search 代数：" << m << '\n';
        string bestBin = best.getBinary();
        double bestFit = best.getFitness();
        //两两交换当前best的两个顶点
        for (int i = 0; i < CHROM_LEN; i ++) {
            if (bestBin[i] == '1') {
                for (int j = 0; j < CHROM_LEN; j++) {
                    if (bestBin[j] == '0') {
                        swap(bestBin[i], bestBin[j]);
                        tempChrom.setBinary(bestBin);
                        //判定对换后的适应度
                        double fitness = packFunc(tempChrom, PACK_WEIGHT, PACK_VAL);
                        tempChrom.setFitness(fitness);
                        fitness = fitness - bestFit;
                        tabuItem item = tabuItem(pair<int, int>(i, j), fitness, tempChrom);
                        tabuTemp.push_back(item);
                        //恢复best
                        bestBin = best.getBinary();
                    }
                    else {
                        continue;
                    }
                }
            }
        }
        //将所有交换后适应度改变的结果从大到小排序
        sort(tabuTemp.begin(), tabuTemp.end(), [](tabuItem a, tabuItem b) {return a.fitDiff > b.fitDiff; });

        //选出不在tabuList中且提升最大的对换策略
        for (auto& x : tabuTemp) {
            if (find(tabuList.begin(), tabuList.end(), x) != tabuList.end()) {
                //突破禁忌
                if (x.fitDiff + bestFit > bestFit) {
                    best = x.newIndv;
                    //将新的选择加入禁忌表
                    tabuList.pop_front();
                    tabuList.push_back(x);
                    break;
                }
                continue;
            }
            else {
                best = x.newIndv;
                //将新的选择加入禁忌表
                tabuList.pop_front();
                tabuList.push_back(x);
                break;
            }
        }
        cout << "---------------------\n";
        cout << "new best = " << best.toString() << endl;
        cout << "---------------------\n";

        /*for (auto& x : tabuList) {
            cout << "tabuList = " << x.toString() << '\n';
        }*/

        /*int flag = 10;
        for (auto& x : tabuTemp) {
            if (flag <= 0) break;
            flag--;
            cout << x.toString();
        }*/
    }

}

/*
* 读取无向图的边数据
* @param edge: 边
* @param path: 图的txt路径
*/
int readEdge(vector<pair<int, int>>& edge, string path) {
    ifstream infile;
    infile.open(path, ios::in);
    if (!infile.is_open())
    {
        cout << "读取文件失败" << endl;
        return 0;
    }
    //读取每一条边，并存入edge
    string line;
    while (getline(infile, line))
    {
        auto it = line.find(' ');
        int first = atoi(line.substr(0, it).c_str()) - 1;
        int second = atoi(line.substr(it + 1).c_str()) - 1;
        edge.push_back(pair<int, int>(first, second));
    }
    return 1;
}

/*
* 打印种群
*/
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
    int flag = 50;
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


        //printPopulation(population);
        cout << "##############################\n";
        cout << "best chrom:\n" << best.toString() << endl;
        cout << "##############################\n";
        flag--;
        if (flag == 0) {
            flag = 50;
            
            tabu(population, best);
            
        }
        
        
        //mutation(population);
    }
    cout << "--------------------------\n";
    //printPopulation(population);
}

int main()
{
    PACK_VAL = { 25, 24, 92, 30, 42, 22, 30, 70, 75, 60, 76, 96, 16, 59, 47, 94, 29, 99,
                    17, 31, 21, 85, 48, 90, 48, 29, 98, 73, 42, 45, 7, 71, 74, 84, 33, 90,
                    80, 87, 71, 42, 7, 20, 65, 96, 67, 42, 90, 36, 48, 34, 63, 28, 49, 97,
                    78, 71, 86, 4, 75, 58, 23, 51, 70, 18, 35, 70, 43, 72, 94, 79, 64, 92,
                    55, 46, 75, 87, 74, 27, 2, 21, 87, 32, 29, 23, 77, 80, 43, 98, 12, 29,
                    25, 90, 16, 83, 53, 56, 20, 51, 57, 25 };

    PACK_WEIGHT = { 34, 21, 56, 30, 20, 8, 39, 46, 24, 77, 30, 47, 18, 10, 15, 56, 72, 12, 89,
                    67, 60, 48, 73, 25, 31, 72, 86, 32, 70, 99, 21, 83, 95, 24, 73, 30, 77,
                    52, 36, 16, 61, 35, 28, 23, 46, 35, 48, 35, 94, 21, 21, 11, 2, 14, 45,
                    2, 14, 83, 35, 56, 66, 18, 92, 91, 2, 39, 71, 86, 79, 35, 7, 53, 60, 40,
                    40, 83, 30, 47, 10, 10, 69, 49, 94, 50, 93, 51, 62, 23, 89, 26, 13,
                    43, 79, 59, 21, 77, 86, 28, 83, 57 };
    /*PACK_WEIGHT = { 4, 6, 3, 8, 6 };
    PACK_VAL = { 2, 5, 8, 3, 5 };*/
    TSP_CITIES = { {16.47, 96.10}, {16.47, 94.44}, {20.09, 92.54}, {22.39, 93.37}, {25.23, 97.24}, {22.00, 96.05}, {20.47, 97.02},
                    {17.20, 96.29}, {16.30, 97.38}, {14.05, 98.12}, {16.53, 97.38}, {21.52, 95.59}, {19.41, 97.13}, {20.09, 94.55} };
    srand(unsigned(time(NULL)));
    quantumAlgorithm();

    std::cout << "Hello World!\n";
}
