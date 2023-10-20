#include "PQGA.h"
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
* 多宇宙机制初始化，采用四宇宙模型
* u1为主宇宙
*/
void initUniverse(vector<Individual>& u1, vector<Individual>& u2, vector<Individual>& u3, vector<Individual>& u4) {
    initPop(POP_SIZE, CHROM_LEN, u1, 0);
    initPop(POP_SIZE, CHROM_LEN, u2, 1);
    initPop(POP_SIZE, CHROM_LEN, u3, 1);
    initPop(POP_SIZE, CHROM_LEN, u4, 1);
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
void qGateAdaptive(vector<Individual>& population, Individual& best, double f_max, double f_min) {
    const double CONST_ARG = (f_max > f_min) ? (K2 - K1) / (f_max - f_min) : 0;

    double bestFit = best.getFitness();                         //最优个体适应度
    for (int i = 0; i < POP_SIZE; i++) {
        //根据当前个体计算旋转角
        double theta = K1 + (population[i].getFitness() - f_min) * CONST_ARG;
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
* 0-1背包问题
* @param indv: 个体解
* @param weight: 物品重量
* @param value: 物品价值
*/
double packFunc(Individual& indv, vector<double> weight, vector<double> value) {
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

    for (int i = 0; i < value.size(); i++) {
        fitness += x[i] * value[i];
    }

    if (weightTemp > PACK_MAX) {
        fitness -= 100 * (weightTemp - PACK_MAX);
    }

    return fitness;
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
* 适应度计算，计算种群中所有个体的适应度值
* @param population: 种群
*/
void calFitness(vector<Individual>& population, Individual& best, double f_max, double f_min) {
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
        //double fitness = objFunc(population[i]);
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
    if (bestFit != best.getFitness()) {
        best = population[bestIdx];
    }
    //根据适应度，从大到小排序种群
    sort(population.begin(), population.end(), [](Individual& a, Individual& b) {return a.getFitness() > b.getFitness(); });
    
    //记录最优个体
    f_max = maxFit;
    //记录最差个体
    f_min = minFit;
}

/*
* 禁忌搜索算子，用于局部搜索
*/
void tabu(vector<Individual>& population, Individual& best) {
    //记录进入禁忌搜索前的最优值
    static Individual preElite;
    static int stayFlag = 0;
    if (preElite.getFitness() != best.getFitness()) {
        stayFlag = 0;
    }
    preElite = best;

    //禁忌搜索实际迭代次数计算
    int gens = TABU_GEN;
    if (stayFlag == 0) {
        gens = TABU_GEN / 4;
    }
    else if (stayFlag == 1) {
        gens = TABU_GEN / 2;
    }
    else if (stayFlag > 5) {
        stayFlag = 0;
    }

    //建立最优记录，保证最后输出的是迭代过程中的最优值
    Individual elite = best;
    //建立禁忌表
    deque<tabuItem> tabuList(35);

    //确定扔掉、添加的物品数量，范围是(-4, 4)
    int flag = -4 + 8 * srand();

    int temp = flag;

    //将best作为初始解，迭代局部最优
    for (int m = 0; m < gens; m++) {
        cout << "Tabu search 代数：" << m << '\n';
        string bestBin = best.getBinary();
        double bestFit = best.getFitness();
        vector<tabuItem> tabuTemp;  //这里用于存储best所有对换位置结果

        cout << "stayFlag = " << stayFlag << '\n';

        cout << "Flag = " << temp << '\n';
        //上次未找到更优解，添加一个物品
        if (stayFlag != 0) {
            for (auto& x : bestBin) {
                //决定是多拿一个物品还是丢掉一个物品
                if (flag > 0) {
                    //添加物品
                    if (x == '0') {
                        x = '1';
                        flag--;
                        if (flag <= 0) {
                            break;
                        }
                    }
                }
                else if (flag < 0) {
                    //丢掉物品
                    if (x == '1') {
                        x = '0';
                        flag++;
                        if (flag >= 0) {
                            break;
                        }
                    }
                }
                else {
                    break;
                }
            }


            best.setBinary(bestBin);
            double fitness = packFunc(best, PACK_WEIGHT, PACK_VAL);
            best.setFitness(fitness);
            bestFit = fitness;
        }
        //两两交换当前best的两个顶点
        for (int i = 0; i < CHROM_LEN; i++) {
            if (bestBin[i] == '0') {
                for (int j = 0; j < CHROM_LEN; j++) {
                    if (bestBin[j] == '1') {
                        Individual tempChrom = Individual();
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
                }
            }
        }
        //将所有交换后适应度改变的结果从大到小排序
        sort(tabuTemp.begin(), tabuTemp.end(), [](tabuItem a, tabuItem b) {return a.fitDiff > b.fitDiff; });

        //选出不在tabuList中且提升最大的对换策略
        for (auto& x : tabuTemp) {
            if (find(tabuList.begin(), tabuList.end(), x) != tabuList.end()) {
                //突破禁忌
                if (x.fitDiff + bestFit > elite.getFitness()) {
                    best = x.newIndv;
                    elite = best;           //更新禁忌最优值
                    //将新的选择加入禁忌表
                    tabuList.pop_front();
                    tabuList.push_back(x);
                    break;
                }
                continue;
            }
            else {
                best = x.newIndv;
                if (best.getFitness() >= elite.getFitness()) {
                    elite = best;
                }
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
        }

        int flag = 10;
        for (auto& x : tabuTemp) {
            if (flag <= 0) break;
            flag--;
            cout << x.toString();
        }*/
    }
    //将迭代中最优解赋值给best
    best = elite;
    //如果迭代后最优解没有改变，则下次进行禁忌搜索时加入物品
    if (preElite.getFitness() == best.getFitness()) {
        stayFlag++;
    }
    else {
        stayFlag = 0;
    }
}

/*
* 移民算子，将从宇宙的最优解移民到u1，并替换掉u1较差的解，同时主宇宙将一些次优解回馈给从宇宙
* @param u1 主宇宙
* @param u2 从宇宙
* @param u3 从宇宙
* @param u4 从宇宙
*/
void migration(vector<Individual>& u1, vector<Individual>& u2, vector<Individual>& u3, vector<Individual>& u4) {
    cout << "Migrating... " << '\n';
    //各个从宇宙的最优解移民道到主宇宙
    int m = POP_SIZE * MIGRATE_RATE - 1;
    int n = POP_SIZE * MIGRATE_RATE - 1;
    int k = POP_SIZE * MIGRATE_RATE - 1;
    for (int i = POP_SIZE - 1; i >= POP_SIZE - POP_SIZE * MIGRATE_RATE * 3; i--) {
        if (i >= POP_SIZE - POP_SIZE * MIGRATE_RATE) {
            u1[i] = u2[m--];
        }
        else if (i < POP_SIZE - POP_SIZE * MIGRATE_RATE * 2) {
            u1[i] = u3[n--];
        }
        else {
            u1[i] = u4[k--];
        }
    }
    //cout << m << " " << n << " " << k << endl;
}

void printPopulation(vector<Individual>& pop) {
    for (int i = 0; i < pop.size(); i++) {
        cout << pop[i].toString() << endl;
    }
}

/*
* 多宇宙并行自适应QGA主函数
*/
void maQuantumAlgorithm() {
    Individual best;        //全局最优解
    vector<Individual> u1;
    vector<Individual> u2;
    vector<Individual> u3;
    vector<Individual> u4;

    //建立线程池
    vector<thread> uThreads(4);

    //初始化各个宇宙
    initUniverse(u1, u2, u3, u4);
    //对种群进行一次测量，得到二进制编码
    collapse(u1);
    collapse(u2);
    collapse(u3);
    collapse(u4);
    /*uThreads[0] = thread(collapse, ref(u2));
    uThreads[1] = thread(collapse, ref(u3));
    uThreads[2] = thread(collapse, ref(u4));
    uThreads[3] = thread(collapse, ref(u1));
    uThreads[0].join();
    uThreads[1].join();
    uThreads[2].join();
    uThreads[3].join();*/
 
    //计算适应度
    double f_max1 = 0, f_min1 = 0;
    double f_max2 = 0, f_min2 = 0;
    double f_max3 = 0, f_min3 = 0;
    double f_max4 = 0, f_min4 = 0;
    calFitness(u1, best, f_max1, f_min1);
    calFitness(u2, best, f_max2, f_min2);
    calFitness(u3, best, f_max3, f_min3);
    calFitness(u4, best, f_max4, f_min4);
    /*uThreads[0] = thread(calFitness, ref(u2), ref(best), f_max2, f_min2);
    uThreads[1] = thread(calFitness, ref(u3), ref(best), f_max3, f_min3);
    uThreads[2] = thread(calFitness, ref(u4), ref(best), f_max4, f_min4);
    uThreads[0].join();
    uThreads[1].join();
    uThreads[2].join();*/

    int flag = 2;
    while (flag--) {
        //各个宇宙进化迭代
        for (int gen = 1; gen <= MAX_GEN; gen++) {
            cout << "当前进化代数： " << gen << endl;
            //测量种群
             
            /*uThreads[0] = thread(collapse, ref(u2));
            uThreads[1] = thread(collapse, ref(u3));
            uThreads[2] = thread(collapse, ref(u4));
            uThreads[3] = thread(collapse, ref(u1));
            uThreads[0].join();
            uThreads[1].join();
            uThreads[2].join();
            uThreads[3].join();*/

            //计算适应度
            double f_max1 = 0, f_min1 = 0;
            double f_max2 = 0, f_min2 = 0;
            double f_max3 = 0, f_min3 = 0;
            double f_max4 = 0, f_min4 = 0;
            calFitness(u1, best, f_max1, f_min1);
            calFitness(u2, best, f_max2, f_min2);
            calFitness(u3, best, f_max3, f_min3);
            calFitness(u4, best, f_max4, f_min4);
            /*uThreads[0] = thread(calFitness, ref(u2), ref(best), f_max2, f_min2);
            uThreads[1] = thread(calFitness, ref(u3), ref(best), f_max3, f_min3);
            uThreads[2] = thread(calFitness, ref(u4), ref(best), f_max4, f_min4);
            uThreads[0].join();
            uThreads[1].join();
            uThreads[2].join();*/
            //量子旋转门
            qGateAdaptive(u2, best, f_max2, f_min2);
            qGateAdaptive(u3, best, f_max3, f_min3);
            qGateAdaptive(u4, best, f_max4, f_min4);
            qGateAdaptive(u1, best, f_max1, f_min1);
            /*uThreads[0] = thread(qGateAdaptive, ref(u2), ref(best), f_max2, f_min2);
            uThreads[1] = thread(qGateAdaptive, ref(u3), ref(best), f_max3, f_min3);
            uThreads[2] = thread(qGateAdaptive, ref(u4), ref(best), f_max4, f_min4);
            uThreads[0].join();
            uThreads[1].join();
            uThreads[2].join();*/
            cout << "best chrom:\n" << best.toString() << endl;
            //移民算子
            
        }
        
        migration(u1, u2, u3, u4);
    }
}

int main()
{
    int testNum = 9;
    switch (testNum)
    {
    case 1:
    {
        /* 1 (10,269) */
        PACK_WEIGHT = { 95, 4, 60, 32, 23, 72, 80, 62, 65, 46 };
        PACK_VAL = { 55, 10, 47, 5, 4, 50, 8, 61, 85, 87 };
    }
    break;
    case 2:
    {
        /* 2 (15,375) */
        PACK_WEIGHT = { 56.358531, 80.874050, 47.987304, 89.596240, 74.66048, 85.894345, 51.353496, 1.498459, 36.445204, 16.589862, 44.56923, 0.4669, 37.788018, 57.118442, 60.716575 };
        PACK_VAL = { 0.125126, 19.330424, 58.500931, 35.029145, 82.284005, 17.410810, 71.050142, 30.399487, 9.140294, 14.731285, 98.852504, 11.908322, 0.891140, 53.166295, 60.176397 };
    }
    break;
    case 3:
    {
        /* 3 (20,878) */
        PACK_WEIGHT = { 92, 4, 43, 83, 84, 68, 92, 82, 6, 44, 32, 18, 56, 83, 25, 96, 70, 48, 14, 58 };
        PACK_VAL = { 44, 46, 90, 72, 91, 40, 75, 35, 8, 54, 78, 40, 77, 15, 61, 17, 75, 29, 75, 63 };
    }
    break;
    case 4:
    {
        /* 4 (23,10000) */
        PACK_WEIGHT = { 983, 982, 981, 980, 979, 978, 488, 976, 972, 486, 486, 972, 972, 485, 485, 969, 966, 483, 964, 963, 961, 958, 959 };
        PACK_VAL = { 981, 980, 979, 978, 977, 976, 487, 974, 970, 485, 485, 970, 970, 484, 484, 976, 974, 482, 962, 961, 959, 958, 857 };
    }
    break;
    case 5:
    {
        /* 5 (50,1000) */
        PACK_WEIGHT = { 80, 82, 85, 70, 72, 70, 66, 50, 55, 25, 50, 55, 40, 48, 50, 32, 22, 60, 30, 32, 40, 38, 35, 32, 25, 28, 30, 22, 50, 30, 45, 30, 60, 50, 20, 65, 20, 25,
30, 10, 20, 25, 15, 10, 10, 10, 4, 4, 2, 1 };
        PACK_VAL = { 220, 208, 198, 192, 180, 180, 165, 162, 160, 158, 155, 130, 125, 122, 120, 118, 115, 110, 105, 101, 100, 100, 98, 96,
95, 90, 88, 82, 80, 77, 75, 73, 72, 70, 69, 66, 65, 63, 60, 58, 56, 50, 30, 20, 15, 10, 8, 5, 3, 1 };
    }
    break;
    case 6:
    {
        /* 6 (50,11258) */
        PACK_WEIGHT = { 438, 754, 699, 587, 789, 912, 819, 347, 511, 287, 541, 784, 676, 198, 572, 914, 988, 4, 355, 569, 144, 272, 531,
556, 741, 489, 321, 84, 194, 483, 205, 607, 399, 747, 118, 651, 806, 9, 607, 121, 370, 999, 494, 743, 967, 718, 397, 589,
193, 369 };
        PACK_VAL = { 72, 490, 651, 833, 883, 489, 359, 337, 267, 441, 70, 934, 467, 661, 220, 329, 440, 774, 595, 98, 424, 37, 807, 320,
501, 309, 834, 851, 34, 459, 111, 253, 159, 858, 793, 145, 651, 856, 400, 285, 405, 95, 391, 19, 96, 273, 152, 473, 448,
231 };
    }
    break;
    case 7:
    {
        /* 7 (60,2400) */
        PACK_WEIGHT = { 135, 133, 130, 11, 128, 123, 20, 75, 9, 66, 105, 43, 18, 5, 37, 90, 22, 85, 9, 80, 70, 17, 60, 35, 57, 35, 61, 40, 8, 50, 32,
40, 72, 35, 100, 2, 7, 19, 28, 10, 22, 27, 30, 88, 91, 47, 68, 108, 10, 12, 43, 11, 20, 37, 17, 4, 3, 21, 10, 67 };
        PACK_VAL = { 350, 310, 300, 295, 290, 287, 283, 280, 272, 270, 265, 251, 230, 220, 215, 212, 207, 203, 202, 200, 198, 196, 190, 182, 181, 175,
160, 155, 154, 140, 132, 125, 110, 105, 101, 92, 83, 77, 75, 73, 72, 70, 69, 66, 60, 58, 45, 40, 38, 36, 33, 31, 27, 23, 20,
19, 10, 9, 4, 1 };
    }
    break;
    case 8:
    {
        /* 8 (80,1173) */
        PACK_WEIGHT = { 40, 27, 5, 21, 51, 16, 42, 18, 52, 28, 57, 34, 44, 43, 52, 55, 53, 42, 47, 56, 57, 44, 16, 2, 12, 9, 40, 23, 56, 3, 39, 16, 54,
36, 52, 5, 53, 48, 23, 47, 41, 49, 22, 42, 10, 16, 53, 58, 40, 1, 43, 56, 40, 32, 44, 35, 37, 45, 52, 56, 40, 2, 23, 49, 50, 26,
11, 35, 32, 34, 58, 6, 52, 26, 31, 23, 4, 52, 53, 19 };
        PACK_VAL = { 199, 194, 193, 191, 189, 178, 174, 169, 164, 164, 161, 158, 157,
154, 152, 152, 149, 142, 131, 125, 124, 124, 124, 122, 119, 116, 114, 113, 111, 110, 109, 100, 97, 94, 91, 82, 82, 81, 80,
80, 80, 79, 77, 76, 74, 72, 71, 70, 69, 68, 65, 65, 61, 56, 55, 54, 53, 47, 47, 46, 41, 36, 34, 32, 32, 30, 29, 29, 26, 25, 23,
22, 20, 11, 10, 9, 5, 4, 3, 1 };
    }
    break;
    case 9:
    {
        /* 9 (100,3820) */
        PACK_WEIGHT = { 54, 95, 36, 18, 4, 71, 83, 16, 27, 84, 88, 45, 94, 64, 14, 80, 4, 23, 75, 36, 90, 20, 77, 32, 58, 6, 14, 86, 84, 59, 71, 21,
30, 22, 96, 49, 81, 48, 37, 28, 6, 84, 19, 55, 88, 38, 51, 52, 79, 55, 70, 53, 64, 99, 61, 86, 1, 64, 32, 60, 42, 45, 34, 22, 49,
37, 33, 1, 78, 43, 85, 24, 96, 32, 99, 57, 23, 8, 10, 74, 59, 89, 95, 40, 46, 65, 6, 89, 84, 83, 6, 19, 45, 59, 26, 13, 8, 26, 5,
9 };
        PACK_VAL = { 297, 295, 293, 292, 291, 289, 284, 284, 283, 283, 281, 280, 279, 277, 276, 275, 273, 264, 260, 257, 250, 236, 236,
    235, 235, 233, 232, 232, 228, 218, 217, 214, 211, 208, 205, 204, 203, 201, 196, 194, 193, 193, 192, 191, 190, 187, 187,
    184, 184, 184, 181, 179, 176, 173, 172, 171, 160, 128, 123, 114, 113, 107, 105, 101, 100, 100, 99, 98, 97, 94, 94, 93, 91,
    80, 74, 73, 72, 63, 63, 62, 61, 60, 56, 53, 52, 50, 48, 46, 40, 40, 35, 28, 22, 22, 18, 15, 12, 11, 6, 5 };
    }
    break;
    case 10:
    {
        /* 10 (100,282) */
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
    }
    break;
    default:
        break;
    }

    clock_t startTime = clock();
    maQuantumAlgorithm();
    clock_t endTime = clock();
    double costTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    cout << "costTime = " << costTime * 1000 << "ms" << endl;

    std::cout << "Hello World!\n";
}
