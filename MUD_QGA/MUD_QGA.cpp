#include "MUD_QGA.h"
using namespace std;

/*
* 生成[0, 1]之间的随机数
*/
double srand() {
    int N = rand() % 999;
    double random = static_cast<double>(N) / 1000.0;;//随机产生0到1的小数
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
* @param population: 种群
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
vector<double> decodeBinary(string binary, vector<range>& bound) {
    vector<double> x(GENE_NUM);
    for (int i = 0; i < binary.size(); i += GENE_LEN) {
        int index = i / GENE_LEN;
        string curGene = binary.substr(i, GENE_LEN);
        double var = stoi(curGene, nullptr, 2);
        //将变量值映射到对应的range区间
        var = bound[index].floor + var / (pow(2, GENE_LEN) - 1) * (bound[index].ceil - bound[index].floor);
        x[index] = static_cast<double>(var);
    }
    return x;
}

/*
* 信道容量计算
* @param indv: 个体解
* @param weight: 物品重量
* @param value: 物品价值
*/
double mudFunc(Individual& indv, vector<Individual>& population, vector<double>& gains) {
    string binary = indv.getBinary();
    vector<double> power(GENE_NUM, 0);

    double fitness = 1.0;

    //将二进制解码为十进制
    vector<range> bounds = {GENE_NUM, BOUND};
    //计算每个节点的功率范围
    power = decodeBinary(binary, bounds);
    indv.setGeneDec(power);

    //计算g * p
    vector<double> gp(power.size());
    for (int i = 0; i < gp.size(); i++) {
        gp[i] = gains[i] * power[i];
    }

    for (int i = 0; i < power.size(); i++) {
        //计算已经被解码的节点的干扰
        double innerNoise = 0.0;
        for (int j = 0; j < power.size(); j++) {
            if (j != i) {
                innerNoise += gp[j];
            }
        }
        //计算每一个节点的网络吞吐量
        fitness *= (1 + gp[i] / (0.0001 + innerNoise));
    }
    fitness = log2(fitness);
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
        double fitness = mudFunc(population[i], population, GAINS);
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
    if (bestFit != best.getFitness()) {
        best = population[bestIdx];
    }
    //根据适应度，从大到小排序种群，并根据移民比率修改flag
    sort(population.begin(), population.end(), [](Individual& a, Individual& b) {return a.getFitness() > b.getFitness(); });

    //记录最优个体
    f_max = maxFit;
    //记录最差个体
    f_min = minFit;
}


/*
* 量子灾变函数
* 初始化种群中部分或全部个体
* @param population: 种群
* @param initLot: 灾变占种群的份额，1为全体灾变
*/
void catastrophe(vector<Individual>& population, int initLot) {
    for (int i = POP_SIZE - 1; i >= POP_SIZE - POP_SIZE * initLot; i--) {
        population[i] = Individual();
    }
}

/*
* 量子变异函数
* 根据种群中所有个体与最优解的海明距离控制变异概率
* @param population: 种群
* @param best: 最优个体
*/
void mutation(vector<Individual>& population, Individual& best) {
    //变异概率从0.001到0.01变化
    double mutaRate = 0.001;
    //评价当前种群的与最优个体间的平均海明距离，最大距离为个体编码长度
    int hammingSum = 0;
    for (auto& x : population) {
        hammingSum += x.getHammingDis(best);
    }
    double avgHamming = static_cast<double>(hammingSum) / POP_SIZE;
    cout << "haming = " << avgHamming << '\n';
    //当平均海明距离小于40% * 最大海明距离时，变异概率为最低
    if (avgHamming < 0.4 * CHROM_LEN) {
        mutaRate = 0.005;
    }
    if (avgHamming < 0.2 * CHROM_LEN) {
        mutaRate = 0.01;
    }
    //当平均海明距离小于10%时，认为算法收敛，进行量子灾变
    if (avgHamming < 0.1 * CHROM_LEN) {
        mutaRate = 0.1;
        catastrophe(population, 1);
        return;
    }

    //根据mutaRate进行随机变异
    for (int i = 0; i < POP_SIZE; i++) {
        double pick = srand();
        if (pick <= mutaRate) {
            for (int j = 0; j < CHROM_LEN; j++) {
                swap(population[i].getChrom()[j].alpha, population[i].getChrom()[j].beta);
            }
        }
    }
}


/*
* 禁忌搜索算子
* @param population: 种群
* @param best: 禁忌搜索的初始解
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
    deque<tabuItem> tabuList(TABU_LIST_SIZE);
    //将best作为初始解，迭代局部最优
    for (int m = 0; m < gens; m++) {
        cout << "Tabu search 代数：" << m << '\n';
        string bestBin = best.getBinary();
        double bestFit = best.getFitness();
        vector<tabuItem> tabuTemp;  //这里用于存储best所有对换位置结果

        //两两交换当前best的两个顶点
        for (int i = 0; i < CHROM_LEN; i++) {
            if (bestBin[i] == '0') {
                for (int j = 0; j < CHROM_LEN; j++) {
                    if (bestBin[j] == '1') {
                        Individual tempChrom = Individual();
                        swap(bestBin[i], bestBin[j]);
                        tempChrom.setBinary(bestBin);
                        //判定对换后的适应度
                        double fitness = mudFunc(population[i], population, GAINS);
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
        //将所有交换后适应度改变的结果从大到小排序，维护禁忌表
        sort(tabuTemp.begin(), tabuTemp.end(), [](tabuItem a, tabuItem b) {return a.fitDiff > b.fitDiff; });

        //选出不在tabuList中且提升最大的对换策略
        for (auto& x : tabuTemp) {
            if (find(tabuList.begin(), tabuList.end(), x) != tabuList.end()) {
                //突破禁忌
                if (x.fitDiff + bestFit > elite.getFitness()) {
                    best = x.newIndv;       //更新禁忌最优值
                    elite = best;           //更新全局最优值
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
        //最后一次迭代，将迭代中最优解elite赋值给best
        if (m == gens - 1) {
            best = elite;
        }
        cout << "---------------------\n";
        cout << "new best = " << best.toString() << endl;
        cout << "---------------------\n";
    }

    //如果迭代后最优解没有改变，则停留Flag加一，调整下次禁忌搜索次数
    if (preElite.getFitness() == best.getFitness()) {
        stayFlag++;
    }
    else {
        stayFlag = 0;
    }
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
        if (!flag) {
            flag = 50;
            tabu(population, best);
        }
        mutation(population, best);
    }
    cout << "--------------------------\n";
    //printPopulation(population);
}

int main()
{
    /*
    * β   = 0.5
    * σ²  = 0.05
    * pmax = 0.1
    */
    double beta = 5;
    double sigma = 0.0001;
    double pmax = 0.1;
    //节点坐标
    vector<pair<double, double>> coords;

    //计算每个节点的功率取值范围
    double Aj = 1 + log(pmax / beta*sigma) / log(beta + 1);
    double floor = beta * sigma * pow((beta + 1), (Aj + 1));
    double ceil = pmax;
    BOUND = range(floor, ceil);
    cout << BOUND.toString() << endl;
    

    //节点到终端的距离
    DISTANCE = { 81.4723686393179,
                90.5791937075619,
                12.6986816293506,
                91.3375856139019,
                63.2359246225410,
                9.75404049994095,
                27.8498218867048,
                54.6881519204984,
                95.7506835434298,
                96.4888535199277 };
    //每个节点的信道增益
    GAINS = {
        5.54740950894548e-06,
        4.03678763845922e-06,
        0.00146502605831966,
        3.93706579298868e-06,
        1.18639589518653e-05,
        0.00323271628650082,
        0.000138884566867956,
        1.83417813620838e-05,
        3.41739835278010e-06,
        3.33956434986324e-06
    };
    srand(unsigned(time(NULL)));
    quantumAlgorithm();

    

    std::cout << "Hello World!\n";
}
