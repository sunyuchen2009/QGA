#include "MUD_QGA.h"
using namespace std;

//����[0, 1]֮��������
double srand() {
    int N = rand() % 999;
    double random = static_cast<double>(N) / 1000.0;;//�������0��1��С��
    //cout << random << endl;
    return random;
}

/*
* ��ʼ����Ⱥ��ʹ�����ӱ��ر���
* @param M ��Ⱥ��С
* @param N ÿ����������ӱ��ر��볤��
* @param strategyFlag ��ʼ�����ԣ�0��Ĭ����������̬�ȸ��ʣ�
                                  1��С������ʼ�����ԣ� (��, ��) = (sqrt(j / n), sqrt(1 - j / n))
*/
void initPop(int M, int N, vector<Individual>& population, int strategyFlag) {
    switch (strategyFlag)
    {
    case 0:
        //����Ⱥ�����е�Ⱦɫ���ʼ��Ϊ���Ŷ���֮һ
        population.resize(M);
        for (int i = 0; i < M; i++) {
            population[i] = Individual();
        }
        break;
    case 1: //С����
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
    case 2: //�����ʼ�������TSP����
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
* ������������ɶ���Ⱥ��һ�β���
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
        //ͬʱ��ʼ��specFlag;
        population[i].setSpecFlag(0);
        //population[i].setBinary("00010010011010010000");
    }
}

/*
* ������ת��QGate Adaptive
*   x	best 	f(x)>=f(best)	delta	    ab>0	ab<0	a=0	    b=0
*   0	0		false			0		    0		0		0	    0
*   0	0		true			0		    0		0		0	    0
*   0	1		false			��   	    1		-1		0	    +-1
*   0	1		true			��   	    -1		1		+-1	    0
*   1	0		false			��   	    -1		1		+-1	    0
*   1	0		true			��   	    1		-1		0	    +-1
*   1	1		false			0		    0		0		0	    0
*   1	1		true			0		    0		0		0	    0
* ���У���ת�Ǧ��ļ��㹫ʽ����:
* �� = K1 + (K2 - K1) * (f_i - f_min) / (f_max - f_min)  [f_min != f_max]
* �� = K1                                                [f_min == f_max]
* f_i��f_max��f_min�ֱ������ǰ�������Ӧ�ȡ���ǰ��Ⱥ���Ÿ�����Ӧ���Լ������Ӧ��
* K1��K2Ϊ������������K1 < K2�����ڿ��������ٶ�
*/
void qGateAdaptive(vector<Individual>& population, Individual& best, double& f_max, double& f_min) {
    if (f_max == f_min) return;
    const double CONST_ARG = (f_max > f_min) ? (K2 - K1) / (f_max - f_min) : 0;

    double bestFit = best.getFitness();                         //���Ÿ�����Ӧ��
    for (int i = 0; i < POP_SIZE; i++) {
        //���ݵ�ǰ���������ת��
        double theta = K1 + (population[i].getFitness() - f_min) * CONST_ARG;
        //���������ʣ�Ĭ��Ϊ80%
        double mutationPick = srand();
        if (mutationPick > 1) {
            continue;
        }
        for (int j = 0; j < CHROM_LEN; j++) {
            double alpha = population[i].getChrom()[j].alpha;   //��
            double beta = population[i].getChrom()[j].beta;     //��
            double delta = 0.0;                                 //��ת�Ǵ�С
            char x = population[i].getBinary()[j];              //��ǰ�����j������
            char b = best.getBinary()[j];                       //���Ÿ����j������
            int s = 0;                                          //��ת�Ƿ��򣬼�������
            double curFit = population[i].getFitness();         //��ǰ������Ӧ��
            //��ǰ���������Ÿ�����ͬ������ת
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
            double e = s * delta;       //��ת��
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
                //H-ep��
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
* ��������Ʊ��룬����ʮ���ƣ���ת�Ƶ����������Ƶ�����
* @param binary: ����Ķ����Ʊ��룬�����������
* @param bound: ������������
* @return: �������򣨱�����������
*/
vector<double> decodeBinary(string binary, vector<range> bound) {
    vector<double> x(GENE_NUM);
    for (int i = 0; i < binary.size(); i += GENE_LEN) {
        int index = i / GENE_LEN;
        string curGene = binary.substr(i, GENE_LEN);
        double var = stoi(curGene, nullptr, 2);
        //������ֵӳ�䵽��Ӧ��range����
        var = bound[index].floor + var / (pow(2, GENE_LEN) - 1) * (bound[index].ceil - bound[index].floor);
        x[index] = static_cast<double>(var);
        //cout << "x" << index << " = " << x[index] << endl;
    }
    return x;
}

/*
* �ŵ���������
* @param indv: �����
* @param weight: ��Ʒ����
* @param value: ��Ʒ��ֵ
*/
double mudFunc(Individual& indv, vector<Individual>& population, vector<double>& gains) {
    string binary = indv.getBinary();
    vector<double> power(GENE_NUM, 0);

    double fitness = 0.0;

    //�������ƽ���Ϊʮ����
    vector<range> bound = {GENE_NUM, range(0, 0.1)};
    //����ÿ���ڵ�Ĺ��ʷ�Χ
    power = decodeBinary(binary, bound);
    indv.setGeneDec(power);

    //����g * p
    vector<double> gp(power.size());
    for (int i = 0; i < gp.size(); i++) {
        gp[i] = gains[i] * power[i];
    }

    for (int i = 0; i < power.size(); i++) {
        //�����Ѿ�������Ľڵ�ĸ���
        double postNoise = 0.0;
        double theta = 0.4;
        for (int j = power.size() - 1; j > i; j--) {
            postNoise += gp[j];
        }
        postNoise *= theta;

        //������δ������Ľڵ����
        double preNoise = 0.0;
        for (int j = 0; j < i; j++) {
            preNoise += gp[j];
        }

        //����ÿһ���ڵ������������
        fitness += (gains[i] * power[i]) / (0.05 * 2 + postNoise + preNoise);
    }
    fitness = 2 * log2(1 + fitness);
    return fitness;
}


/*
* ��Ӧ�ȼ��㣬������Ⱥ�����и������Ӧ��ֵ
* @param population: ��Ⱥ
*/
void calFitness(vector<Individual>& population, Individual& best, double& f_max, double& f_min) {
    double bestFit = -DBL_MAX;
    double maxFit = -DBL_MAX;
    double minFit = DBL_MAX;
    if (best.getFitness() != 0) {
        bestFit = best.getFitness();
    }
    int bestIdx = 0; //���Ÿ���λ��
    int maxIdx = 0;  //f_maxλ��
    int minIdx = 0;  //f_minλ��
    for (int i = 0; i < population.size(); i++) {
        double fitness = mudFunc(population[i], population, GAINS);
        //��¼���Ÿ���
        bestIdx = (fitness >= bestFit) ? i : bestIdx;
        bestFit = max(bestFit, fitness);
        population[i].setFitness(fitness);
        //��¼���������Ÿ���
        maxIdx = (fitness > maxFit) ? i : maxIdx;
        maxFit = max(maxFit, fitness);
        //��¼������������
        minIdx = (fitness < minFit) ? i : minIdx;
        minFit = min(minFit, fitness);
    }
    //���������Ⱥ���б����Ÿ�����Ӧ�ȸ��ߵĸ��壬����best
    //cout << "best indv = " << best.getFitness() << endl;
    //cout << "bestFit = " << bestFit << endl;
    if (bestFit != best.getFitness()) {
        best = population[bestIdx];
    }
    //������Ӧ�ȣ��Ӵ�С������Ⱥ����������������޸�flag
    sort(population.begin(), population.end(), [](Individual& a, Individual& b) {return a.getFitness() > b.getFitness(); });

    //��¼���Ÿ���
    f_max = maxFit;
    //population[maxIdx].setSpecFlag(1);
    //��¼������
    f_min = minFit;
    //population[minIdx].setSpecFlag(2);
}


/*
* �����ֱ�
* ��ʼ����Ⱥ�в��ֻ�ȫ������
*/
void catastrophe(vector<Individual>& population) {
    //for (int i = POP_SIZE - 1; i >= POP_SIZE - POP_SIZE / 2; i--) {
    for (int i = POP_SIZE - 1; i >= 0; i--) {
        population[i] = Individual();
    }
}

/*
* ���ӱ���
*/
void mutation(vector<Individual>& population, Individual& best) {
    //������ʴ�0.001��0.01�仯
    double mutaRate = 0.001;
    //���۵�ǰ��Ⱥ�������Ÿ�����ƽ���������룬������Ϊ������볤��
    int hammingSum = 0;
    for (auto& x : population) {
        hammingSum += x.getHammingDis(best);
    }
    double avgHamming = static_cast<double>(hammingSum) / POP_SIZE;
    cout << "haming = " << avgHamming << '\n';
    //��ƽ����������С��40% * ���������ʱ���������Ϊ���
    if (avgHamming < 0.4 * CHROM_LEN) {
        mutaRate = 0.005;
    }
    if (avgHamming < 0.2 * CHROM_LEN) {
        mutaRate = 0.01;
    }
    if (avgHamming < 0.1 * CHROM_LEN) {
        mutaRate = 0.1;
        catastrophe(population);
        return;
    }
    //��ƽ����������С��10%ʱ����Ϊ�㷨���������������ֱ�
    /*if (avgHamming <= 0.15 * CHROM_LEN) {
        catastrophe(population);
        return;
    }*/
    //catastrophe(population);


    for (int i = 0; i < POP_SIZE; i++) {
        double pick = srand();
        //cout << "muta pick = " << pick << endl;
        if (pick <= mutaRate) {
            /*vector<qubit> chrom(CHROM_LEN);
            for (int j = 0; j < CHROM_LEN; j++) {
                double seed = srand();
                qubit initQ = qubit(sqrt(seed), sqrt(1 - seed));
                chrom[j] = initQ;
            }
            population[i] = Individual(chrom, 0, "");*/
            for (int j = 0; j < CHROM_LEN; j++) {
                swap(population[i].getChrom()[j].alpha, population[i].getChrom()[j].beta);
            }

        }


        /*for (int j = 0; j < CHROM_LEN; j++) {
            swap(population[i].getChrom()[j].alpha, population[i].getChrom()[j].beta);
            cout << i << " a = " << population[i].getChrom()[j].alpha << endl;
            cout << i << " b = " << population[i].getChrom()[j].beta << endl;
        }*/
    }
}


/*
* �����������ӣ����ھֲ�����
*/
void tabu(vector<Individual>& population, Individual& best) {
    //��¼�����������ǰ������ֵ
    static Individual preElite;
    static int stayFlag = 0;
    if (preElite.getFitness() != best.getFitness()) {
        stayFlag = 0;
    }
    preElite = best;

    //��������ʵ�ʵ�����������
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

    //�������ż�¼����֤���������ǵ��������е�����ֵ
    Individual elite = best;
    //�������ɱ�
    deque<tabuItem> tabuList(35);

    //ȷ���ӵ�����ӵ���Ʒ��������Χ��(-4, 4)
    int flag = -4 + 8 * srand();

    int temp = flag;

    //��best��Ϊ��ʼ�⣬�����ֲ�����
    for (int m = 0; m < gens; m++) {
        cout << "Tabu search ������" << m << '\n';
        string bestBin = best.getBinary();
        double bestFit = best.getFitness();
        vector<tabuItem> tabuTemp;  //�������ڴ洢best���жԻ�λ�ý��

        cout << "stayFlag = " << stayFlag << '\n';

        cout << "Flag = " << temp << '\n';
        //�ϴ�δ�ҵ����Ž⣬���һ����Ʒ
        if (stayFlag != 0) {
            for (auto& x : bestBin) {
                //�����Ƕ���һ����Ʒ���Ƕ���һ����Ʒ
                if (flag > 0) {
                    //�����Ʒ
                    if (x == '0') {
                        x = '1';
                        flag--;
                        if (flag <= 0) {
                            break;
                        }
                    }
                }
                else if (flag < 0) {
                    //������Ʒ
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
            double fitness = mudFunc(best, population, GAINS);
            best.setFitness(fitness);
            bestFit = fitness;
        }
        //����������ǰbest����������
        for (int i = 0; i < CHROM_LEN; i++) {
            if (bestBin[i] == '0') {
                for (int j = 0; j < CHROM_LEN; j++) {
                    if (bestBin[j] == '1') {
                        Individual tempChrom = Individual();
                        swap(bestBin[i], bestBin[j]);
                        tempChrom.setBinary(bestBin);
                        //�ж��Ի������Ӧ��
                        double fitness = mudFunc(population[i], population, GAINS);
                        tempChrom.setFitness(fitness);
                        fitness = fitness - bestFit;
                        tabuItem item = tabuItem(pair<int, int>(i, j), fitness, tempChrom);
                        tabuTemp.push_back(item);
                        //�ָ�best
                        bestBin = best.getBinary();
                    }
                }
            }
        }
        //�����н�������Ӧ�ȸı�Ľ���Ӵ�С����
        sort(tabuTemp.begin(), tabuTemp.end(), [](tabuItem a, tabuItem b) {return a.fitDiff > b.fitDiff; });

        //ѡ������tabuList�����������ĶԻ�����
        for (auto& x : tabuTemp) {
            if (find(tabuList.begin(), tabuList.end(), x) != tabuList.end()) {
                //ͻ�ƽ���
                if (x.fitDiff + bestFit > elite.getFitness()) {
                    best = x.newIndv;
                    elite = best;           //���½�������ֵ
                    //���µ�ѡ�������ɱ�
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
                //���µ�ѡ�������ɱ�
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
    //�����������Ž⸳ֵ��best
    best = elite;
    //������������Ž�û�иı䣬���´ν��н�������ʱ������Ʒ
    if (preElite.getFitness() == best.getFitness()) {
        stayFlag++;
    }
    else {
        stayFlag = 0;
    }
}

/*
* ��ȡ����ͼ�ı�����
* @param edge: ��
* @param path: ͼ��txt·��
*/
int readEdge(vector<pair<int, int>>& edge, string path) {
    ifstream infile;
    infile.open(path, ios::in);
    if (!infile.is_open())
    {
        cout << "��ȡ�ļ�ʧ��" << endl;
        return 0;
    }
    //��ȡÿһ���ߣ�������edge
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
* ��ӡ��Ⱥ
*/
void printPopulation(vector<Individual>& pop) {
    for (int i = 0; i < pop.size(); i++) {
        cout << pop[i].toString() << endl;
    }
}

/*
* QGA������
*/
void quantumAlgorithm() {
    Individual best;
    vector<Individual> population;

    //��ʼ����Ⱥ
    initPop(POP_SIZE, CHROM_LEN, population, 0);
    //����Ⱥ����һ�β������õ������Ʊ���
    collapse(population);
    //������Ӧ�ȣ��ҳ����Ÿ�����ٵ�������
    double f_max = 0;
    double f_min = 0;
    calFitness(population, best, f_max, f_min);
    //��������
    int flag = 50;
    for (int gen = 0; gen < MAX_GEN; gen++) {
        cout << "��ǰ���������� " << gen << endl;
        //������Ⱥ
        collapse(population);
        //������Ӧ��
        double f_max = 0;
        double f_min = 0;
        calFitness(population, best, f_max, f_min);
        //������ת��
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
            //tabu(population, best);
        }
        mutation(population, best);
    }
    cout << "--------------------------\n";
    //printPopulation(population);
}

int main()
{
    //�ڵ�����
    vector<pair<double, double>> coords;

    //�ڵ㵽�ն˵ľ���
    DISTANCE = { 16.2182308193243,
                 79.4284540683907,
                 31.1215042044805,
                 52.8533135506213,
                 16.5648729499781,
                 60.1981941401637,
                 26.2971284540144,
                 65.4079098476782,
                 68.9214503140008,
                 74.8151592823709 };
    GAINS = {
        0.130052951147251,
        0.0807473766577741,
        0.106955489777367,
        0.0912430068595769,
        0.129230439684949,
        0.0877498162562804,
        0.112499024440491,
        0.0855918026909817,
        0.0842587345626009,
        0.0822099479526409
    };
    srand(unsigned(time(NULL)));
    quantumAlgorithm();

    std::cout << "Hello World!\n";
}
