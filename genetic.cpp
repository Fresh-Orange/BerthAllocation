/*author：赖贤城
creatTime：2017/6/4
function:使用遗传算法优化港口的船舶调度问题*/

#include<iostream>
#include<fstream>
#include<cmath>
#include<set>
#include<vector>
#include<queue>
#include<string>
#include<functional>
#include<map>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<algorithm>

using namespace std;

//位置类型
struct position{
	int row;
	int col;
	bool operator<(const position& b)const {
		return col < b.col || (col == b.col && row < b.row);
	}
};

//定义染色体类型和种群类型
typedef vector<int> chromType;
typedef vector< vector<int> > populationType;

//船的到达时间，服务时间，船长度等信息
struct vesselInfo{
	vector<int> arrTime;
	vector<int> serTime;
	vector<int> berOccu;
};

//泊位的“大小”，行数表示泊位的长度，列数表示总的服务时间
struct berthSize{
	int rowSize;
	int colSize;
};

//包含了船的放置位置，以及lastDeparureTime，totalWaitingTime的类型
struct mySolution{
	vector<position> place;
	int lastDeparureTime;
	int totalWaitingTime;
};


/*判断如果把vesselID这艘船放在pos这个位置可行不可行。这跟船的长度，需要的服务时间相关*/
bool isAva(position& pos, berthSize& berSize, vesselInfo& info, map<position, int>& usedBerth, int vesselID){
	bool ava = true;
	//如果放在pos位置需要占用到的位置，对这些位置进行检测，看是否被占用，是否超出泊位长度或泊位总的服务时间
	for (int c = pos.col; c < pos.col + info.serTime[vesselID]; c++) {
		for (int r = pos.row; r < pos.row + info.berOccu[vesselID]; r++) {
			if (r > berSize.rowSize || c > berSize.colSize || usedBerth[position{ r,c }] || c < info.arrTime[vesselID]+1) {
				ava = false;
				break;
			}
		}
		if (!ava)
			break;
	}
	return ava;
}

/*标记vesselID这艘船已经占有的区域*/
void markUsed(position& pos, berthSize& berSize, vesselInfo& info, map<position, int>& usedBerth, int vesselID) {
	for (int c = pos.col; c < pos.col + info.serTime[vesselID]; c++) {
		for (int r = pos.row; r < pos.row + info.berOccu[vesselID]; r++) {
			usedBerth[position{ r,c }] = 1;
		}
	}
}


/*贪心的方法对染色体进行解码，得到染色体对应的船摆放位置，总等待时间等等结果
berSize：		泊位的“大小”，即服务时间和泊位长度
info：			所有船的等待时间，到达时间，船长等信息
chrom:			要进行贪心操作的染色体*/
mySolution greedy(berthSize& berSize, vesselInfo& info, vector<int>& chrom) {
	position last = position{ 0,1 };
	vector<position> solution;
	set<position> posOrder;
	vector<int> unplaced;
	map<position, int> pos2gene;
	int lastDeparureTime = 0;
	int totalWaitingTime = 0;
	map<position, int> usedBerth;
	for (int i = 0; i < chrom.size(); i++) {
		bool placed = false;
		for (int c = last.col; c <= berSize.colSize; c++) {
			for (int r = last.row + 1; r <= berSize.rowSize; r++) {
				if (isAva(position{ r,c }, berSize, info, usedBerth, chrom[i])) {
					//cout << "find" << chrom.size() << "  " << berSize.colSize << " " << berSize.rowSize << endl;
					posOrder.insert(position{ r,c });
					pos2gene[position{ r,c }] = chrom[i];
					solution.push_back(position{ r,c });
					int departTime = c + info.serTime[chrom[i]] - 1;
					lastDeparureTime = lastDeparureTime < departTime ? departTime : lastDeparureTime;
					totalWaitingTime += (c - info.arrTime[chrom[i]] - 1);
					markUsed(position{ r,c }, berSize, info, usedBerth, chrom[i]);
					placed = true;
					break;
				}
			}
			if (placed)
				break;
		}
		
	}
	solution.clear();
	set<position>::iterator it;
	for (it = posOrder.begin(); it != posOrder.end(); it++)
	{
		solution.push_back(*it);
	}
	return mySolution{ solution, lastDeparureTime, totalWaitingTime };
}


/*通过某一个序列进行贪心，得到这个序列进行贪心时候按照船摆放位置进行排序的船序列（也就是染色体），方法与上一方法相同，只是最后结果类型不同
berSize：		泊位的“大小”，即服务时间和泊位长度
info：			所有船的等待时间，到达时间，船长等信息
chrom:			要进行贪心操作的染色体*/
vector<int>  greedyChrom(berthSize& berSize, vesselInfo& info, vector<int>& chrom) {
	position last = position{ 0,1 };
	vector<position> solution;
	set<position> posOrder;
	vector<int> unplaced;
	map<position, int> pos2gene;
	int lastDeparureTime = 0;
	int totalWaitingTime = 0;
	map<position, int> usedBerth;
	for (int i = 0; i < chrom.size(); i++) {
		bool placed = false;
		for (int c = last.col; c <= berSize.colSize; c++) {
			for (int r = last.row + 1; r <= berSize.rowSize; r++) {
				if (isAva(position{ r,c }, berSize, info, usedBerth, chrom[i])) {
					posOrder.insert(position{ r,c });
					pos2gene[position{ r,c }] = chrom[i];
					solution.push_back(position{ r,c });
					int departTime = c + info.serTime[chrom[i]] - 1;
					lastDeparureTime = lastDeparureTime < departTime ? departTime : lastDeparureTime;
					totalWaitingTime += (c - info.arrTime[chrom[i]] - 1);
					markUsed(position{ r,c }, berSize, info, usedBerth, chrom[i]);
					placed = true;
					break;
				}
			}
			if (placed)
				break;
		}
		if (!placed)
			unplaced.push_back(chrom[i]);
	}
	vector<int> gChrom;
	set<position>::iterator it;
	for (it = posOrder.begin(); it != posOrder.end(); it++)
	{
		gChrom.push_back(pos2gene[*it]);
	}
	for (int i = 0; i < unplaced.size(); i++) {
		gChrom.push_back(unplaced[i]);
	}
	return gChrom;
}



//计算染色体所对应的fx值的函数
int getFx(vector<int>& chrom, berthSize& berSize, vesselInfo& info){
	//调用贪心的方法对染色体进行解码，得到染色体对应的总等待时间等等信息
	mySolution msolu = greedy(berSize, info, chrom);
	vector<position> solution = msolu.place;
	int unassigned = chrom.size() - solution.size();
	int lastDeparureTime = msolu.lastDeparureTime;
	int totalWaitingTime = msolu.totalWaitingTime;
	return 100 * unassigned + 2 * totalWaitingTime + lastDeparureTime;
}

//定义几种染色体的变异方式，用于遗传算法函数（GA）的最后一个参数method
const int swapGeneMethod = 1;
const int reverseMethod = 2;
const int swapChromMethod = 3;
const int insertMethod = 4;
const int hybridMethod = 5;

/*遗传算法核心函数
vesselCount：	船的总数，也即一条染色体内的基因数目
initSize：		种群的初始大小。由于此处使用每一代都相同个体数的方法，因此每一代的种群大小都为initSize
maxGeneration：	繁殖maxGeneration代之后不再繁殖，返回结果
variationRate：	变异率，在 0~1之间，变异率越高，子代染色体的变异可能性越大
berSize：		泊位的“大小”，即服务时间和泊位长度
info：			所有船的等待时间，到达时间，船长等信息
best_chrom：		当前找到的最好结果，传入的染色体当遗传结束之后会被修改为最优染色体
killGap：		每隔killGap代数，就把种群全部杀掉，换成当前最好的个体
method：			变异的方法，有交换基因法，反转染色体片段法，交换染色体片段法，插入基因法，以及混合法（基于概率随机选择前四种方法）
*/
mySolution GA(int vesselCount, int initSize, int maxGeneration, double variationRate, berthSize& berSize, 
	vesselInfo& info, vector<int>& best_chrom, int killGap, int method) {
	populationType population;
	map<int, int> isGeneGenerated;
	map< chromType , int> isChromGenerated;
	mySolution best_solution;
	int minFitness = 999999999;
	
	//初始种群全部都是贪心得到的个体
	int generation = 0;
	while (generation < maxGeneration) {
		if (generation % killGap == 0) {
			population.clear();
			for (int s = 0; s < initSize*2; s++) {
				population.push_back(best_chrom);
			}
		}
		generation++;
		//计算每个个体的fx值，为之后的轮盘赌做数据的准备
		vector<int> frontier;
		int last_frontier = 0;
		for (int i = 0; i < population.size(); i++) {
			int fitness = getFx(population[i], berSize, info);
			if (minFitness > fitness) {
				minFitness = fitness;
				cout << "minFitness = " << minFitness << endl;
				best_solution = greedy(berSize, info, population[i]);
				best_chrom = population[i];
				best_chrom = greedyChrom(berSize, info, best_chrom);
				cout << endl;
				for (int i = 0; i < best_chrom.size(); i++)
					cout << best_chrom[i] << " ";
				cout << endl;
			}
			
			frontier.push_back(last_frontier);
			last_frontier += fitness;
		}
		//轮盘赌进行对劣势个体的清除，清除数量为原来种群大小。相当于挑选出一半的优良个体
		set<int> abandon;
		for (int i = 0; i < initSize; i++)
		{
			int ran = rand() % last_frontier;
			for (int f = frontier.size() - 1; f >= 0; f--)
			{
				
				if (frontier[f] < ran) {
					if (abandon.count(f)) {
						i--;
						break;
					}
					abandon.insert(f);
					break;
				}
			}
		}

		populationType goodPopulation = population;

		int crossNum = 0;//已经“当过父母”的个体数
		set<int> crossParents;
		populationType childPopulation;
		while (crossNum < goodPopulation.size()) {

			//交叉产生后代
			//随机获得交叉的边界
			crossNum += 2;
			int crossPoint1 = rand() % vesselCount;
			int crossPoint2 = rand() % vesselCount;
			while (crossPoint2 == crossPoint1) {
				crossPoint2 = rand() % vesselCount;
			}
			if (crossPoint1 > crossPoint2)
				swap(crossPoint1, crossPoint2);

			//随机获得本次要当父母的个体，但不取之前已经当过父母的,已经当过父母的放入crossParents里面
			int chrom1 = rand() % goodPopulation.size();
			while (crossParents.count(chrom1)) {
				chrom1 = rand() % goodPopulation.size();
			}
			crossParents.insert(chrom1);
			int chrom2 = rand() % goodPopulation.size();
			while (crossParents.count(chrom2) && crossNum < goodPopulation.size()-1) {
				chrom2 = rand() % goodPopulation.size();
			}
			crossParents.insert(chrom2);


			set<int> section1;//取出要交叉的片段
			set<int> section2;
			for (int i = crossPoint1; i <= crossPoint2; i++)
			{
				section1.insert(goodPopulation[chrom1][i]);
				section2.insert(goodPopulation[chrom2][i]);
			}
			vector<int> childChorm1;
			vector<int> childChorm2;
			//复制前面的部分
			for (int i = 0; i < crossPoint1; i++)
			{
				childChorm1.push_back(goodPopulation[chrom1][i]);
				childChorm2.push_back(goodPopulation[chrom2][i]);
			}
			//进行所选片段的杂交
			for (int p = 0; p < vesselCount; p++)
			{
				if (section1.count(goodPopulation[chrom2][p]))
					childChorm1.push_back(goodPopulation[chrom2][p]);
				if (section2.count(goodPopulation[chrom1][p]))
					childChorm2.push_back(goodPopulation[chrom1][p]);
			}
			//复制后面的部分
			for (int i = crossPoint2 + 1; i < vesselCount; i++)
			{
				childChorm1.push_back(goodPopulation[chrom1][i]);
				childChorm2.push_back(goodPopulation[chrom2][i]);
			}
			//将产生的两个子代个体放入子代的种群当中
			childPopulation.push_back(childChorm1);
			childPopulation.push_back(childChorm2);
		}
		


		//变异的操作
		int kind = method;
		if (method == hybridMethod)
			kind = (rand() % 4) + 1;
		if (kind == 1)//交换基因的变异法
			for (int i = 0; i < childPopulation.size(); i++)
			{
				if ((rand() % 1000)*1.0 / 1000 < variationRate) {
					int crossPoint1 = rand() % vesselCount;
					int crossPoint2 = rand() % vesselCount;
					while (crossPoint2 == crossPoint1) {
						crossPoint2 = rand() % vesselCount;
					}
					swap(childPopulation[i][crossPoint1], childPopulation[i][crossPoint2]);
				}
			}
		else if (kind == 2)//反转某段染色体片段的变异法
			for (int i = 0; i < childPopulation.size(); i++)
			{
				if ((rand() % 1000)*1.0 / 1000 < variationRate) {
					int crossPoint1 = rand() % vesselCount;
					int crossPoint2 = rand() % vesselCount;
					while (crossPoint2 == crossPoint1) {
						crossPoint2 = rand() % vesselCount;
					}
					if (crossPoint1 > crossPoint2)
						swap(crossPoint1, crossPoint2);
					reverse(childPopulation[i].begin() + crossPoint1, childPopulation[i].begin() + crossPoint2);
				}
			}
		else if (kind == 3)//交换左右两段染色体的变异方法
			for (int i = 0; i < childPopulation.size(); i++)
			{
				if ((rand() % 1000)*1.0 / 1000 < variationRate) {
					int crossPoint = rand() % vesselCount;
					vector<int> childCopy = childPopulation[i];
					childPopulation[i].clear();
					for (int k = crossPoint; k < childCopy.size(); k++) {
						childPopulation[i].push_back(childCopy[k]);
					}
					for (int k = 0; k < crossPoint; k++) {
						childPopulation[i].push_back(childCopy[k]);
					}
				}
			}
		else if (kind == 4) //将某个基因插入到别处的变异方法
			for (int i = 0; i < childPopulation.size(); i++)
			{
				if ((rand() % 1000)*1.0 / 1000 < variationRate) {
					int crossPoint1 = rand() % vesselCount;
					int crossPoint2 = rand() % vesselCount;
					while (crossPoint2 == crossPoint1) {
						crossPoint2 = rand() % vesselCount;
					}
					int gene2 = childPopulation[i][crossPoint2];
					childPopulation[i].erase(childPopulation[i].begin() + crossPoint2);
					childPopulation[i].insert(childPopulation[i].begin() + crossPoint1, gene2);					
				}
			}

		//将子代和父代放到同一个种群，然后while进行下一次循环，也就是进行下一次繁殖
		population.clear();
		population = childPopulation;
		for (int i = 0; i < goodPopulation.size(); i++)
		{
			population.push_back(goodPopulation[i]);
			if (population.size() >= initSize*2)
				break;
		}
	}
	return best_solution;
}

int main() {
	
	string fileID;
	while (cout << "------------------------------------------\n请输入要玩的游戏编号： ", cin >> fileID) {
		//使用文件流进行输入
		string fileName = "F:\\C++\\GAtest\\" + fileID + ".txt";
		ifstream cif(fileName);
		srand((unsigned)time(NULL));
		int c, r, n;
		int arr, ser, ber;
		cif >> c;
		cif >> r;
		cif >> n;
		berthSize berSize = berthSize{ r,c };
		vesselInfo info;
		info.arrTime.push_back(arr);
		info.serTime.push_back(ser);
		info.berOccu.push_back(ber);
		for (int i = 0; i < n; i++) {
			cif >> arr >> ser >> ber;
			info.arrTime.push_back(arr);
			info.serTime.push_back(ser);
			info.berOccu.push_back(ber);
		}

		chromType chrom;
		for (int i = 1; i < n + 1; i++) {
			chrom.push_back(i);
		}
		clock_t startTime = clock();
		mySolution msolu = GA(n, 3, 5, 0.9, berSize, info, chrom, 5, swapGeneMethod);
		clock_t endTime = clock();
		vector<position> solution = msolu.place;
		int unassigned = chrom.size() - solution.size();
		int lastDeparureTime = msolu.lastDeparureTime;
		int totalWaitingTime = msolu.totalWaitingTime;
		cout << "best_chrom: ";
		for (int i = 0; i < chrom.size(); i++) {
			cout << chrom[i] << ", ";
		}
		cout << endl;

		//定义一个新类型，其大小是编号的大小，用于将得到的解按照其编号的顺序排序
		struct standAns {
			int index;
			position pos;
			bool operator<(const standAns& b) const {
				return index < b.index;
			}
		};
		set<standAns> ans;
		int i;
		for (i = 0; i < solution.size(); i++) {
			ans.insert(standAns{ chrom[i], solution[i] });
		}
		for (int k = i; k < chrom.size(); k++) {
			ans.insert(standAns{ chrom[k], position{ 0,0 } });//下一个for循环中减一之后变成-1
		}

		//输出按照其编号的顺序排序后的船摆放位置
		set<standAns>::iterator it;
		for (it = ans.begin(); it != ans.end(); it++) {
			cout << it->pos.row - 1 << "," << it->pos.col - 1 << ";";
		}

		cout << endl;
		cout << "unassigned = " << unassigned << " lastDeparureTime = " << lastDeparureTime << " totalWaitingTime = " << totalWaitingTime << endl;
		cout << "Time used:" << endTime - startTime << endl;
	}
}

