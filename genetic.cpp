/*author�����ͳ�
creatTime��2017/6/4
function:ʹ���Ŵ��㷨�Ż��ۿڵĴ�����������*/

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

//λ������
struct position{
	int row;
	int col;
	bool operator<(const position& b)const {
		return col < b.col || (col == b.col && row < b.row);
	}
};

//����Ⱦɫ�����ͺ���Ⱥ����
typedef vector<int> chromType;
typedef vector< vector<int> > populationType;

//���ĵ���ʱ�䣬����ʱ�䣬�����ȵ���Ϣ
struct vesselInfo{
	vector<int> arrTime;
	vector<int> serTime;
	vector<int> berOccu;
};

//��λ�ġ���С����������ʾ��λ�ĳ��ȣ�������ʾ�ܵķ���ʱ��
struct berthSize{
	int rowSize;
	int colSize;
};

//�����˴��ķ���λ�ã��Լ�lastDeparureTime��totalWaitingTime������
struct mySolution{
	vector<position> place;
	int lastDeparureTime;
	int totalWaitingTime;
};


/*�ж������vesselID���Ҵ�����pos���λ�ÿ��в����С�������ĳ��ȣ���Ҫ�ķ���ʱ�����*/
bool isAva(position& pos, berthSize& berSize, vesselInfo& info, map<position, int>& usedBerth, int vesselID){
	bool ava = true;
	//�������posλ����Ҫռ�õ���λ�ã�����Щλ�ý��м�⣬���Ƿ�ռ�ã��Ƿ񳬳���λ���Ȼ�λ�ܵķ���ʱ��
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

/*���vesselID���Ҵ��Ѿ�ռ�е�����*/
void markUsed(position& pos, berthSize& berSize, vesselInfo& info, map<position, int>& usedBerth, int vesselID) {
	for (int c = pos.col; c < pos.col + info.serTime[vesselID]; c++) {
		for (int r = pos.row; r < pos.row + info.berOccu[vesselID]; r++) {
			usedBerth[position{ r,c }] = 1;
		}
	}
}


/*̰�ĵķ�����Ⱦɫ����н��룬�õ�Ⱦɫ���Ӧ�Ĵ��ڷ�λ�ã��ܵȴ�ʱ��ȵȽ��
berSize��		��λ�ġ���С����������ʱ��Ͳ�λ����
info��			���д��ĵȴ�ʱ�䣬����ʱ�䣬��������Ϣ
chrom:			Ҫ����̰�Ĳ�����Ⱦɫ��*/
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


/*ͨ��ĳһ�����н���̰�ģ��õ�������н���̰��ʱ���մ��ڷ�λ�ý�������Ĵ����У�Ҳ����Ⱦɫ�壩����������һ������ͬ��ֻ����������Ͳ�ͬ
berSize��		��λ�ġ���С����������ʱ��Ͳ�λ����
info��			���д��ĵȴ�ʱ�䣬����ʱ�䣬��������Ϣ
chrom:			Ҫ����̰�Ĳ�����Ⱦɫ��*/
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



//����Ⱦɫ������Ӧ��fxֵ�ĺ���
int getFx(vector<int>& chrom, berthSize& berSize, vesselInfo& info){
	//����̰�ĵķ�����Ⱦɫ����н��룬�õ�Ⱦɫ���Ӧ���ܵȴ�ʱ��ȵ���Ϣ
	mySolution msolu = greedy(berSize, info, chrom);
	vector<position> solution = msolu.place;
	int unassigned = chrom.size() - solution.size();
	int lastDeparureTime = msolu.lastDeparureTime;
	int totalWaitingTime = msolu.totalWaitingTime;
	return 100 * unassigned + 2 * totalWaitingTime + lastDeparureTime;
}

//���弸��Ⱦɫ��ı��췽ʽ�������Ŵ��㷨������GA�������һ������method
const int swapGeneMethod = 1;
const int reverseMethod = 2;
const int swapChromMethod = 3;
const int insertMethod = 4;
const int hybridMethod = 5;

/*�Ŵ��㷨���ĺ���
vesselCount��	����������Ҳ��һ��Ⱦɫ���ڵĻ�����Ŀ
initSize��		��Ⱥ�ĳ�ʼ��С�����ڴ˴�ʹ��ÿһ������ͬ�������ķ��������ÿһ������Ⱥ��С��ΪinitSize
maxGeneration��	��ֳmaxGeneration��֮���ٷ�ֳ�����ؽ��
variationRate��	�����ʣ��� 0~1֮�䣬������Խ�ߣ��Ӵ�Ⱦɫ��ı��������Խ��
berSize��		��λ�ġ���С����������ʱ��Ͳ�λ����
info��			���д��ĵȴ�ʱ�䣬����ʱ�䣬��������Ϣ
best_chrom��		��ǰ�ҵ�����ý���������Ⱦɫ�嵱�Ŵ�����֮��ᱻ�޸�Ϊ����Ⱦɫ��
killGap��		ÿ��killGap�������Ͱ���Ⱥȫ��ɱ�������ɵ�ǰ��õĸ���
method��			����ķ������н������򷨣���תȾɫ��Ƭ�η�������Ⱦɫ��Ƭ�η���������򷨣��Լ���Ϸ������ڸ������ѡ��ǰ���ַ�����
*/
mySolution GA(int vesselCount, int initSize, int maxGeneration, double variationRate, berthSize& berSize, 
	vesselInfo& info, vector<int>& best_chrom, int killGap, int method) {
	populationType population;
	map<int, int> isGeneGenerated;
	map< chromType , int> isChromGenerated;
	mySolution best_solution;
	int minFitness = 999999999;
	
	//��ʼ��Ⱥȫ������̰�ĵõ��ĸ���
	int generation = 0;
	while (generation < maxGeneration) {
		if (generation % killGap == 0) {
			population.clear();
			for (int s = 0; s < initSize*2; s++) {
				population.push_back(best_chrom);
			}
		}
		generation++;
		//����ÿ�������fxֵ��Ϊ֮������̶������ݵ�׼��
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
		//���̶Ľ��ж����Ƹ����������������Ϊԭ����Ⱥ��С���൱����ѡ��һ�����������
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

		int crossNum = 0;//�Ѿ���������ĸ���ĸ�����
		set<int> crossParents;
		populationType childPopulation;
		while (crossNum < goodPopulation.size()) {

			//����������
			//�����ý���ı߽�
			crossNum += 2;
			int crossPoint1 = rand() % vesselCount;
			int crossPoint2 = rand() % vesselCount;
			while (crossPoint2 == crossPoint1) {
				crossPoint2 = rand() % vesselCount;
			}
			if (crossPoint1 > crossPoint2)
				swap(crossPoint1, crossPoint2);

			//�����ñ���Ҫ����ĸ�ĸ��壬����ȡ֮ǰ�Ѿ�������ĸ��,�Ѿ�������ĸ�ķ���crossParents����
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


			set<int> section1;//ȡ��Ҫ�����Ƭ��
			set<int> section2;
			for (int i = crossPoint1; i <= crossPoint2; i++)
			{
				section1.insert(goodPopulation[chrom1][i]);
				section2.insert(goodPopulation[chrom2][i]);
			}
			vector<int> childChorm1;
			vector<int> childChorm2;
			//����ǰ��Ĳ���
			for (int i = 0; i < crossPoint1; i++)
			{
				childChorm1.push_back(goodPopulation[chrom1][i]);
				childChorm2.push_back(goodPopulation[chrom2][i]);
			}
			//������ѡƬ�ε��ӽ�
			for (int p = 0; p < vesselCount; p++)
			{
				if (section1.count(goodPopulation[chrom2][p]))
					childChorm1.push_back(goodPopulation[chrom2][p]);
				if (section2.count(goodPopulation[chrom1][p]))
					childChorm2.push_back(goodPopulation[chrom1][p]);
			}
			//���ƺ���Ĳ���
			for (int i = crossPoint2 + 1; i < vesselCount; i++)
			{
				childChorm1.push_back(goodPopulation[chrom1][i]);
				childChorm2.push_back(goodPopulation[chrom2][i]);
			}
			//�������������Ӵ���������Ӵ�����Ⱥ����
			childPopulation.push_back(childChorm1);
			childPopulation.push_back(childChorm2);
		}
		


		//����Ĳ���
		int kind = method;
		if (method == hybridMethod)
			kind = (rand() % 4) + 1;
		if (kind == 1)//��������ı��취
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
		else if (kind == 2)//��תĳ��Ⱦɫ��Ƭ�εı��취
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
		else if (kind == 3)//������������Ⱦɫ��ı��췽��
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
		else if (kind == 4) //��ĳ��������뵽�𴦵ı��췽��
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

		//���Ӵ��͸����ŵ�ͬһ����Ⱥ��Ȼ��while������һ��ѭ����Ҳ���ǽ�����һ�η�ֳ
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
	while (cout << "------------------------------------------\n������Ҫ�����Ϸ��ţ� ", cin >> fileID) {
		//ʹ���ļ�����������
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

		//����һ�������ͣ����С�Ǳ�ŵĴ�С�����ڽ��õ��Ľⰴ�����ŵ�˳������
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
			ans.insert(standAns{ chrom[k], position{ 0,0 } });//��һ��forѭ���м�һ֮����-1
		}

		//����������ŵ�˳�������Ĵ��ڷ�λ��
		set<standAns>::iterator it;
		for (it = ans.begin(); it != ans.end(); it++) {
			cout << it->pos.row - 1 << "," << it->pos.col - 1 << ";";
		}

		cout << endl;
		cout << "unassigned = " << unassigned << " lastDeparureTime = " << lastDeparureTime << " totalWaitingTime = " << totalWaitingTime << endl;
		cout << "Time used:" << endTime - startTime << endl;
	}
}

