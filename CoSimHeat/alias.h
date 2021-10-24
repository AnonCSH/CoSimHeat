#ifndef ALIAS_H
#define ALIAS_H

#include <random>
#include <algorithm>
#include <stack>
#include "Random.h"
using namespace std;

class Alias{
public:
	double* p;
	int* h;
	pair<int, int>* map1;
	int n;
	Alias(vector<pair<pair<int, int>, double> > pi){
		double sum = 0;
		n = pi.size();
		stack<int> xiao;
		stack<int> big;
		p = new double[n];
		h = new int[n];
		map1 = new pair<int, int>[n];
		for(int i = 0; i < n; i++){
			sum += pi[i].second;
			map1[i] = pi[i].first;
		}
		for(int i = 0; i < n; i++){
			p[i] = pi[i].second * n / sum;
			if(p[i] > 1)
				big.push(i);
			else
				xiao.push(i);
		}
		while(!xiao.empty() && !big.empty()){
			int xiaoId = xiao.top();
			xiao.pop();
			int bigId = big.top();
			h[xiaoId] = bigId;
			p[bigId] -= (1-p[xiaoId]);
			if(p[bigId] < 1){
				xiao.push(bigId);
				big.pop();
			}
		}
	}

	~Alias(){
		delete[] p;
		delete[] h;
		delete[] map1;
	}
	pair<int, int> generateRandom(Random& R){
		int firstId = R.drand() * n;
		pair<int, int> answer = R.drand() < p[firstId] ? map1[firstId] : map1[h[firstId]];
		return answer;
	}
};

#endif