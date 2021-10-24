#pragma once
#ifndef EXACTSIM_H
#define EXACTSIM_H
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unordered_map>
#include "alias.h"
#include "stdafx.h"
#include "algorithm.h"
using namespace std;
using namespace Eigen;

class ExactSim : public algorithm{
public:
	const string name = "ExactSim";

	int kmax;
	double c;
	string resultFileEN;
	string resultFileDP;
	SparseMatrix<double> mat;

	PNGraph graph;
	int vert;
	double level;
	Random R;
	double eps;
	double delta;
	int* vecpos;
	double* reserve;
	double* residue;
	double* newReserve;
	double* newResidue;
	bool* isInArray;
	long power_walk_num;
	int power_ite_num;
	unordered_map<int, unordered_map<int, double> > P_lw;//transition probablity used in function cal_Zw()
	unordered_map<int, unordered_map<int, unordered_map<int, double> > > P_ly;//<final_node<level,<levelnode,pr>>>
	unordered_map<int, unordered_map<int, vector<pair<int, double> > > > P_ally;//<level <node <final_node, pr> > >
	unordered_map<int, unordered_map<int, double> > Fl;//used in function cal_Zw()
	unordered_map<int, unordered_map<int, double> > answer;
	vector<vector<int> > forwardNode;
	vector<int> pathnum_set;
	vector<int> new_pathnum_set;
	vector<int> candidate_set_ford;
	vector<int> new_candidate_set_ford;

	string getName() { return name; }
	string getResultFile(string scene);

	ExactSim(SparseMatrix<double>& M, PNGraph& graph, double decay, int iterations, bool isWeigthed);
	~ExactSim();

	//generate transfer matrix
	void getTransferMatrix(SparseMatrix<double>& mat);
	//single-source
	VectorXd singleSource(int u);

	// SYNONYM EXTRACTION
	VectorXd SYN_SingleSource(int id);

	//calculate transition probablity
	unordered_map<int, unordered_map<int, double> > forwardPushByLevel(int u, int maxLevel, int& deterlevel,
		long walk_num, vector<vector<int> >& forwardNode);
	void cal_Zw(int w, int& deterlevel, int maxlevel, long total_walknum, double& pl_value);
	double sampleD(int nodeId, double ppr_nodeId, long walk_num, int maxlevel, int& deterlevel, double para_r);
	SparseMatrix<double> estimationMatrixD(int u, double para_r);
};

ExactSim::ExactSim(SparseMatrix<double>& M, PNGraph& graph, double decay, int iterations, bool isWeighted) {
	mat = M;
	kmax = iterations;
	c = decay;
	resultFileEN = "result/ExactSim/ts68.txt";
	resultFileDP = "result/ExactSim/DP.txt";

	this->graph = graph;
	vert = mat.rows();
	level = 1000;
	eps = 1e-5;
	delta = 0.01;
	reserve = new double[vert];
	residue = new double[vert];
	newReserve = new double[vert];
	newResidue = new double[vert];
	isInArray = new bool[vert];
	vecpos = new int[vert];
	power_walk_num = (long)((c * 0.00001) * (log(vert / delta) / log(2)) / (eps * eps));
	power_ite_num = (int)((log(eps * (1 - c)) / log(c) - 1) + 1);
	for (int i = 0; i < vert; i++) {
		isInArray[i] = false;
		residue[i] = 0;
		newResidue[i] = 0;
		reserve[i] = 0;
		newReserve[i] = 0;
		vecpos[i] = 0;
	}

	if (!isWeighted) {
		getTransferMatrix(mat);
	}
}

ExactSim::~ExactSim() {
	delete[] reserve;
	delete[] residue;
	delete[] newReserve;
	delete[] newResidue;
	delete[] isInArray;
	delete[] vecpos;
	vector<vector<int> >().swap(forwardNode);
	vector<int>().swap(pathnum_set);
	vector<int>().swap(new_pathnum_set);
	vector<int>().swap(candidate_set_ford);
	vector<int>().swap(new_candidate_set_ford);
}

string ExactSim::getResultFile(string scene) {
	if (scene == "EN") {
		return resultFileEN;
	}
	return resultFileDP;
}

void ExactSim::getTransferMatrix(SparseMatrix<double>& mat) {
	for (int k = 0; k < mat.outerSize(); ++k) {
		double sum = mat.col(k).sum();
		for (SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
		{
			it.valueRef() /= sum;
		}
	}
}

//calculate transition probablity
unordered_map<int, unordered_map<int, double> > ExactSim::forwardPushByLevel(int u, int maxLevel, int& deterlevel,
	long walk_num, vector<vector<int> >& forwardNode) {
	answer.clear();
	forwardNode.clear();
	pathnum_set.clear();
	new_pathnum_set.clear();

	residue[u] = 1;
	candidate_set_ford.push_back(u);
	pathnum_set.push_back(1);
	int tempLevel = 0;
	bool deter_flag = false;

	while (1) {
		long candidateSetSize = candidate_set_ford.size();
		if (candidateSetSize <= 0) {
			deterlevel = tempLevel - 1;
			break;
		}
		long current_pathnum = 0;
		long residue_count = 0;
		if (tempLevel > 1) {
			for (int j = 0; j < candidateSetSize; j++) {
				int tempNode = candidate_set_ford[j];
				if (graph->IsNode(tempNode)) {
					TNGraph::TNodeI NI_tmp = graph->GetNI(tempNode);
					int inSize = NI_tmp.GetInDeg();
					current_pathnum += pathnum_set[j] * inSize;
					if (current_pathnum * current_pathnum >= walk_num) {
						deterlevel = tempLevel;
						deter_flag = true;
						break;;
					}
				}
			}
		}
		for (int j = 0; j < candidateSetSize; j++) {
			int tempNode = candidate_set_ford[j];
			answer[tempLevel][tempNode] = residue[tempNode];
			if (tempLevel == maxLevel || deter_flag == true) {
				residue[tempNode] = 0;
				continue;
			}
			TNGraph::TNodeI NI_tmp = graph->GetNI(tempNode);
			int inSize = NI_tmp.GetInDeg();
			for (int k = 0; k < inSize; k++) {
				int newNode = NI_tmp.GetInNId(k);
				newResidue[newNode] += residue[tempNode] / (double)inSize;
				if (!isInArray[newNode]) {
					isInArray[newNode] = true;
					new_candidate_set_ford.push_back(newNode);
					new_pathnum_set.push_back(pathnum_set[j]);
					vecpos[newNode] = residue_count;
					residue_count += 1;
				}
				else {
					new_pathnum_set[vecpos[newNode]] += pathnum_set[j];
				}
			}
			residue[tempNode] = 0;
		}

		long new_candidateSetSize = new_candidate_set_ford.size();
		for (int j = 0; j < new_candidateSetSize; j++) {
			int tempNode = new_candidate_set_ford[j];
			residue[tempNode] = newResidue[tempNode];
			newResidue[tempNode] = 0;
			isInArray[tempNode] = false;
			vecpos[tempNode] = 0;
		}

		forwardNode.push_back(candidate_set_ford);

		candidate_set_ford.swap(new_candidate_set_ford);
		new_candidate_set_ford.clear();
		pathnum_set.swap(new_pathnum_set);
		new_pathnum_set.clear();

		tempLevel++;
		if (tempLevel > maxLevel) {
			deterlevel = maxLevel;
			break;
		}
		if (deter_flag == true) {
			break;
		}
	}

	candidate_set_ford.clear();
	pathnum_set.clear();

	return answer;
}

//calculate Z(w,q)
void ExactSim::cal_Zw(int w, int& deterlevel, int maxlevel, long total_walknum, double& pl_value) {
	pl_value = 0.0;
	P_lw.clear();//unordered_map <int,unordered_map<int,double> >
	P_ly.clear();//unordered_map <int,unordered_map<int,double> >
	P_ally.clear();//unordered_map <int,unordered_map<int, vector<pair<int,double> > > >
	Fl.clear();//unordered_map <int,unordered_map<int,double> >
	forwardNode.clear();

	P_lw = forwardPushByLevel(w, maxlevel, deterlevel, total_walknum, forwardNode);

	if (deterlevel <= 1) {
		deterlevel = 1;
		if (graph->IsNode(w)) {
			TNGraph::TNodeI NI_w = graph->GetNI(w);
			int inSize = NI_w.GetInDeg();
			if (inSize > 0)
				pl_value = c / inSize;
		}
		return;
	}
	for (int level = 1; level <= deterlevel; level++) {
		for (int ilevel = level; ilevel >= 1; ilevel--) {
			for (int itey = 0; itey < forwardNode[ilevel].size(); itey++) {
				int y = forwardNode[ilevel][itey];
				if (ilevel == level) {
					P_ally[level - ilevel][y].push_back(make_pair(y, 1));
					P_ly[y][level - ilevel][y] = 1;
				}
				else {
					unordered_map<int, int> isinvec;
					int cnt = 0;
					TNGraph::TNodeI NI_y = graph->GetNI(y);
					int inSize = NI_y.GetInDeg();
					for (int k = 0; k < inSize; k++) {
						int nextNode = NI_y.GetInNId(k);
						for (int veci = 0; veci < P_ally[level - ilevel - 1][nextNode].size(); veci++) {
							double prob = P_ally[level - ilevel - 1][nextNode][veci].second / NI_y.GetInDeg();
							int finalNode = P_ally[level - ilevel - 1][nextNode][veci].first;
							if (isinvec.find(finalNode) == isinvec.end()) {
								P_ally[level - ilevel][y].push_back(make_pair(finalNode, prob));
								isinvec[finalNode] = P_ally[level - ilevel][y].size() - 1;
								cnt += 1;
								P_ly[finalNode][level - ilevel][y] = prob;
							}
							else {
								int pos = isinvec[finalNode];
								P_ally[level - ilevel][y][pos].second += prob;
								P_ly[finalNode][level - ilevel][y] += prob;
							}
						}
					}
				}
			}
		}
		for (int itex = 0; itex < forwardNode[level].size(); itex++) {
			int x = forwardNode[level][itex];
			double Pr_repeatmeet = 0.0;
			for (int jlevel = 1; jlevel <= level - 1; jlevel++) {
				for (int itey = 0; itey < forwardNode[jlevel].size(); itey++) {
					int y = forwardNode[jlevel][itey];
					if (P_ly[x][level - jlevel].find(y) != P_ly[x][level - jlevel].end()) {
						double tmp_ply = P_ly[x][level - jlevel][y];
						Pr_repeatmeet += pow(c, level - jlevel) * tmp_ply * tmp_ply * Fl[jlevel][y];
					}
				}
			}
			Fl[level][x] = pow(c, level) * P_lw[level][x] * P_lw[level][x] - Pr_repeatmeet;
			if (Fl[level][x] < 0.01 * eps) {
				Fl[level][x] = 0.0;
			}
			pl_value += Fl[level][x];
		}
		P_ly.clear();
		P_ally.clear();
	}

	return;
}

double ExactSim::sampleD(int nodeId, double ppr_nodeId, long walk_num, int maxlevel, int& deterlevel, double para_r) {
	double meet = 0;
	double pl = 0.0;
	deterlevel = 1;

	if (walk_num == 0.0) {
		deterlevel = 1;
		TNGraph::TNodeI NI_n = graph->GetNI(nodeId);
		int inSize = NI_n.GetInDeg();
		if (inSize > 0)
			return 1 - c / inSize;
	}
	else {
		bool first_lflag = false;
		long firstlevelnum = 0;
		TNGraph::TNodeI NI_n = graph->GetNI(nodeId);
		int winSize = NI_n.GetInDeg();
		if (winSize * winSize >= para_r * walk_num) {
			deterlevel = 1;
			first_lflag = true;
		}
		else {
			for (int li = 0; li < winSize; li++) {
				int tmpN = NI_n.GetInNId(li);
				TNGraph::TNodeI NI_tN = graph->GetNI(tmpN);
				firstlevelnum += NI_tN.GetInDeg();
				if (firstlevelnum * firstlevelnum >= para_r * walk_num) {
					deterlevel = 1;
					first_lflag = true;
					break;
				}
			}
			if (first_lflag == false) {
				pl = 0.0;
				cal_Zw(nodeId, deterlevel, maxlevel, (long)(para_r * walk_num), pl);
			}
		}
		for (long i = 0; i < walk_num; i++) {
			bool flag = true;
			int u_newNode = nodeId, v_newNode = nodeId, u_nextNode, v_nextNode;
			int count = 0;
			if (winSize < 2) {
				break;
			}
			u_newNode = NI_n.GetInNId(R.generateRandom() % winSize);
			v_newNode = NI_n.GetInNId(R.generateRandom() % winSize);
			for (int l = 1; l < deterlevel; l++) {
				if (u_newNode == v_newNode) {
					flag = false;
					break;
				}
				TNGraph::TNodeI NI_u = graph->GetNI(u_newNode);
				TNGraph::TNodeI NI_v = graph->GetNI(v_newNode);
				if ((NI_u.GetInDeg() == 0) || (NI_v.GetInDeg() == 0)) {
					flag = false;
					break;
				}
				int u_tmpNode, v_tmpNode;
				u_tmpNode = NI_u.GetInNId((int)(R.generateRandom() % NI_u.GetInDeg()));
				v_tmpNode = NI_v.GetInNId((int)(R.generateRandom() % NI_v.GetInDeg()));
				u_newNode = u_tmpNode;
				v_newNode = v_tmpNode;
			}
			if ((u_newNode != v_newNode) && (flag == true)) {
				while (R.drand() < c) {
					TNGraph::TNodeI NI_u = graph->GetNI(u_newNode);
					TNGraph::TNodeI NI_v = graph->GetNI(v_newNode);
					int length = NI_u.GetInDeg();
					if (length == 0)
						break;
					int r = R.generateRandom() % length;
					u_nextNode = NI_u.GetInNId(r);
					length = NI_v.GetInDeg();
					if (length == 0)
						break;
					r = R.generateRandom() % length;
					v_nextNode = NI_v.GetInNId(r);
					if (u_nextNode == v_nextNode) {
						meet += 1.0;
						break;
					}
					u_newNode = u_nextNode;
					v_newNode = v_nextNode;
				}
			}
		}

		meet /= walk_num;

		if (deterlevel <= 1) {
			return (1 - c / NI_n.GetInDeg() - c * meet);
		}
		else {
			return (1 - pl - (pow(c, deterlevel)) * meet);
		}
	}
}

SparseMatrix<double> ExactSim::estimationMatrixD(int u, double para_r) {
	SparseMatrix <double> D(vert, vert);

	int targetlevel = power_ite_num;//iteration nums
	double* scores = new double[vert];
	double* nextScores = new double[vert];
	double** ulevelscores = new double* [targetlevel];
	double alpha = 1 - sqrt(c);//teleport probability in PPR computation
	long nr_pre = power_walk_num;//the sampling number in total
	int deterlevel;//l(k)
	int max_deterlevel = 5;//upper bound of l(k)  
	double avg_ppr = 0.0;//
	int ppr_nnz = 0;//nonzero of PPR 
	for (int i = 0; i < vert; i++) {
		scores[i] = 0;
		nextScores[i] = 0;
	}

	//PPR computation for source node u
	//
	scores[u] = 1.0;
	for (int i = 0; i < targetlevel; i++) {
		for (int j = 0; j < vert; j++) {
			TNGraph::TNodeI NI_j = graph->GetNI(j);
			for (int k = 0; k < NI_j.GetInDeg(); k++) {
				int tempNode = NI_j.GetInNId(k);
				nextScores[tempNode] += (1 - alpha) * scores[j] / NI_j.GetInDeg();
			}
		}
		nextScores[u] += alpha;
		for (int j = 0; j < vert; j++) {
			TNGraph::TNodeI NI_j = graph->GetNI(j);
			scores[j] = nextScores[j];
			nextScores[j] = 0.0;
			if ((i == targetlevel - 1) && (scores[j] > 0.0) && (NI_j.GetInDeg() > 1)) {
				avg_ppr += scores[j];
				ppr_nnz += 1;
			}
		}
	}
	if (ppr_nnz > 0) {
		avg_ppr /= ppr_nnz;
	}

	//Matrix D's estimation
	//
	for (int j = 0; j < vert; j++) {
		TNGraph::TNodeI NI_j = graph->GetNI(j);
		if (NI_j.GetInDeg() == 0) {
			D.coeffRef(j, j) = 1.0;
			continue;
		}
		else if ((scores[j] == 0) || (NI_j.GetInDeg() == 1)) {
			D.coeffRef(j, j) = 1 - c / (float)NI_j.GetInDeg();
			continue;
		}
		else {
			long nrj_pre = (long)nr_pre * scores[j] / (avg_ppr * ppr_nnz);//the sample number from node vj
			D.coeffRef(j, j) = sampleD(j, scores[j], nrj_pre, max_deterlevel, deterlevel, para_r);
		}
	}

	delete[] scores;
	delete[] nextScores;

	return D;
}

VectorXd ExactSim::SYN_SingleSource(int id) {
	// construct v0 -- vk
	SparseMatrix<double> D = estimationMatrixD(id, level);

	VectorXd ej = VectorXd::Zero(vert);
	ej(id) = 1;
	vector<VectorXd> v(kmax + 1);
	v[0] = ej;
	for (int k = 1; k <= kmax; k++) {
		v[k] = mat * v[k - 1];
	}

	for (int k = 0; k <= kmax; k++) {
		v[k] = D * v[k];
	}

	// compute uk as single-source
	int iter = kmax;
	VectorXd u = v[iter];
	while (iter > 0) {
		u = c * (u.transpose() * mat).transpose() + v[iter - 1];
		iter--;
	}

	return u;
}

VectorXd ExactSim::singleSource(int u) {
	VectorXd scores(vert);
	if (graph->IsNode(u)) {
		TNGraph::TNodeI NI_w = graph->GetNI(u);
		int inSize = NI_w.GetInDeg();
		if (inSize == 0)
			return scores;

		clock_t t_start = clock();

		int targetlevel = power_ite_num;//iteration nums
		double* nextScores = new double[vert];
		double* d = new double[vert];//matrix D
		double** ulevelscores = new double* [targetlevel];
		double* simrank = new double[vert];
		double alpha = 1 - sqrt(c);//teleport probability in PPR computation
		long nr_pre = power_walk_num;//thesampling number in total
		int deterlevel;//l(k)
		int max_deterlevel = 5;//upper bound of l(k)  
		double avg_ppr = 0.0;//
		int ppr_nnz = 0;//nonzero of PPR 
		for (int i = 0; i < vert; i++) {
			simrank[i] = 0;
			scores[i] = 0;
			nextScores[i] = 0;
		}

		//PPR computation for source node u
		scores[u] = 1.0;
		for (int i = 0;i < targetlevel;i++) {
			for (int j = 0;j < vert;j++) {
				TNGraph::TNodeI NI_j = graph->GetNI(j);
				for (int k = 0;k < NI_j.GetInDeg();k++) {
					int tempNode = NI_j.GetInNId(k);
					nextScores[tempNode] += (1 - alpha) * scores[j] / NI_j.GetInDeg();
				}
			}
			nextScores[u] += alpha;
			for (int j = 0;j < vert;j++) {
				TNGraph::TNodeI NI_j = graph->GetNI(j);
				scores[j] = nextScores[j];
				nextScores[j] = 0.0;
				if ((i == targetlevel - 1) && (scores[j] > 0.0) && (NI_j.GetInDeg() > 1)) {
					avg_ppr += scores[j];
					ppr_nnz += 1;
				}
			}
		}
		if (ppr_nnz > 0) {
			avg_ppr /= ppr_nnz;
		}

		//Matrix D's estimation
		for (int j = 0; j < vert; j++) {
			TNGraph::TNodeI NI_j = graph->GetNI(j);
			if (NI_j.GetInDeg() == 0) {
				d[j] = 1;
				continue;
			}
			else if ((scores[j] == 0) || (NI_j.GetInDeg() == 1)) {
				d[j] = 1 - c / (double)NI_j.GetInDeg();
				continue;
			}
			else {
				long nrj_pre = (long)nr_pre * scores[j] / (avg_ppr * ppr_nnz);//the sample number from node vj
				d[j] = sampleD(j, scores[j], nrj_pre, max_deterlevel, deterlevel, this->level);
			}
		}

		//SimRank vectors' linearization
		//    
		for (int i = 0;i < vert;i++) {
			nextScores[i] = 0.0;
		}
		nextScores[u] = 1.0;
		for (int level = 0;level <= targetlevel;level++) {
			ulevelscores[level] = new double[vert];
			for (int i = 0;i < vert;i++) {
				ulevelscores[level][i] = d[i] * nextScores[i];
				scores[i] = 0.0;
			}
			for (int j = 0;j < vert;j++) {
				TNGraph::TNodeI NI_j = graph->GetNI(j);
				for (int k = 0;k < NI_j.GetOutDeg();k++) {
					int tempNode = NI_j.GetOutNId(k);
					TNGraph::TNodeI NI_tmp = graph->GetNI(tempNode);
					if (NI_tmp.GetInDeg() > 0) {
						scores[j] += nextScores[tempNode] / NI_tmp.GetInDeg();
					}
				}
			}
			for (int j = 0;j < vert;j++) {
				nextScores[j] = scores[j];
			}
		}

		for (int j = 0;j < vert;j++) {
			simrank[j] = ulevelscores[targetlevel][j];
		}

		for (int level = 1;level <= targetlevel;level++) {
			for (int i = 0; i < vert; i++) {
				nextScores[i] = 0;
				scores[i] = 0;
			}
			for (int j = 0;j < vert;j++) {
				TNGraph::TNodeI NI_j = graph->GetNI(j);
				for (int k = 0;k < NI_j.GetInDeg();k++) {
					int tempNode = NI_j.GetInNId(k);
					if (NI_j.GetInDeg() > 0) {
						nextScores[j] += simrank[tempNode] / NI_j.GetInDeg();
					}
				}
				int tmplevel = targetlevel - level;
				nextScores[j] = nextScores[j] * c;
				scores[j] = nextScores[j] + ulevelscores[tmplevel][j];
			}
			for (int i = 0;i < vert;i++) {
				simrank[i] = scores[i];
			}
		}
		clock_t t_end = clock();
		double querytime = (t_end - t_start) / (double)CLOCKS_PER_SEC;
		cout << "ExactSim singleSource : " << querytime << endl;

		delete[] nextScores;
		delete[] d;
		delete[] simrank;
		for (int i = 0;i <= targetlevel;i++) {
			delete[] ulevelscores[i];
		}
	}
	return scores;
}

#endif // !EXACTSIM_H

