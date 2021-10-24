#pragma once
#include <iostream>
#include <vector>
#include "tool.h"
#include <cstdlib>
#include "RPCoSim.h"
using namespace std;

vector<int> RP_singleQuery(MatrixXd& S, int vert, string request, string correct, int& synNum) {

	vector<int> rank;
	cout << request << "  begin==================" << endl;

	//get single-sourse
	int id = stoi(request);
	VectorXd res = S.col(id);

	//sort single-source
	vector<pair<int, float> > row_sim;
	for (int i = 0; i < vert; i++) {
		if (res(i) > 0.000000000001)
			row_sim.push_back(pair<int, float>(i, res(i)));
	}
	if (row_sim.size() == 0) {
		rank.push_back(vert);
		return rank;
	}
	sort(row_sim.begin(), row_sim.end(), cmp);

	// compare
	int i = 0;
	int pos = 0;

	vector<string> split_correct;
	splitString(correct, split_correct, ",");
	synNum = split_correct.size();
	for (i = 0; i < row_sim.size(); i++) {
		string word = to_string(row_sim[i].first);

		for (int j = 0; j < split_correct.size(); j++) {
			if (word == split_correct[j]) {
				cout << "request = " << request << "  word = " << word << "  i = " << i - 1 << endl;
				rank.push_back(i - 1);
			}
		}
		if (rank.size() == split_correct.size()) {
			break;
		}
	}
	if (rank.size() == 0) {
		i = row_sim.size() < 10 ? 11 : row_sim.size();
	}

	return rank;
}

void metricDP(vector<vector<int> > allResult, int testsetSize, vector<int> synNumVec) {
	float ndcg = 0;
	float idcg = 0;
	float notFound = 0;
	float kendall = 0;
	float SPCC = 0;
	float MAP = 0;

	for (int i = 0;i < allResult.size();i++) {
		vector<int> result = allResult[i];
		int synNum = synNumVec[i];

		if (result.size() == 0) {
			continue;
		}
		if (result[0] == -1) {
			notFound++;
			continue;
		}

		float tmp = 0;
		int n = 20;
		for (int j = 0;j < result.size();j++) {
			if (result[j] < n)
				tmp += ((float)(n + j - result[j]) / n);
		}
		kendall += tmp / result.size();

		tmp = 0;
		for (int j = 0;j < result.size();j++) {
			if (result[j] < n)
				tmp += pow(result[j] - j, 2);
			else
				tmp += pow(n - j, 2);
		}
		SPCC += 1 - 6 * tmp / (n * (pow(n, 2) - 1));

		float dcg = 0;
		for (int j = 0; j < result.size(); j++) {
			if (result[j] < n)
				dcg += log(2.0) / log(result[j] + 2.0);
		}
		float idcg = 0;
		for (int j = 0; j < synNum; j++) {
			idcg += log(2.0) / log(j + 2.0);
		}
		if (idcg != 0)
			ndcg += dcg / idcg;

		float AP = 0;
		for (int j = 0; j < result.size(); j++) {
			AP += (j + 1.0) / (result[j] + 1.0);
		}
		if (synNum > 0) {
			MAP += AP / (float)synNum;
		}

	}

	cout << "Result:" << endl;
	cout << "MAP: " << setprecision(3) << (MAP / testsetSize) << "/" << setprecision(3) << (MAP / (testsetSize - notFound)) << endl;
	cout << "NDCG: " << setprecision(3) << (ndcg / testsetSize) << "/" << setprecision(3) << (ndcg / (testsetSize - notFound)) << endl;
	cout << "kendall: " << setprecision(3) << (kendall / testsetSize) << "/" << setprecision(3) << (kendall / (testsetSize - notFound)) << endl;
	cout << "SPCC: " << setprecision(3) << (SPCC / testsetSize) << "/" << setprecision(3) << (SPCC / (testsetSize - notFound)) << endl;

}

void RP_DP(RPCoSim rp, vector<string> testset, int vert) {
	cout << endl << "============ " << rp.getName() << " ============" << endl;
	cout << "Accuracy of DBLP" << endl;

	MatrixXd S = rp.allPairs();
	int num = 0;
	int testsetSize = testset.size() / 2;
	vector<vector<int> > allResult;
	vector<int> synNumVec;

	//time
	clock_t total_start = clock();

	for (int i = num * 2; i < 2 * testsetSize; i += 2)
	{
		int synNum = 0;
		string request = testset[i];
		string correct = testset[i + 1];

		vector<int> result = RP_singleQuery(S, vert, request, correct, synNum);
		synNumVec.push_back(synNum);
		allResult.push_back(result);
	}

	metricDP(allResult, testsetSize, synNumVec);

	clock_t total_end = clock();
	float querytime = (total_end - total_start) / (float)CLOCKS_PER_SEC;
	cout << "total querytime = " << querytime << endl;
}

vector<string> testset2dblp(FILE* file)
{
	int length = 100;
	vector<string> testset = vector<string>(length * 2);

	for (int i = 0; i < length; i++)
	{
		string request = "                                                  ";
		string correct = "                                                                                                                                                      ";

		fscanf(file, "%s\t%s\n", &request[0], &correct[0]);

		testset[(2 * i)] = request.erase(request.find_last_not_of(" "));
		testset[(2 * i) + 1] = correct.erase(correct.find_last_not_of(" "));
	}

	return testset;
}

vector<int> singleQuery(algorithm& algo, int vert, string request, string correct, int& synNum) {
	//store in file
	ofstream fout(algo.getResultFile("DP"), ios::app);

	vector<int> rank;
	cout << request << "  begin==================" << endl;

	//get single-sourse
	int id = stoi(request);
	VectorXd res = algo.singleSource(id);

	//sort single-source
	vector<pair<int, float> > row_sim;
	for (int i = 0; i < vert; i++) {
		if (res(i) > 0.000000000001)
			row_sim.push_back(pair<int, float>(i, res(i)));
	}
	if (row_sim.size() == 0) {
		fout << request << "  all-zero " << vert << "\n";
		rank.push_back(vert);
		return rank;
	}
	sort(row_sim.begin(), row_sim.end(), cmp);

	// compare
	int i = 0;
	int pos = 0;

	vector<string> split_correct;
	splitString(correct, split_correct, ",");
	synNum = split_correct.size();
	for (i = 0; i < row_sim.size(); i++) {
		string word = to_string(row_sim[i].first);

		for (int j = 0; j < split_correct.size(); j++) {
			if (word == split_correct[j]) {
				cout << "request = " << request << "  word = " << word << "  i = " << i - 1 << endl;
				fout << request << " " << word << " " << i - 1 << "\n";
				rank.push_back(i - 1);
			}
		}
		if (rank.size() == split_correct.size()) {
			break;
		}
	}
	if (rank.size() == 0) {
		i = row_sim.size() < 10 ? 11 : row_sim.size();
		fout << request << " zero-similarity " << i << "\n";
	}
	fout.close();

	return rank;
}

void allQueries(algorithm& algo, vector<string> testset, int vert) {
	int num = 0;
	int testsetSize = testset.size() / 2;
	vector<vector<int> > allResult;
	vector<int> synNumVec;

	//time
	clock_t total_start = clock();

	for (int i = num * 2; i < 2 * testsetSize; i += 2)
	{
		int synNum = 0;
		string request = testset[i];
		string correct = testset[i + 1];

		vector<int> result = singleQuery(algo, vert, request, correct, synNum);
		synNumVec.push_back(synNum);
		allResult.push_back(result);
	}

	metricDP(allResult, testsetSize, synNumVec);

	clock_t total_end = clock();
	float querytime = (total_end - total_start) / (float)CLOCKS_PER_SEC;
	cout << "total querytime = " << querytime << endl;
}

void mainAllQueries(algorithm& algo, vector<string> testset, int vert) {
	cout << endl << "============ " << algo.getName() << " ============" << endl;
	cout << "Accuracy of DBLP" << endl;

	// do all queries
	allQueries(algo, testset, vert);
}

void acc_DP(double c, double theta, int kmax) {
	cout << "Start loading data" << endl;
	string test_filename = "../Datasets/DP/DP_ground_truth.txt";
	FILE* testsetFile = fopen(test_filename.c_str(), "r");
	vector<string> testset = testset2dblp(testsetFile);
	
	PNGraph G = TNGraph::New();
	SparseMatrix<double> mat;
	vector<Triplet<float> > tripletList;
	string infile = "../Datasets/DP/dblp.txt";
	int vert = 0;
	ReadFile(infile, 0, 0, G, tripletList, vert);
	cout << "Data loading completed" << endl;

	cout << "Start building the matrix" << endl;
	generateMat(G, mat, tripletList, vert);
	cout << "Matrix construction completed" << endl;

	vector<int> M, N;
	CoSimHeat csh = CoSimHeat(G, mat, c, kmax, theta, M, N);

	// CoSimRank
	CoSR cosr = CoSR(mat, c, kmax, false);

	// ExactSim
	ExactSim es = ExactSim(mat, G, c, kmax, false);

	// SimRank*
	SRS srs = SRS(mat, c, kmax, false);

	// ASCOS
	ASCOS ascos = ASCOS(mat, c, kmax, false);

	mainAllQueries(csh, testset, mat.rows());

	// RPCoSim
	/*double delta = 1.2;
	double pf = 1.0 / mat.rows();
	RPCoSim rp = RPCoSim(mat, c, delta, pf);
	RP_DP(rp, testset, mat.rows());*/
}