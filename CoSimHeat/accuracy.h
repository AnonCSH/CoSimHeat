#pragma once
#include <iostream>
#include "construct.h"
#include "testset.h"
#include "CoSimHeat.h"
#include "CoSR.h"
#include "ExactSim.h"
#include "SRS.h"
#include "ASCOS.h"
#include "algorithm.h"

using namespace std;

vector<string> createTestFile(string test_filename, char& from, char& to)
{
	cout << "Testset file: " << test_filename << endl;
	FILE* testsetFile = fopen(test_filename.c_str(), "r");
	vector<string> testset = ConstructTestset(testsetFile, from, to);
	fclose(testsetFile);

	return testset;
}

Construct* createMatrix(string A, string A_dic, PNGraph G, bool transpose)
{
	Construct* aMatrix = new Construct();

	cout << "Link file A: " << A << endl;

	FILE* aFile = fopen(A.c_str(), "r");

	aMatrix->AddEdges(aFile, transpose, G);

	fclose(aFile);

	cout << "Dictionary A file: " << A_dic.c_str() << endl;
	FILE* aDictionaryFile = fopen(A_dic.c_str(), "r");
	aMatrix->ConstructDictionarys(aDictionaryFile);
	fclose(aDictionaryFile);

	return aMatrix;
}

int findId(Construct& Input, string request, string request_type) {
	for (int id = 1; id <= Input.rows; ++id)
	{
		if ((request_type.find(Input.typeDictionary[id]) != string::npos) && Input.dictionary[id] == request)
		{
			return id;
		}
	}

	return 0;
}

vector<int> singleQuery(algorithm& algo, Construct& Input, string request_type, string request, string correct, int& synNum) {
	//store in file
	ofstream fout(algo.getResultFile("EN"), ios::app);

	vector<int> rank;
	// find request id
	int id = findId(Input, request, request_type);
	if (id == 0) {
		cout << "request = " << request << "  not found" << endl;
		fout << request << " NotFound " << -1 << "\n";
		rank.push_back(-1);
		return rank;
	}

	char request_typeC = Input.typeDictionary[id];
	cout << Input.dictionary[id] << "  begin==================" << endl;

	//get single-sourse
	VectorXd res = algo.SYN_SingleSource(id - 1);

	//sort single-source
	vector<pair<int, float> > row_sim;
	for (int i = 0; i < Input.rows; i++) {
		if (res(i) > 1e-13)
			row_sim.push_back(pair<int, float>(i, res(i)));
	}
	if (row_sim.size() == 0) {
		fout << request << "  all-zero " << Input.rows << "\n";
		rank.push_back(Input.rows);
		return rank;
	}
	sort(row_sim.begin(), row_sim.end(), cmp);

	// compare with correct word
	int i = 0;
	int pos = 0;

	vector<string> split_correct;
	splitString(correct, split_correct, ",");
	synNum = split_correct.size();
	int skipNum = 0;
	for (i = 0; i < row_sim.size(); i++) {
		string word = Input.dictionary[row_sim[i].first + 1];
		char type = Input.typeDictionary[row_sim[i].first + 1];
		// if type is wrong, skip
		if (type != request_typeC || word == request) {
			skipNum++;
			continue;
		}
		for (int j = 0; j < split_correct.size(); j++) {
			if (word == split_correct[j]) {
				cout << "request = " << request << "  word = " << word << "  i = " << i - skipNum << endl;
				fout << request << " " << word << " " << i - skipNum << "\n";
				rank.push_back(i - skipNum);
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

void metricEN(vector<vector<int> > allResult, int testsetSize, vector<int> synNumVec) {
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

		int n = 15;
		if (result[0] < n) {
			kendall += ((float)(n - result[0]) / n);
		}

		if (result[0] < n) {
			SPCC += (1 - 6 * (float)pow(result[0], 2) / (n * (pow(n, 2) - 1)));
		}

		float dcg = 0;
		for (int j = 0; j < result.size(); j++) {
			if (result[j] < 100)
				dcg += log(2.0) / log(result[j] + 2.0);
		}
		float idcg = 0;
		for (int j = 0; j < synNum; j++) {
			idcg += log(2.0) / log(j + 2.0);
		}
		if (idcg != 0)
			ndcg += dcg / idcg;

		float AP = 0;
		for (int j = 0; j < synNum; j++) {
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

void allQueries(algorithm& algo, vector<string> testset, Construct& Input) {
	int num = 0;
	int testsetSize = testset.size() / 3;
	vector<vector<int> > allResult;
	vector<int> synNumVec;

	//time
	clock_t total_start = clock();
	
	for (int i = num * 3; i < 3 * testsetSize; i += 3)
	{
		int synNum = 0;
		string request_type = testset[i];
		string request = testset[i + 1];
		string correct = testset[i + 2];

		vector<int> result = singleQuery(algo, Input, request_type, request, correct, synNum);
		synNumVec.push_back(synNum);
		allResult.push_back(result);
	}

	metricEN(allResult, testsetSize, synNumVec);

	clock_t total_end = clock();
	float querytime = (total_end - total_start) / (float)CLOCKS_PER_SEC;
	cout << "total query time = " << querytime << endl;
}

void mainAllQueries(algorithm& algo, vector<string> testset, Construct& Input) {
	cout << endl << "============ " << algo.getName() << " ============" << endl;
	cout << "SYNONYM EXTRACTION" << endl;

	// do all queries
	allQueries(algo, testset, Input);
}

void acc_EN(double c, double theta, int kmax) {
	// generate datasets
	string A = "../Datasets/EN/edges.mtx";
	string A_dic = "../Datasets/EN/dictionary.txt";
	string test_filename = "../Datasets/EN/EN_TS68_ground_truth.txt";

	PNGraph G = TNGraph::New();
	char from, to;
	vector<string> testset;

	cout << "Start loading data" << endl;
	testset = createTestFile(test_filename, from, to);
	Construct* fromMatrix;
	fromMatrix = createMatrix(A, A_dic, G, false);
	cout << "Data loading completed" << endl;

	// CoSimHeat
	CoSimHeat csh = CoSimHeat(fromMatrix->mat, c, kmax, theta);

	// CoSimRank
	//CoSR cosr = CoSR(fromMatrix->mat, c, kmax, true);

	// ExactSim
	//ExactSim es = ExactSim(fromMatrix->mat, G, c, kmax, true);

	// SimRank*
	//SRS srs = SRS(fromMatrix->mat, c, kmax, true);

	// ASCOS++
	//ASCOS ascos = ASCOS(fromMatrix->mat, c, kmax, true);
	
	mainAllQueries(csh, testset, *fromMatrix);
}

