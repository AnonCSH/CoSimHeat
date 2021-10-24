#pragma once
#include <iostream>
#include "construct.h"
#include "testset.h"
#include "CoSimHeat.h"
#include "CoSR.h"
#include "RPCoSim.h"
#include "ExactSim.h"
#include "SRS.h"
#include "ASCOS.h"
#include "tool.h"
using namespace std;

vector<int> ReadQuery(string queries) {
	//read file
	ifstream fin(queries);
	string str;
	vector<int> q_vec;
	while (getline(fin, str)) {
		stringstream ss;
		ss << str;
		int query;
		ss >> query;
		q_vec.push_back(query);
	}
	fin.close();
	return q_vec;
}

void loadData(PNGraph& G, SparseMatrix<double>& mat, vector<int>& M, vector<int>& N) {
	cout << "Start loading data" << endl;
	string infile = "sandbox/graphs/time/uk-2002/UK.mtx";
	vector<Triplet<float> > tripletList;
	int vert = 0;
	ReadFile(infile, 0, 0, G, tripletList, vert);

	// generate queries
	string M_query = "sandbox/graphs/time/uk-2002/M/M_500.txt";
	string N_query = "sandbox/graphs/time/uk-2002/N/N_50.txt";
	M = ReadQuery(M_query);
	N = ReadQuery(N_query);
	cout << "Data loading completed" << endl;

	cout << "Start building the matrix" << endl;
	generateMat(G, mat, tripletList, vert);
	cout << "Matrix construction completed" << endl;
}

void timeTest(double c, double theta, int kmax) {
	PNGraph G = TNGraph::New();
	vector<int> M, N;
	SparseMatrix<double> mat;

	loadData(G, mat, M, N);

	//start time
	clock_t start, end;
	float duration;
	start = clock();

	// CoSimHeat
	CoSimHeat csh = CoSimHeat(G, mat, c, kmax, theta, M, N);
	csh.partial_pairs();
	//csh.singleSource(1);

	
	// CoSimRank
	//CoSR cosr = CoSR(mat, c, kmax, false);
	//cosr.partial_pairs(M, N);

	// RPCoSim
	/*double delta = 1.2;
	double pf = 1.0 / mat.rows();
	RPCoSim rp = RPCoSim(mat, c, delta, pf);
	rp.allPairs();*/

	// ExactSim
	//ExactSim es = ExactSim(mat, G, c, kmax, false);
	//es.singleSource(N[40]);
	/*for (int i = 0;i < N.size();i++) {
		es.singleSource(N[i]);
	}*/

	// SimRank*
	//SRS srs = SRS(mat, c, kmax, false);
	//srs.singleSource(N[0]);
	/*for (int i = 0;i < N.size();i++) {
		srs.singleSource(N[i]);
	}*/

	// ASCOS
	//ASCOS ascos = ASCOS(mat, c, kmax, false);
	//ascos.singleSource(N[1]);
	/*for (int i = 0;i < N.size();i++) {
		ascos.singleSource(N[i]);
	}*/

	//end time
	end = clock();
	duration = (float)(end - start) / CLOCKS_PER_SEC;
	cout << "Total Time£º" << duration << endl;
	cout << endl;
}

void errorTest(double c, double theta, int kmax) {
	PNGraph G = TNGraph::New();
	vector<int> M, N;
	SparseMatrix<double> mat;

	loadData(G, mat, M, N);

	// Use the result of 200 iterations as the true value
	CoSimHeat csh_200 = CoSimHeat(G, mat, c, 30, theta, M, N);
	CoSimHeat csh_k = CoSimHeat(G, mat, c, kmax, theta, M, N);

	double err = 0.;
	for (int i = 0; i < N.size(); i++) {
		VectorXd truth = csh_200.singleSource(N[i]);
		VectorXd res = csh_k.singleSource(N[i]);
		// The first norm
		err = max((truth - res).lpNorm<1>(), err);
		cout << "iterate :" << i << endl;
	}

	cout << "err = " << err << endl;
}
