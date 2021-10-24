#pragma once
#ifndef COSR_H
#define COSR_H
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "construct.h"
#include <fstream>
#include <cmath>
#include "algorithm.h"
using namespace std;
using namespace Eigen;

class CoSR : public algorithm{
public:
	const string name = "CoSimRank";

	int kmax;
	double c;
	string resultFileEN;
	string resultFileDP;
	SparseMatrix<double> mat;
	vector<SparseVector<double> > i_pagerank;

	CoSR() {}
	CoSR(SparseMatrix<double>& M, double decay, int iterations, bool isWeighted);
	~CoSR() {}

	string getName() { return name; }
	string getResultFile(string scene);

	//generate transfer matrix
	void getTransferMatrix(SparseMatrix<double>& mat);
	void cal_pagerank(int id);
	double similarity(int i, int j);
	VectorXd singleSource(int id);
	void partial_pairs(vector<int> M, vector<int> N);

	// SYNONYM EXTRACTION
	VectorXd SYN_SingleSource(int id);
};

CoSR::CoSR(SparseMatrix<double>& M,  double decay, int iterations, bool isWeighted) {
	c = decay;
	kmax = iterations;
	resultFileEN = "result/CoSR/ts68.txt";
	resultFileDP = "result/CoSR/DP.txt";
	mat = M;

	if (!isWeighted) {
		getTransferMatrix(mat);
	}
}

string CoSR::getResultFile(string scene) {
	if (scene == "EN") {
		return resultFileEN;
	}
	return resultFileDP;
}

void CoSR::getTransferMatrix(SparseMatrix<double>& mat) {
	for (int k = 0; k < mat.outerSize(); ++k) {
		double sum = mat.col(k).sum();
		for (SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
		{
			it.valueRef() /= sum;
		}
	}
}

void CoSR::cal_pagerank(int id) {
	i_pagerank.clear();
	SparseVector<double> p(mat.rows());
	p.insert(id) = 1;
	i_pagerank.push_back(p);

	for (int k = 1; k < kmax; k++) {
		p = mat * p;
		i_pagerank.push_back(p);
	}
}

double CoSR::similarity(int i, int j) {
	//calculate j's pagerank and sum with i's pagerank
	SparseVector<double> p(mat.rows());
	p.insert(j) = 1;
	double score = 0;
	for (int k = 0; k < kmax; k++) {
		score += pow(c, k) * (i_pagerank[k].dot(p));
		p = mat * p;
	}

	return score;
}

VectorXd CoSR::singleSource(int id) {
	VectorXd res(mat.rows());

	//calculate id's pagerank
	cal_pagerank(id);
	int schedule = 0;
	printf("%d / %d\r", schedule, mat.rows());

	for (int i = 0; i < mat.rows(); i++) {
		double score = similarity(id, i);
		res(i) = score;
		schedule += 1;
		printf("%d / %d\r", schedule, mat.rows());
	}

	return res * (1 - c);
}

VectorXd CoSR::SYN_SingleSource(int id) {
	// construct v0 -- vk
	VectorXd ej = VectorXd::Zero(mat.rows());
	ej(id) = 1;
	vector<VectorXd> v(kmax + 1);
	v[0] = ej;
	for (int k = 1; k <= kmax; k++) {
		v[k] = mat * v[k - 1];
	}

	// compute uk as single-source
	int iter = kmax;
	VectorXd u = v[iter];
	while (iter > 0) {
		u = c * (u.transpose() * mat).transpose() + v[iter - 1];
		iter--;
	}

	return u * (1 - c);
}

void CoSR::partial_pairs(vector<int> M, vector<int> N) {
	int schedule = 0;
	printf("%d / %d\r", schedule, M.size());
	for (int i = 0;i < M.size();i++) {
		//calculate id's pagerank
		cal_pagerank(M[i]);
		for (int j = 0;j < N.size();j++) {
			similarity(M[i], N[j]);
		}
		schedule += 1;
		printf("%d / %d\r", schedule, M.size());
	}
}

#endif // !COSR_H

