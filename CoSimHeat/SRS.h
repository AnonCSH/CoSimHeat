#pragma once
#ifndef SRS_H
#define SRS_H
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "construct.h"
#include <fstream>
#include <cmath>
#include "algorithm.h"
using namespace std;
using namespace Eigen;

class SRS : public algorithm{
public:
	const string name = "SimRank*";

	int kmax;
	double c;
	string resultFileEN;
	string resultFileDP;
	SparseMatrix<double> mat;

	SRS() {}
	SRS(SparseMatrix<double>& M, double decay, int iterations, bool isWeighted);
	~SRS() {}

	string getName() { return name; }
	string getResultFile(string scene);

	//generate transfer matrix
	void getTransferMatrix(SparseMatrix<double>& mat);
	VectorXd singleSource(int id);
	// SYNONYM EXTRACTION naive single-source
	VectorXd SYN_SingleSource(int id);

};

SRS::SRS(SparseMatrix<double>& M, double decay, int iterations, bool isWeighted) {
	mat = M;
	c = decay;
	kmax = iterations;
	resultFileEN = "result/SRS/ts68.txt";
	resultFileDP = "result/SRS/DP.txt";

	if (!isWeighted) {
		getTransferMatrix(mat);
	}
}

string SRS::getResultFile(string scene) {
	if (scene == "EN") {
		return resultFileEN;
	}
	return resultFileDP;
}

void SRS::getTransferMatrix(SparseMatrix<double>& mat) {
	for (int k = 0; k < mat.outerSize(); ++k) {
		double sum = mat.col(k).sum();
		for (SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
		{
			it.valueRef() /= sum;
		}
	}
}

VectorXd SRS::SYN_SingleSource(int id) {
	//u[0] ... u[k]
	VectorXd ej = VectorXd::Zero(mat.rows());
	ej(id) = 1;
	vector<VectorXd> u(kmax + 1);
	u[0] = ej;
	for (int i = 1; i <= kmax; i++) {
		u[i] = (c / 2.0) * u[i - 1];
	}

	for (int a = 0; a < kmax; a++) {
		int iter = kmax;
		for (int i = a; i < kmax; i++) {
			u[iter - 1] = u[iter - 1] + mat * u[iter];
			iter--;
		}
	}

	//v[0]...v[k - 1]
	vector<VectorXd> v(kmax);
	int iter = kmax;
	v[0] = u[iter--];
	for (int i = 1; i < kmax; i++) {
		v[i] = mat.transpose() * v[i - 1] + u[iter--];
	}

	//finally
	VectorXd res = (1 - c) * ((v[kmax - 1].transpose() * mat).transpose() + ej);

	return res;
}

VectorXd SRS::singleSource(int id) {
	//start time
	clock_t start, end;
	double duration;
	start = clock();

	VectorXd res = SYN_SingleSource(id);

	//end time
	end = clock();
	duration = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "SimRank* singleSource£º" << duration << endl;
	cout << endl;

	return res;
}

#endif // !SRS_H

