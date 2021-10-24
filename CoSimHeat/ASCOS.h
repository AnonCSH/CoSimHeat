#pragma once
#ifndef ASCOS_H
#define ASCOS_H
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "construct.h"
#include <fstream>
#include <cmath>
#include "algorithm.h"
using namespace std;
using namespace Eigen;

class ASCOS : public algorithm{
public:
	const string name = "ASCOS++";

	int kmax;
	double c;
	string resultFileEN;
	string resultFileDP;
	SparseMatrix<double> mat;

	ASCOS() {}
	ASCOS(SparseMatrix<double>& M, double decay, int iterations, bool isWeighted);
	~ASCOS() {}

	string getName() { return name; }
	string getResultFile(string scene);

	//generate transfer matrix
	void getTransferMatrix(SparseMatrix<double>& mat);
	VectorXd singleSource(int id);
	// SYNONYM EXTRACTION naive single-source
	VectorXd SYN_SingleSource(int id);

};

ASCOS::ASCOS(SparseMatrix<double>& M, double decay, int iterations, bool isWeighted) {
	mat = M;
	c = decay;
	kmax = iterations;
	resultFileEN = "result/ASCOS/ts68.txt";
	resultFileDP = "result/ASCOS/DP.txt";

	if (!isWeighted) {
		getTransferMatrix(mat);
	}
}

string ASCOS::getResultFile(string scene) {
	if (scene == "EN") {
		return resultFileEN;
	}
	return resultFileDP;
}

void ASCOS::getTransferMatrix(SparseMatrix<double>& mat) {
	for (int k = 0; k < mat.outerSize(); ++k) {
		double sum = mat.col(k).sum();
		for (SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
		{
			it.valueRef() /= sum;
		}
	}
}

VectorXd ASCOS::singleSource(int id) {
	//start time
	clock_t start, end;
	double duration;
	start = clock();

	VectorXd res = SYN_SingleSource(id);

	//end time
	end = clock();
	duration = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "ASCOS++ singleSource£º" << duration << endl;
	cout << endl;

	return res;
}

VectorXd ASCOS::SYN_SingleSource(int id) {
	// construct v0 -- vk

	VectorXd ej = VectorXd::Zero(mat.rows());
	ej(id) = 1;

	// compute uk as single-source
	int iter = kmax;
	VectorXd u = ej;
	while (iter > 0) {
		u = c * (u.transpose() * mat).transpose() + ej;
		iter--;
	}

	VectorXd res = u * (1 - c);

	return res;
}

#endif // !ASCOS_H

