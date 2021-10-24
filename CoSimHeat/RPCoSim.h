#pragma once
#ifndef RPCOSIM_H
#define RPCOSIM_H
#include <iostream>
#include "stdafx.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <random>
using namespace std;
using namespace Eigen;

class RPCoSim{
public:
	const string name = "RPCoSim";

	double c;
	double delta;
	double eps;
	double pf;
	int vert;
	SparseMatrix<double> mat;

	RPCoSim() {}
	RPCoSim(SparseMatrix<double>& M, double c, double delta, double pf);
	~RPCoSim() {}

	string getName() { return name; }
	// TernarySearch
	double fun(double x);
	double TernarySearch();

	//generate transfer matrix
	void getTransferMatrix(SparseMatrix<double>& P);

	// similarity
	MatrixXd allPairs();
};

RPCoSim::RPCoSim(SparseMatrix<double>& M, double c, double delta, double pf) {
	mat = M;
	this->c = c;
	this->delta = delta;
	this->pf = pf;
	this->vert = mat.rows();
	eps = TernarySearch();
	cout << "eps = " << eps << endl;

	//generate transfer matrix
	getTransferMatrix(mat);
}

void RPCoSim::getTransferMatrix(SparseMatrix<double>& P) {
	for (int k = 0; k < P.outerSize(); ++k) {
		double sum = P.col(k).sum();
		for (SparseMatrix<double>::InnerIterator it(P, k); it; ++it)
		{
			it.valueRef() /= sum;
		}
	}
}

double RPCoSim::fun(double x) {
	double res = log(1 - (c - (1 - c) * delta) / (c * (1 - x))) / ((x - log(1 + x)) * log(c));
	return res;
}

double RPCoSim::TernarySearch() {
	double epsL = (1 - c) / (4 * c) * delta;
	double epsU = (1 - c) / c * delta;
	for (int i = 0; i < 10; i++) {
		double eL = epsL + (epsU - epsL) / 3;
		double eU = epsU - (epsU - epsL) / 3;
		if (fun(eL) < fun(eU))
			epsL = eL;
		else
			epsU = eU;
	}
	return (epsL + epsU) / 2;
}

MatrixXd RPCoSim::allPairs() {
	//start time
	clock_t start, end;
	double duration;
	start = clock();

	int t = ceil(log(1 - (c - (1 - c) * delta) / (c * (1 - eps))) / log(c));
	int d = ceil(2 * log(pow(vert, 2) / (2 * pf)) / (eps - log(1 + eps)));
	cout << "t = " << t << endl;
	cout << "d = " << d << endl;

	static default_random_engine e(time(0));
	static normal_distribution<double> n(0, 1);
	MatrixXd T = MatrixXd::Zero(vert, d).unaryExpr([](double dummy) {return n(e); });

	//Q
	MatrixXd Q = (1.0 / sqrt(d)) * T.transpose() * mat;
	//H  R(d * n)
	MatrixXd H = sqrt(c) * Q;
	//I
	//MatrixXd I = MatrixXd::Identity(vert, vert);
	SparseMatrix<double> I(vert, vert);
	for (int i = 0; i < vert; i++) {
		I.coeffRef(i, i) = 1.0;
	}
	//S
	cout << "S" << endl;
	MatrixXd S = I + H.transpose() * H;

	for (int i = 1; i < t; i++) {
		cout << "i = " << i << endl;
		H = sqrt(c) * H * mat;
		S = S + H.transpose() * H;
	}

	//end time
	end = clock();
	duration = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "RP-CoSim all-pairs£º" << duration << endl;
	cout << endl;

	return S * (1 - c);
}

#endif // !RPCOSIM_H

