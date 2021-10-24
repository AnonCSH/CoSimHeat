#pragma once
#ifndef COSIMHEAT_H
#define COSIMHEAT_H
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "construct.h"
#include <fstream>
#include "algorithm.h"
#include "tool.h"
using namespace std;
using namespace Eigen;

class CoSimHeat : public algorithm {
public:
	const string name = "CoSimHeat";

	int kmax;
	double c;
	double theta;
	int edges;
	int vert;
	SparseMatrix<double> mat;
	PNGraph G;
	string resultFileEN;
	string resultFileDP;
	string doutFile;
	VectorXd dout;
	SparseMatrix<double> D, I, W;
	MatrixXd T;
	vector<int> M, N;

	CoSimHeat() {}
	CoSimHeat(SparseMatrix<double>& mat, double decay, int iterations, double theta);
	CoSimHeat(PNGraph& G, SparseMatrix<double>& mat, double decay,
		int iterations, double theta, vector<int> M, vector<int> N);
	~CoSimHeat() {}

	string getName() { return name; }
	string getResultFile(string scene);

	VectorXd singleSource(int id);
	MatrixXd partial_pairs();
	// SYNONYM EXTRACTION 
	VectorXd SYN_SingleSource(int id);

	void getOutDegree(VectorXd& dout, int& edges, string doutFile);
};

CoSimHeat::CoSimHeat(SparseMatrix<double>& M, double decay, int iterations, double theta) {
	mat = M;
	c = decay;
	kmax = iterations;
	this->theta = theta;
	
	resultFileEN = "result/CoSimHeat/ts68.txt";
	doutFile = "../Datasets/EN/outDegree.txt";

	dout.resize(mat.rows());
	SparseMatrix<double> I(mat.rows(), mat.rows());
	for (int i = 0; i < mat.rows(); i++) {
		dout[i] = 0.0;
		I.coeffRef(i, i) = 1.0;
	}
	getOutDegree(dout, edges, doutFile);
	mat = mat - I;
}

CoSimHeat::CoSimHeat(PNGraph& G, SparseMatrix<double>& mat, double decay,
	int iterations, double theta, vector<int> M, vector<int> N) {
	this->mat = mat;
	vert = mat.rows();
	c = decay;
	this->G = G;
	kmax = iterations;
	edges = G->GetEdges();
	this->theta = theta;
	this->M = M;
	this->N = N;
	resultFileDP = "result/CoSimHeat/DP.txt";

	D.resize(vert, vert);
	W.resize(vert, vert);
	I.resize(vert, vert);
	dout.resize(vert);
	for (int i = 0; i < vert; i++) {
		TNGraph::TNodeI NI_j = G->GetNI(i);
		D.coeffRef(i, i) = NI_j.GetOutDeg();
		I.coeffRef(i, i) = 1.0;
		dout[i] = NI_j.GetOutDeg();
	}
}

string CoSimHeat::getResultFile(string scene) {
	if (scene == "EN") {
		return resultFileEN;
	}
	return resultFileDP;
}

void CoSimHeat::getOutDegree(VectorXd& dout, int& edges, string doutFile) {
	//read file
	ifstream fin(doutFile);
	string str;
	while (getline(fin, str)) {
		stringstream ss;
		ss << str;
		int pos, val;
		ss >> pos;
		ss >> val;
		dout[pos - 1] = val;
		edges += val;
	}
	fin.close();
}

VectorXd CoSimHeat::SYN_SingleSource(int id) {
	// CoSimHeat can decompose two parts

	// compute right vector
	VectorXd ej = VectorXd::Zero(mat.rows());
	ej(id) = 1;
	VectorXd u = ej;
	VectorXd v = mat * u;
	float coeff = 1.0;
	for (int k = 1; k <= kmax; k++) {
		coeff = (c / 2) * coeff * (1.0 / k);
		u += coeff * v;
		v = mat * v;
	}

	// T0
	VectorXd tmp = (1 - theta) * u / mat.rows();
	u = dout.transpose() * u;
	u = dout * u;
	u = theta * u / pow(edges, 2) + tmp;

	// compute left vector
	v = (u.transpose() * mat).transpose();
	coeff = 1.0;
	for (int k = 1; k <= kmax; k++) {
		coeff = (c / 2) * coeff * (1.0 / k);
		u += coeff * v;
		v = (v.transpose() * mat).transpose();
	}

	return u;
}

VectorXd CoSimHeat::singleSource(int id) {
	//start time
	clock_t start, end;
	float duration;
	start = clock();

	// W = AT * D.inv() - I
	mat = mat.transpose();
	W = mat * D.cwiseInverse() - I;
	// compute right vector
	VectorXd ej = VectorXd::Zero(vert);
	ej(id) = 1;
	VectorXd u = ej;
	VectorXd v = W.transpose() * u;
	float coeff = 1.0;
	for (int k = 1; k <= kmax; k++) {
		coeff = (c / 2) * coeff * (1.0 / k);
		u += coeff * v;
		v = W.transpose() * v;
	}

	// T0
	VectorXd tmp = (1 - theta) * u / vert;
	u = dout.transpose() * u;
	u = dout * u;
	u = theta * u / pow(edges, 2) + tmp;

	// compute left vector
	v = W * u;
	coeff = 1.0;
	for (int k = 1; k <= kmax; k++) {
		coeff = (c / 2) * coeff * (1.0 / k);
		u += coeff * v;
		v = W * v;
	}

	//end time
	end = clock();
	duration = (float)(end - start) / CLOCKS_PER_SEC;
	cout << "CoSimHeat single-source£º" << duration << endl;
	cout << endl;
	return u;
}

MatrixXd CoSimHeat::partial_pairs() {
	// W = AT * D.inv() - I
	W = (c / 2) * (mat.transpose() * D.cwiseInverse() - I);
	VectorXd z = sqrt(theta) / edges * dout;
	VectorXd v = z;

	for (int k = 1; k <= kmax; k++) {
		v = (1.0 / k) * W * v;
		z = z + v;
	}

	VectorXd z_M(M.size()), z_N(N.size());
	for (int i = 0;i < M.size();i++) {
		z_M[i] = z[M[i]];
	}
	for (int i = 0;i < N.size();i++) {
		z_N[i] = z[N[i]];
	}

	SparseMatrix<double> speye(vert, N.size());
	for (int i = 0;i < N.size();i++) {
		speye.coeffRef(N[i], i) = 1.0;
	}

	SparseMatrix<double> r = (1 - theta) / vert * speye;
	SparseMatrix<double> u = r;

	for (int k = 1; k <= kmax; k++) {
		u = (1.0 / k) * W.transpose() * u;
		r = r + u;
	}

	SparseMatrix<double> p = r, w = r;
	for (int k = 1; k <= kmax; k++) {
		w = (1.0 / k) * W * w;
		p = p + w;
	}
	MatrixXd p_M = MatrixXd::Zero(M.size(), N.size());
	for (int i = 0;i < M.size();i++) {
		p_M.row(i) = p.row(M[i]);
	}

	MatrixXd res = z_M * z_N.transpose() + p_M;

	return res;
}

#endif // !COSIMHEAT_H

