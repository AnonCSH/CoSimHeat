#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
using namespace std;

bool cmp(const pair<int, float>& a, const pair<int, float>& b) {
	return a.second > b.second;
}

bool cmp2(const pair<int, int>& a, const pair<int, int>& b) {
	return a.second > b.second;
}

bool cmp1(const pair<string, int>& a, const pair<string, int>& b) {
	return a.second > b.second;
}

void splitString(const string& s, vector<string>& v, const string& c)
{
	string::size_type pos1, pos2;
	pos2 = s.find(c);
	pos1 = 0;
	while (string::npos != pos2)
	{
		v.push_back(s.substr(pos1, pos2 - pos1));
		pos1 = pos2 + c.size();
		pos2 = s.find(c, pos1);
	}
	if (pos1 != s.length())
		v.push_back(s.substr(pos1));
}

void generateMat(PNGraph& G, SparseMatrix<double>& mat, vector<Triplet<float> >& tripletList, int vert) {
	mat.resize(vert, vert);
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
}

void readMtx(string fileName, PNGraph& G, vector<Triplet<float> >& tripletList, int& vert)
{
	ifstream fin(fileName);
	int M, N, L;
	while (fin.peek() == '%')
	{
		fin.ignore(2048, '\n');
	}

	fin >> M >> N >> L;
	vert = M;

	for (int i = 0; i < L; ++i)
	{
		int m, n;
		fin >> m >> n;
		tripletList.push_back(Triplet<float>(m - 1, n - 1, 1));

		if (!G->IsNode(m-1)) {
			G->AddNode(m-1);
		}
		if (!G->IsNode(n-1)) {
			G->AddNode(n-1);
		}
		G->AddEdge(m-1, n-1);

	}
	fin.close();

}

void ReadFile(string infile, int begin, int first, PNGraph& G, vector<Triplet<float> >& tripletList, int& vert) {
	string suffix_str = infile.substr(infile.find_last_of('.') + 1);
	if (suffix_str == "mtx") {
		readMtx(infile, G, tripletList, vert);
		return;
	}

	//read file
	ifstream fin(infile);
	string str;
	while (getline(fin, str)) {
		stringstream ss;
		stringstream toNum;
		ss << str;
		string s;
		int from, to;
		for (int i = 0;i < begin;i++)
			ss >> str;
		ss >> s;
		if (s[0] < '0' || s[0] > '9')
			continue;
		toNum << s;
		toNum >> from;
		ss >> to;
		if (first == 0)
			tripletList.push_back(Triplet<float>(from, to, 1.0));
		else
			tripletList.push_back(Triplet<float>(from - 1, to - 1, 1.0));
		
		if (!G->IsNode(from)) {
			G->AddNode(from);
		}
		if (!G->IsNode(to)) {
			G->AddNode(to);
		}
		G->AddEdge(from, to);
	}

	vert = G->GetNodes();
	fin.close();
}



