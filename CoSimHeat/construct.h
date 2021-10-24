#pragma once
#ifndef CONSTRUCT_H
#define CONSTRUCT_H
#include <iostream>
#include <vector>
#include <Eigen/Sparse>
#include "mmio.h"
#include "stdafx.h"
#include <fstream>
using namespace std;
using namespace Eigen;

class Construct {
public:
	int rows, columns, edges;
	vector<string> dictionary;
	vector<char> typeDictionary;
	SparseMatrix<double> mat;

	Construct(){}
	~Construct(){}
	void AddEdges(FILE* file, bool transpose, PNGraph G);
	void ConstructDictionarys(FILE* file);
};

void Construct::AddEdges(FILE* file, bool transpose, PNGraph G) {
	//Re-point the position pointer inside the file to the beginning of a stream (data stream/file).
	rewind(file);
	MM_typecode matcode;

	int from, to, rowsThis, columnsThis, nonZerosThis;
	float value;

	if (mm_read_banner(file, &matcode) != 0)
	{
		printf("Could not process Matrix Market banner for A.\n");
		exit(1);
	}
	if (!(mm_is_matrix(matcode) || mm_is_coordinate(matcode) || mm_is_real(matcode) || mm_is_general(matcode)))
	{
		printf("Wrong Matrix Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(1);
	}
	if (mm_read_mtx_crd_size(file, &rowsThis, &columnsThis, &nonZerosThis) != 0)
	{
		printf("Could not read Matrix size\n");
		exit(1);
	}
	if ((columns != 0 && columns != columnsThis) || (rows != 0 && rows != rowsThis))
	{
		printf("Matrix size does not match\n");
		exit(1);
	}

	rows = rowsThis;
	columns = columnsThis;
	if (!transpose)
		mat.resize(columns, rows);
	else
		mat.resize(rows, columns);

	vector<Triplet<float> > tripletList;
	// reading file line by line
	edges = nonZerosThis;
	for (int i = 0; i < nonZerosThis; i++)
	{
		fscanf(file, "%d %d %f\n", &from, &to, &value);
		if (!transpose) {
			tripletList.push_back(Triplet<float>(to - 1, from - 1, value));
		}
		else {
			tripletList.push_back(Triplet<float>(from - 1, to - 1, value));
		}

		if (!G->IsNode(from)) {
			G->AddNode(from);
		}
		if (!G->IsNode(to)) {
			G->AddNode(to);
		}
		G->AddEdge(from, to);
	}
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
}

void Construct::ConstructDictionarys(FILE* file)
{
	dictionary.resize(rows + 1);
	typeDictionary.resize(rows + 1);

	for (int i = 0; i < rows; i++)
	{
		int id;
		char type;
		string word = "                                                ";

		fscanf(file, "%d %s %c\n", &id, &word[0], &type);

		if (id > rows)
		{
			cout << "Dictionary entry not found: " << id << endl;
			continue;
		}

		dictionary[id] = word.erase(word.find_last_not_of(" "));
		typeDictionary[id] = type;
	}
}



#endif // !CONSTRUCT_H

