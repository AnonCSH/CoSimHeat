#pragma once
#include <vector>
#include <iostream>
using namespace std;

class algorithm {
public:
	virtual string getName() = 0;
	virtual string getResultFile(string scene) = 0;
	virtual VectorXd SYN_SingleSource(int id) = 0;
	virtual VectorXd singleSource(int id) = 0;
};
