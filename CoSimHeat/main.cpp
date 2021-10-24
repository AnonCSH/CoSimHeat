#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <map>
#include "accuracy.h"
#include "acc_dblp.h"
#include "Test.h"
using namespace std;
using namespace Eigen;

int main() {
	double c = 0.6;       // Damping factor
	double theta = 0.1;   // Type weight
	int kmax = 20;        // Number of iterations

	
	#pragma region Search Quality.
		// Search Quality on EN
		acc_EN(c, theta, kmax);
		// Search Quality on DP
		//acc_DP(c, theta, kmax);
	#pragma endregion 

	#pragma region CPU Time.
		//timeTest(c, theta, kmax);
	#pragma endregion 

	#pragma region Iterative Error.
		//errorTest(c, theta, kmax);
	#pragma endregion 

	return 0;
}