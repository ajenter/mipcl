namespace KNAPSACK {
	double intKnapsack(int n, double *c, int *a, int b, int*x, double *dpMem=0);
// solves integer knapsack 
// max{sum_{j=1}^n c_j x_j : sum_{j=1}^n a_j x_j \le b, x_j\in Z_+}

	void getMemForBinKnapsack(int n, int b, double* &dpF, int* &ipBt);
// allocates memory for solving 0,1-knapsack problems
// \sum_{j=1}^n a_j x_j \le b
// for n \le m_iTskNum and  b \le q
// frees memory previously allocated by GetMemForBinKnapsack
	double binKnapsack(int n, double *c, int* a, int b, int *x, double *dpF, int *ipBt);
// solves 0,1-knapsack problem
// max{sum_{j=1}^n c_j x_j : sum_{j=1}^n a_j x_j \le b, x_j\in {0,1}}
// It is assumed that a_j \le b for all j
};
