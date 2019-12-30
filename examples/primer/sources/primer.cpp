#include <cmip.h>
#include <except.h>
#include <iostream>

int variant1()
{
	try {
		CMIP prob("MIPCLtest"); // 1
		prob.openMatrix(2,2,4); // 2
		prob.addVar(0,CMIP::VAR_INT, 100.0,0.0, CLP::VAR_INF); // 3
		prob.addVar(1,CMIP::VAR_INT,  64.0,0.0, CLP::VAR_INF); // 4
		prob.addCtr(0,0,-CLP::INF,250); // 5
		prob.addCtr(1,0,-CLP::INF,4); // 6
		prob.addEntry(50.0,0,0); // 7.1
		prob.addEntry(31.0,0,1); // 7.2
		prob.addEntry(-3.0,1,0); // 7.3
		prob.addEntry( 2.0,1,1); // 7.4
		prob.closeMatrix(); // (8)
		prob.optimize(); // (9)
		prob.printSolution("primer.sol"); // (10)
	}
	catch(CException* pe) {
		std::cerr << pe->what() << std::endl;
		delete pe;
		return 1;
	}
	return 0;
} // end of variant1
/////////////////////////////////////////////////////////////////////////////

// Constants that represent our test IP.
const int n=2, m=2, nz=4;
double c[] = {100,64}; // cost vector

double A[][n] = { // matrix
		{50,31},
		{-3,2}
};
double b[] = {250,4}; // right-hand side vector

int ind[] = {0,1}; // array of indices

int variant2()
{
	try {
		CMIP prob("MIPCLtest"); // create new MIP problem
		prob.openMatrix(n,m,nz); // open matrix

		for (int j=0; j < n; ++j) // adds n variables:
			prob.addVar(j,CMIP::VAR_INT, c[j],0.0, CLP::VAR_INF);

		for (int i=0; i < m; ++i) // add m rows (constraints):
			prob.addRow(i,0,-CLP::INF,b[i],n,A[i],ind);

		prob.closeMatrix(); // close matrix
		prob.optimize(); // solve problem

		prob.printSolution("primer.sol"); // print solution
	}
	catch(CException* pe) {
		std::cerr << pe->what() << std::endl;
		delete pe;
		return 1;
	}
	return 0;
}  // end of variant2


int main(int argc, char *argv[])
{
	int v=1;
	if (argc > 1)
		if (argv[1][0] == '2')
			v=2;
	return (v == 1)? variant1(): variant2();
}
