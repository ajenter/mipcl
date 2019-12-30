#include <iostream>
#include "subgraph.h"

int main(int argc, char* argv[])
{
	try {
		CSubGraph prob("subgraph");
		prob.readData(argv[1]);
		prob.buildMatrix();
		prob.optimize();
		prob.printSolution("test.sol");
	}
	catch(CException* pe) {
		std::cerr << pe->what() << std::endl;
		delete pe;
		return 1;
	}
	return 0;
}
