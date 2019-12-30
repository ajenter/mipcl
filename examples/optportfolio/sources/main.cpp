#include <iostream>
#include "optportfolio.h"

int main(int argc, const char *argv[])
{
	try {
		Coptportfolio prob(argv[1]);
		prob.model(0.9,1.05);
		prob.optimize();
		prob.printSolution(argv[1]);
	}
	catch(CException* pe) {
		std::cerr << pe->what() << std::endl;
		delete pe;
		return 1;
	}
	return 0;
}

