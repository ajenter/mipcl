#include <iostream>
#include <cstring>
#include "fcnf.h"

int main(int argc, char *argv[])
{
	if (argc < 2) {
		std::cerr << "File name is omitted!\n";
		return 1;
	}
	try {
		CFCNF prob("FCFN");
		prob.readNet(argv[1]);
		prob.model();
		prob.optimize();

		char name[128];
		strcpy(name,argv[1]);
		strcat(name,".sol");
		prob.printSolution(name);
	}
	catch(CException* pe) {
		std::cerr << pe->what() << std::endl;
		delete pe;
		return 1;
	}
	return 0;
}
