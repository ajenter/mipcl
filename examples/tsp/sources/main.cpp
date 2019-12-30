#include <iostream>
#include <except.h>
#include "tsp.h"

int main(int argc, char* argv[])
{
	if (argc < 2) {
		std::cerr << "Enter file name!\n";
		return 1;
	}
	try {
		CTsp gr(argv[1]);
		gr.solve();
		gr.printSolution(argv[1]);
	}
	catch(CException* pe) {
		std::cerr << pe->what() << std::endl;
		delete pe;
		return 2;
	}
  return 0;
}
