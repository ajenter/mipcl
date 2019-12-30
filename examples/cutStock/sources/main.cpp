#include <iostream>
#include <except.h>
#include <cstring>
#include <fstream>

#include "cutstock.h"

void readData(char *fileName, int &m, int &L, int* &l, int* &q)
{
	std::ifstream fin(fileName);
	if (!fin.is_open()) {
		throw new CFileException("readData",fileName);
	}
	fin >> m >> L;
	if (m <= 0 || m > 100000 || L <= 0 || L > 100000) {
		char msg[256];
		sprintf(msg,"Parameters m=%d and L=%d are out of range!",m,L);
		throw new CDataException(msg);
	}
	l = new(std::nothrow) int[m];
	q = new(std::nothrow) int[m];
	if (!l || !q)
		throw new CMemoryException("readData()");

	for (int i=0; i < m; i++) {
		fin >> l[i] >> q[i];
	}
	fin.close();
} // end of readData()

int main(int argc, char** argv)
{
	if (argc < 2) {
        std::cerr << "Enter file name!\n";
        return 1;
	}
	int m,L, *l,*q;
	try {
		readData(argv[1],m,L,l,q);
		CCutStock prob("CutStock");
		prob.init(m,l,q,L);
		prob.optimize();
		char sol[256];
		strcpy(sol,argv[1]);
		strcat(sol,".sol");
		prob.printSolution(sol);
	}
	catch(CException* pe) {
		std::cerr << pe->what() << std::endl;
		delete pe;
		return 1;
	}
	return 0;
} // end of main()
