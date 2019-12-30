#include <iostream>
#include <cstring>
#include <except.h>
#include "genAssign.h"


void readData(const char* fileName, int &m, int &n,
		int* &ipMachCap, int* &ipCost, int* &ipProcTime)
{
	int k;
	std::ifstream fin(fileName);
	if (!fin.is_open()) {
		throw new CFileException("readData",fileName);
	}
	fin >> m >> n;
	try {
		ipMachCap=new int[n];
		ipCost=new int[k=m*n];
		ipProcTime=new int[k];
	}
	catch(std::bad_alloc *pe) {
		delete pe;
		throw new CMemoryException("readData");
	}
	for (int i=0; i < n; i++)
		fin >> ipMachCap[i];

	for (int i=0; i < m; i++) {
		for (int j=0; j < n; j++)
			fin >> ipCost[j*m+i];
	}

	for (int i=0; i < m; i++) {
		for (int j=0; j < n; j++)
			fin >> ipProcTime[j*m+i];
	}
	fin.close();
}


int main(int argc, char* argv[])
{
	int m,n, *ipMachCap,*ipCost,*ipProcTime;
	readData(argv[1],m,n,ipMachCap,ipCost,ipProcTime);
	try {
		CGenAssign prob("genAssign",m,n,ipMachCap,ipCost,ipProcTime);
		prob.optimize();
		prob.printSolution(argv[1]);
	}
	catch(CException* pe) {
		std::cerr << pe->what() << std::endl;
		delete pe;
		return 1;
	}
	delete[] ipMachCap;
	delete[] ipCost;
	delete[] ipProcTime;
	return 0;
}
