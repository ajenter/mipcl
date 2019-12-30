#include <iostream>
#include "fcnf.h"

void CFCNF::readNet(const char* fileName)
{
	std::ifstream fin(fileName);
	if (!fin.is_open()) {
		throw CFileException("readNet",fileName);
	}

	int n,m;
	fin >> n >> m;
	m_iVertNum=n; m_iEdgeNum=m;

	if (!(m_ipTail = new(std::nothrow) int[5*m+n])) {
		fin.close();
		throw CMemoryException("readNet: ");
	}
	m_bMemory=true;

	m_ipDem=(m_ipCost=(m_ipFxCost=(m_ipCap=(m_ipHead=m_ipTail+m)+m)+m)+m)+m;

	for (int i=0; i < n; ++i) {
		fin >> m_ipDem[i];
	}

	for (int i=0; i < m; ++i) {
		fin >> m_ipTail[i] >> m_ipHead[i] >> m_ipCap[i]
			>> m_ipFxCost[i] >> m_ipCost[i];
	}

	fin.close();
} // end of CFCNF::readNet

void CFCNF::model()
{
	double w,dpVal[3];
	int ipRow[3],n=m_iVertNum, m=m_iEdgeNum;

	openMatrix(m+n,2*m,4*m);
	setObjSense(false);

	for (int i=0; i < n; ++i) {
		addCtr(i,0,m_ipDem[i],m_ipDem[i]);
	}
	for (int i=0; i < m; ++i) {
		addCtr(n+i,0,-INF,0.0);
	}
	dpVal[0]=-1.0; dpVal[1]=dpVal[2]=1.0;
	for (int r,i=0; i < m; ++i) {
		ipRow[0]=m_ipTail[i]; ipRow[1]=m_ipHead[i]; ipRow[2]=r=n+i;
		addColumn(i,0,m_ipCost[i],0.0,VAR_INF,3,dpVal,ipRow);
		w=-m_ipCap[i];
		addColumn(i+m,VAR_BIN,m_ipFxCost[i],0.0,1.0,1,&w,&r);
	}
	closeMatrix();
} // end of CFCNF::buildProb

void CFCNF::printSolution(const char* name)
{
	double *dpX;
	int m,m2,*ipHd;
	m=m_iEdgeNum; m2=2*m;
	std::ofstream fout(name);
	if (isSolution()) {
		getSolution(dpX,ipHd);
		fout << "Nonzero flows:\n";
		for (int i=0; i < m2; ++i) {
			if (ipHd[i] < m) {
				if (dpX[i] > 0.5) {
					fout << "f(" << m_ipTail[ipHd[i]] << "," 
						<< m_ipHead[ipHd[i]] << ")="
						<< dpX[i] << std::endl;
				}
			}
		}
	}
	else fout << "Problem has no solution!\n";
	fout.close();
} // end of CFCNF::printSolution
