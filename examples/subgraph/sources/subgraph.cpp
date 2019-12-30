#include "subgraph.h"
#include <cmath>

#ifndef __ONE_THREAD_
CSubGraph::CSubGraph(const CSubGraph &other, int thread): CMIP(other,thread)
{
	m_iVertNum=other.m_iVertNum;
	m_iEdgeNum=other.m_iEdgeNum;
	m_iK=other.m_iK;
	m_ipHead=other.m_ipHead;
	m_ipTail=other.m_ipTail;
	sprintf(m_sWarningMsg,"cut%d",thread);
	cutInit(m_sWarningMsg);
}

CMIP* CSubGraph::clone(const CMIP *pMip, int thread)
{
	return static_cast<CMIP*>(new CSubGraph(*static_cast<CSubGraph*>(const_cast<CMIP*>(pMip)),thread));
}
#endif

CSubGraph::~CSubGraph()
{
#ifndef __ONE_THREAD_
	if (!m_iThread)
#endif
		if (m_bMemory && m_ipTail)
			delete[] m_ipTail;
	if (cut)
		delete cut;
}

void CSubGraph::readData(const char* fileName)
{
	std::ifstream fin(fileName);
	if (!fin.is_open()) {
		throw CFileException("readData",fileName);
	}
	int m,n;
	fin >> n >> m >> m_iK;
	m_iVertNum=n; m_iEdgeNum=m;
	if (!(m_ipTail= new(std::nothrow) int[3*m+n])) {
		fin.close();
		throw CMemoryException("readData");
	}
	m_ipDeg=(m_ipCost=(m_ipHead=(m_ipTail+m)+m)+m);
	for (int e=0; e < m; e++) {
		fin >> m_ipTail[e] >> m_ipHead[e] >> m_ipCost[e];
	}
	for (int v=0; v < n; v++) {
		fin >> m_ipDeg[v];
	}
	fin.close();
} // end of CSubGraph::readData

void CSubGraph::buildMatrix()
{
	double dpVal[2];
	int n,m,k,ipRow[2];
	cut = 0;
	n=m_iVertNum; m=m_iEdgeNum; k=m_iK;
	openMatrix(3*m,m,5*m);

	for (int v=0; v < n; ++v) {
		addCtr(v,0,-INF,m_ipDeg[v]);
	}
	dpVal[0]=dpVal[1]=1.0;
	for (int e=0; e < m; ++e) {
		ipRow[0]=m_ipTail[e];
		ipRow[1]=m_ipHead[e];
		addColumn(e,VAR_BIN,-m_ipCost[e],0.0,1.0,2,dpVal,ipRow);
	}
	preprocOff();
	closeMatrix();
	cutInit("cut0");
} // end of CSubGraph::buildMatrix

void CSubGraph::cutInit(const char* name)
{
	int n=m_iVertNum, m=m_iEdgeNum;
	int *t=m_ipTail, *h=m_ipHead, ipCol[3];
	double dpVal[3];
	cut = new CLP(name);
	cut->openMatrix(2*m,n+m,3*m);
	cut->addVar(0,0,0.0,1.0,1.0); // variable p(0)
	for (int v=1; v < n; ++v) { // variables p(v), v=1,...,n
		cut->addVar(v,0,0.0,0.0,1.0);
	}
	for (int e=0; e < m; ++e) { // variables gamma(e)
		cut->addVar(n+e,0,0.0,0.0,1.0);
	}
	dpVal[0]=dpVal[2]=1.0; dpVal[1]=-1.0;
	for (int e=0; e < m; ++e) {
		ipCol[0]=t[e]; ipCol[1]=h[e]; ipCol[2]=n+e;
		cut->addRow(2*e,0,0.0,INF,3,dpVal,ipCol);
		ipCol[0]=h[e]; ipCol[1]=t[e];
		cut->addRow(2*e+1,0,0.0,INF,3,dpVal,ipCol);
	}
	cut->preprocOff();
	cut->setScaling(SCL_NO); // switch off scaling
	cut->CLP::switchLpInfoMsg(false);
	cut->closeMatrix();
} // end of CSubGraph::cutInit

bool CSubGraph::separate(int varNum, const double* x, const tagHANDLE* ipColHd, bool bGenFlag)
{
	int n=m_iVertNum, m=m_iEdgeNum;
	bool flag=false;
	double dK=m_iK-0.001;
	for (int e=0; e < m; e++) {
		cut->setObjCoeff(n+e,-x[e]);
	}
	for (int v=1; !flag && v < n; ++v) {
		cut->setVarUpBound(v,0.0);
		cut->optimize();
		if (-cut->CLP::getObjVal() < dK) {
			int sz=0, *hd, *h=m_ipHead, *t=m_ipTail;
			double *p;
			cut->getSolution(p,hd);
			for (int e=0; e < m; ++e) {
				if (fabs(p[t[e]]-p[h[e]]) > 0.5) {
					m_dpArray[sz]=1.0;
					m_ipArray[sz++]=e;
				}
			}
			flag=true;
			if (bGenFlag)
				addCut(-1,CTR_INT,m_iK,INF,sz,m_dpArray,m_ipArray);
		}
		cut->setVarUpBound(v,1.0);
	}
	return flag;
} // end of CSubGraph::separate
