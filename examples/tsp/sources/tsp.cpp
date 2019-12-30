//////////////////////////////////////////////////////// 
// Ctsp.cpp: implementation of the Ctsp class //
////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <cstring>
#include <cstdio>
#include <except.h>
#include "tspPool.h"
#include "tsp.h"

#define BAD_CUT 1.999
#define VERY_BAD_CUT 1.5

#define MAX_DEGREE 20
#define INIT_DEGREE 10

#ifndef MAX_STR_LEN
#define MAX_STR_LEN 255
#endif

#define MAX_ROW_NUM (m_iPointNum+10*m_MaxEdgeNum)

// Construction/Destruction

CTsp::CTsp(const char *name): CMIP("tsp")
{
	int k;
	m_pTspPool=0;
	m_pNet=0;
	m_dpCoordX=0;
	m_ipNextOnTour=0;
	readPoints(name);
	
	k=m_iPointNum*MAX_DEGREE;
	m_pNet= new CFlowNet<double>(m_iPointNum+k,2*k,CFlowNet<double>::m_iUpCapMsk);
	m_pNet->COMB_AND_CUT_getMem();
	m_pTspPool = new CTspPool(m_iPointNum);
	setMIP();
} // end of CTsp::CTsp

#ifndef __ONE_THREAD_
CTsp::CTsp(const CTsp &other, int thread): CMIP(other,thread)
{
	m_enType=other.m_enType;
	m_iPointNum=other.m_iPointNum;
	m_dpCoordX=other.m_dpCoordX;
	m_dpCoordY=other.m_dpCoordY;
	m_dMaxDist=other.m_dMaxDist;
	m_ipNextOnTour=other.m_ipNextOnTour;
	m_pTspPool=other.m_pTspPool;
	int u=m_iPointNum*MAX_DEGREE;
	m_pNet= new CFlowNet<double>(m_iPointNum+u,2*u,CFlowNet<double>::m_iUpCapMsk);
	m_pNet->COMB_AND_CUT_getMem();
} // end of CTsp::CTsp(const CTsp &other, int thread)

CMIP* CTsp::clone(const CMIP *pMip, int thread)
{
	return static_cast<CMIP*>(new CTsp(*static_cast<CTsp*>(const_cast<CMIP*>(pMip)),thread));
}
#endif

CTsp::~CTsp()
{
#ifndef __ONE_THREAD_
	if (!m_iThread) {
#endif
	if (m_dpCoordX) 
		delete[] m_dpCoordX;
	if (m_ipNextOnTour)
		delete[] m_ipNextOnTour;
	if (m_pTspPool)
		delete m_pTspPool;
#ifndef __ONE_THREAD_
	}
#endif
	if (m_pNet)
		delete m_pNet;
} // end of CTsp::~CTsp()

double CTsp::dist(int i, int j)
{
	double d,dx,dy;
	dx=m_dpCoordX[i]-m_dpCoordX[j];
	dy=m_dpCoordY[i]-m_dpCoordY[j];
	d=sqrt(dx*dx+dy*dy);
	if (m_enType==ATT) {
		int l;
		if ((l=static_cast<int>(floor(d))) < d+1.0e-8)
			++l;
		d=l;
	}
	return d;
} // end of CTsp::dist

double CTsp::goToNearest(int s, int* tau)
{
	int v,w,n=m_iPointNum;
	double d,q, length=0;
	bool* bpVisited=(bool*)(m_pTspPool->getBufPtr());
	memset(bpVisited,0,n*sizeof(char));
	bpVisited[v=s]=true;
	for (int i=1; i < n; ++i) {
		d=std::numeric_limits<double>::max();
		for (int j=0; j < n; ++j)
			if (!bpVisited[j])
			if (d > (q=dist(v,j))) {
				w=j;
				d=q;
			}
		tau[v]=w;
		bpVisited[v=w]=1;
		length+=d;	
	}
	tau[v]=s;
	length+=dist(v,s);
	return length;
}  // end of CTsp::goToNearest

void CTsp::approximate()
{
	double Length;
	char str[64];
	int* tau=(m_pTspPool->getBufPtr())+m_iPointNum;
	m_dTourLength=std::numeric_limits<double>::max();
	for (int s=0; s < m_iPointNum; ++s) {
		Length=goToNearest(s,tau);
//		if (m_iPointNum > 4) 
//			Hc.KL(Length,tau);
		if (Length < m_dTourLength) {
			m_dTourLength=Length;
			memcpy(m_ipNextOnTour,tau,m_iPointNum*sizeof(int));
		}
	}
	sprintf(str,"%s%lf","Go-to-Nearest Tour Length = ",m_dTourLength);
	infoMessage(str);
}  // end of CTsp::approximate()

void CTsp::allocMemForCoords()
{
	m_dpCoordX = new double[2*m_iPointNum];
	m_ipNextOnTour= new int[m_iPointNum];
	if (!m_dpCoordX || !m_ipNextOnTour) {
		throw new CMemoryException("CTsp::allocMemForCoords");
	}
	m_dpCoordY=m_dpCoordX+m_iPointNum;
}

void CTsp::readPoints(const char* fileName)
{
	char *ch, st[MAX_STR_LEN+1];
	std::ifstream in_stream(fileName);
	if (!in_stream) {
		throw new CFileException("CTsp::readPoints",fileName);
	}
	in_stream.getline(st,MAX_STR_LEN);
	for (ch=st; *ch != ':'; ch++);
	for (++ch; *ch == ' '; ch++);
	setProblemName(ch);
	in_stream.getline(st,MAX_STR_LEN);
	in_stream.getline(st,MAX_STR_LEN);
	in_stream.getline(st,MAX_STR_LEN);
	for (ch=st; *ch != ':'; ch++);
	sscanf(++ch,"%d",&m_iPointNum);
	in_stream.getline(st,MAX_STR_LEN);
	for (ch=st; *ch != ':'; ch++);
	for (++ch; *ch == ' '; ch++);
	if (!strncmp(ch,"GEO",3))
		m_enType=GEO;
	else if (!strncmp(ch,"ATT",3))
		m_enType=ATT;
	else if (!strncmp(ch,"EUC_2D",6))
		m_enType=EUC_2D;
	else m_enType=ATT;

	while (strncmp(st,"NODE_COORD_SECTION",18))
		in_stream.getline(st,MAX_STR_LEN);
	allocMemForCoords();

	int Nom;
	double d, maxDist;
	for (int v=0; v < m_iPointNum; ++v) {
		in_stream >> Nom;
		in_stream >> m_dpCoordX[v];
		in_stream >> m_dpCoordY[v];
//		in_stream.getline(st,MAX_STR_LEN-1);
	}
	in_stream.close();
	maxDist=-1.0;
	for (int w=1; w < m_iPointNum; ++w) {
		for (int v=0; v < w; ++v)
			if (maxDist < (d=dist(v,w)))
				maxDist=d;
	}
	m_dMaxDist=maxDist+1.0;
} // end of CTsp::readPoints()

/////////////////////////////////////////////////////
//                     MIP part                    //
/////////////////////////////////////////////////////
void CTsp::setMIP()
{
	approximate();
	buildActiveGraph();
	double dpVal[2];
	int n=m_pNet->getVertNum(), m=m_pNet->getEdgeNum();
	int hd,ipRow[2];
	openMatrix(n*3,m*2,m*10,true,true);
	for (int i=0; i < n; ++i) {
		addCtr(i | 0xE0000000,0,-INF,2.0);
	}
	dpVal[0]=dpVal[1]=1.0;
	for (int hd,j=0; j < m; ++j) {
		hd=(ipRow[0]=m_pNet->getTail(j)) << 16;
		hd|=(ipRow[1]=m_pNet->getHead(j));
		addColumn(hd,VAR_BIN,m_dMaxDist-dist(ipRow[0],ipRow[1]),0.0,1.0,2,dpVal,ipRow);
	}
	preprocOff();
	setScaling(SCL_NO);
	closeMatrix();
} // end of CTsp::setMIP

void CTsp::buildActiveGraph()
{
	int n,u,w,*ipDeg, *ipFlag;
	n=m_iPointNum;
	if (!(ipDeg = new int[2*n])) {
		throw new CMemoryException("CTsp::BuildActiveGraph");
	}
	ipFlag=ipDeg+n;
	m_pNet->setVertNum(n);
	m_pNet->setEdgeNum(0);
	for (int v=0; v < n; ++v) {
		ipDeg[v]=2;
		m_pNet->addEdge(v,w=m_ipNextOnTour[v]);
	}
	int e=m_iPointNum;
	for (int v=0; v < n; ++v) {
		memset(ipFlag,0,m_iPointNum*sizeof(int));
		ipFlag[v]=1;
		for (int j=0; j < e; ++j) {
			w=m_pNet->getTail(j);
			u=m_pNet->getHead(j);
			if (w == v)
				ipFlag[u]=1;
			else if (u == v)
				ipFlag[w]=1;
		}
		for (; ipDeg[v] < INIT_DEGREE; ++e) {
			double q,d=std::numeric_limits<double>::max();
			for (int i=0; i < n; ++i) {
				if (!ipFlag[i])
					if (d > (q=dist(v,i))) {
						d=q; w=i;
					}
			}
			ipFlag[w]=1;
			++ipDeg[v]; ++ipDeg[w];
			m_pNet->addEdge(v,w);
		}
	}
	m_pNet->buildEdgeList();
	delete[] ipDeg;
}  // end of CTsp::buildActiveGraph

bool CTsp::getRow(tagHANDLE hd, int n, const tagHANDLE* ipColHd,
                        int& type, double& b1, double& b2,
                        int& sz, double* dpVal, int* ipCol, bool &bScaled)
{
	type=CTR_INT;
	b1=-INF;
	sz=m_pTspPool->buildRow(hd,n,ipColHd,dpVal,ipCol,b2);
	m_pTspPool->lockCtr(hd);
	bScaled=true;
	return true;
} // end of CTsp::getRow()

bool CTsp::getColumn(tagHANDLE hd, int m, const tagHANDLE* ipRowHd,
                int& type, double& cost, double& d1, double& d2,
                int& sz, double* dpVal, int* ipRow)
{
	int v,w;
	v=hd & (0x0000FFFF);
	w=hd >> 16;
	d1=0.0; d2=1.0;
	type=VAR_BIN;
	cost=m_dMaxDist-dist(v,w);

	dpVal[0]=dpVal[1]=1.0;
	ipRow[0]=v; ipRow[1]=w;
	sz=2;
//	sz=m_pTspPool->buildColumn(v,w,m,ipRowHd,dpVal,ipRow);
	return true;
} // end of CTsp::getColumn()

void CTsp::lockCtr(tagHANDLE hd)
	{m_pTspPool->lockCtr(hd);}

void CTsp::unlockCtr(tagHANDLE hd)
	{m_pTspPool->unlockCtr(hd);}

/////// S E P A R A T I O N ////////
bool CTsp::separate(int n, const double* dpX, const tagHANDLE* ipColHd, bool bGenFlag)
{
	int poolCuts=0, cutCuts=0;
	if (!(poolCuts=separateFromPool(n,dpX,ipColHd,bGenFlag))) {
		buildSupportGraph(n,dpX,ipColHd);
		cutCuts=cutSeparate(n,dpX,ipColHd,bGenFlag);
	}
	return (poolCuts+cutCuts > 0)? true: false;
} // end of CTsp::separate()

int CTsp::separateFromPool(int n, const double* dpX, const int* ipColHd,
							bool bGenFlag)
{
	double rhs;
	int sz, m, cutNum, thread;
	cutNum=0;
	thread=m_iThread;
	m=m_iM;
	m_pTspPool->wrLockPool();
	for (int i=m_iM0; i < m; ++i) {
		if (m_ipRowHd[i] >= 0)
			m_pTspPool->markCtr(m_ipRowHd[i],thread);
	}
	m_pTspPool->wrUnlockPool();
	for (int hd=0; sz= m_pTspPool->getNextCut(thread,hd,n,dpX,ipColHd,m_ipArray,m_dpArray,rhs); ++hd) {
		++cutNum;
		if (!bGenFlag)
			break;
		safeAddCut(hd,CTR_INT,-INF,rhs,sz,m_dpArray,m_ipArray);
		m_pTspPool->lockCtr(hd);
	}
	m_pTspPool->wrLockPool();
	for (int i=m_iM0; i < m; ++i) {
		if (m_ipRowHd[i] >= 0)
			m_pTspPool->unmarkCtr(m_ipRowHd[i],thread);
	}
	m_pTspPool->wrUnlockPool();
	return cutNum;
} // end of CTsp::separateFromPool

void CTsp::buildSupportGraph(int n, const double* dpX, const tagHANDLE* ipColHd)
{
	double tol=getIntTol();
	int hd;
	m_pNet->reset(m_iPointNum);
	for (int e=0; e < n; ++e)
		if (dpX[e] > tol) {
			hd=(int)ipColHd[e];
			m_pNet->addEdge(hd >> 16,hd & 0x0000FFFF,dpX[e]);
		}
	m_pNet->buildEdgeList();
} // end of CTsp::buildSupportGraph()

int CTsp::cutSeparate(int varNum, const double* dpX, const tagHANDLE* ipColHd, bool genFlag)
{
	m_pNet->MC_MinCut(VERY_BAD_CUT);
	int* ipCut=m_pNet->MC_getCut();

	if (m_pNet->MC_getCutValue() > BAD_CUT)
		return 0;

	if (genFlag) {
		double w1,w2;
		int hd,k,n,b;
		n=m_iPointNum;
		for (int v=k=0; v < n; ++v) {
			if (ipCut[v])
				++k;
		}
		w1=w2=-0.001;
		for (int e=0; e < varNum; ++e) {
			hd=static_cast<int>(ipColHd[e]);
			if (ipCut[hd >> 16] == ipCut[hd & 0x0000FFFF]) {
				if (ipCut[hd & 0x0000FFFF] == 1)
					w1+=dpX[e];
				else
					w2+=dpX[e];
			}
		}
		if (k <= n-k) {
			if (w1 > (b=k-1))
				k=1;
			else k=(w2 > (b=n-k-1))? 0: -1;
		}
		else k=(w2 > (b=n-k-1))? 0: -1;
		if (k >= 0) {
			if (!k) {
				for(int v=0; v < n; ++v)
					ipCut[v]=1-ipCut[v];
			}
			for (int e=k=0; e < varNum; ++e) {
				hd=ipColHd[e];
				if  (n=CTspPool::getCoefficient(ipCut[hd >> 16],ipCut[hd & 0x0000FFFF])) {
					m_dpArray[k]=n;
					m_ipArray[k++]=e;
				}
			}
			hd=m_pTspPool->addCut(b,ipCut);
			safeAddCut(hd,CTR_INT,-INF,b,k,m_dpArray,m_ipArray);
		}
	}
	return 1;
} // end of CTsp::cutSeparate()

// It is assumed that the support graph has been built in `sepatate()`
bool CTsp::genCut1(int n, const double* dpX, const tagHANDLE* ipColHd)
{
	int blossomCuts=BLOSSOM_separate(n,dpX,ipColHd);
	return (blossomCuts > 0)? true: false;
} // end of CTsp::separate()

int CTsp::BLOSSOM_separate(int n, const double* dpX, const tagHANDLE* ipColHd)
{
	int *ipComb;
	int k;
	if (k=m_pNet->maxBlossom()) {
		int hdSize, thNum;
		for (int i=0; i < k; ++i) {
			m_pNet->getMaxBlossom(i,hdSize,thNum,ipComb);
			addCombCut(hdSize + thNum/2,ipComb, n,dpX,ipColHd);
		}
	}
	else { // k=0
		int hdSize, thNum;
		m_pNet->blossom();
		for (int s=1; thNum=m_pNet->getNextBlossom(s,0.001,hdSize,ipComb);) {
			addCombCut(hdSize + thNum/2,ipComb, n,dpX,ipColHd);
			++k;
		}
	}
	return k;
} // end of CTsp::BLOSSOM_separate()
//	if (!m_iNode)
//		cutInfo(MISC::getTimeInMlSec()-m_lStartTime,
//				m_iAutoCutRound-1,CLP::getObjVal(),m_iFracNum,cutNum);
//

void CTsp::addCombCut(int b, int* ipComb,
			int n, const double* dpX, const tagHANDLE* ipColHd)
{
	double b0;
	int hd=m_pTspPool->addCut(b,ipComb);
	int sz=m_pTspPool->buildRow(hd,n,ipColHd,m_dpArray,m_ipArray,b0);
	safeAddCut(hd,0,-INF,b0,sz,m_dpArray,m_ipArray);
} // end of CTsp::addCombCut()

bool CTsp::generateColumns(int m, const tagHANDLE* ipRowHd, const double* dpY)
{
	double redCost,cost, bestRedCost,cq;
	int q,sz,num=0;
	buildSupportGraph();
	for (int i=1; i < m_iPointNum; ++i) {
		bestRedCost=0.01; // GetRedCostTol();
		q=-1;
		for (int j=0; j < i; ++j) {
			if (m_pNet->getEdgeNo(i,j) < 0) {
				sz=m_pTspPool->buildColumn(i,j,m,ipRowHd,m_dpArray,m_ipArray);
				redCost=(cost = m_dMaxDist-dist(i,j));
				for (int k=0; k < sz; ++k) 
					redCost-=m_dpArray[k]*dpY[m_ipArray[k]];
				if (redCost > bestRedCost) {
					q=j;
					bestRedCost=redCost;
					cq=cost;
				}	// if (RedCost > BestRedCost)
			}	// if (m_pNet->GetEdgeNo(i,j) == NIL)
		}  // for (int j=i+1; ...
		if (q >= 0) {
			sz=m_pTspPool->buildColumn(i,q,m,ipRowHd,m_dpArray,m_ipArray);
			addNewColumn((i<<16) | q,VAR_BIN,cq,0.0,1.0,
                    sz,m_dpArray,m_ipArray,false,false,0,true);
			++num;
		}	// if (BestRedCost > TolDual)
	}  // for (int i=0;
	return (num)? true: false;
}  // end of CTsp::generateColumns()

void CTsp::buildSupportGraph()
{
	m_pNet->reset(m_iPointNum);
	int hd,n=getVarNum();
	for (int i=0; i < n; ++i) {
		hd=getVarHandle(i);
		m_pNet->addEdge(hd >> 16,hd & 0x0000FFFF);
	}
	m_pNet->buildEdgeList();
}	// end of CTsp::buildSupportGraph()

void CTsp::changeRecord(double dObjVal,
		int n, const double* dpX, const tagHANDLE* ipHd)
{
	int v,w;
	int* ipFirst=m_ipArray;
	int* ipSecond=reinterpret_cast<int*>(m_dpArray);
	for (int i=0; i < m_iPointNum; ++i)
		ipFirst[i]=ipSecond[i]=-1;
	for (int i=0; i < n; ++i) {
		if (dpX[i] > 0.5) {
			v=ipHd[i] >> 16;
			w=ipHd[i] & (0x0000FFFF);
			if (ipFirst[v] >= 0)
				ipSecond[v]=w;
			else ipFirst[v]=w;
			if (ipFirst[w] >= 0)
				ipSecond[w]=v;
			else ipFirst[w]=v;
		}
	}
	for (m_ipNextOnTour[w=0]=v=ipFirst[0]; v; v=m_ipNextOnTour[w=v])
		m_ipNextOnTour[v]=(ipFirst[v] == w)? ipSecond[v]: ipFirst[v];
} // end of CTsp::changeRecord()

void CTsp::solve()
{
	setAutoCutPattern(-1,-1);
	changeObjBound(m_iPointNum*m_dMaxDist-m_dTourLength);
	optimize();
	m_dTourLength=m_dMaxDist*m_iPointNum-getObjVal();
} // end of CTsp::solve()

void CTsp::printSolution(const char* fileName)
{
	int ct;
	char name[128];
	if (!fileName)
		getProblemName(name);
	else
		strcpy(name,fileName);
	strcat(name,".sol");
	std::ofstream fout(name);
	if (!fout.is_open())
		throw new CFileException("CTsp::printSolution",name);

	m_dTourLength=m_dMaxDist*m_iPointNum-getObjVal();

	fout << "Length " <<  m_dTourLength << std::endl;
	fout << "Tour:\n";
	fout << "1, ";
	ct=1;
	for (int v=m_ipNextOnTour[0]; v != 0; v=m_ipNextOnTour[v]) {
		fout << v+1 << ", ";
		if (!(++ct % 10))
			fout << std::endl;
	}
	fout << "1\n";
	fout.close();
} // end of CTsp::printSolution()



/*
int CTsp::Propagate(int e, int* Vars, double* d1, double* d2)
{
	int num;
	if (num=GetCycleEdges(e,m_dpX,m_ipInd) {
		Vars=m_ipInd;
		d2=(d1=m_dpA)+num;
		for (register int i=0; i < num; i++) {
			d1[i]=d2[i]=0.0; 
		}
	}
	return num;
  return 0;
}
*/
